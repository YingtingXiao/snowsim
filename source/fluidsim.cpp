#include "fluidsim.h"

#include "array3_utils.h"
#include "levelset_util.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "makelevelset3.h"

#define EPSILON 0.000001f

float get_weight(float dx, int i, int j, int k, Vec3f pos);
Eigen::Vector3f get_weight_gradient(float dx, int i, int j, int k, Vec3f pos);

//Parameters for computing grid forces
float compression = 0.0025;   //critical compression
float stretch = 0.00075;      //critical stretch
float harden = 10;            //hardening coefficient
float young = 140000;         //Young's modulus
float poisson = 0.2;          //Poisson's ratio
float friction = 0.2;           //dynamic friction
float alpha = 0.95;           //weight parameter used in updating particle velocities

void FluidSim::initialize(float width, int ni_, int nj_, int nk_) {
   ni = ni_;
   nj = nj_;
   nk = nk_;
   dx = width / (float)ni;
   u.resize(ni,nj,nk); temp_u.resize(ni,nj,nk);
   v.resize(ni,nj,nk); temp_v.resize(ni,nj,nk);
   w.resize(ni,nj,nk); temp_w.resize(ni,nj,nk);
 
   particle_indices.resize(ni, nj, nk);
   particle_indices.set_zero();
   particle_mass = 1;   //initial particle mass

   u.set_zero();
   v.set_zero();
   w.set_zero();
   temp_u.set_zero();
   temp_v.set_zero();
   temp_w.set_zero();
   nodal_solid_phi.resize(ni,nj,nk);
   liquid_phi.resize(ni,nj,nk);
   mass.resize(ni,nj,nk);
   mass.set_zero();

   fx.resize(ni, nj, nk);
   fx.set_zero();
   fy.resize(ni, nj, nk);
   fy.set_zero();
   fz.resize(ni, nj, nk);
   fz.set_zero();
}

void FluidSim::initialize_phi(string filename, Array3f& phi, bool solid) {
   if(filename.size() < 5 || filename.substr(filename.size()-4) != string(".obj")) {
      cout << filename << endl;
      cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
      return;
   }

   int size = ni;
   string outname = filename.substr(0, filename.size()-4) + "_" + to_string(size) + string(".sdf");
   if (!load_levelset(outname, phi)) {
      ifstream infile(filename);
      if(!infile) {
         std::cerr << "obj file doesn't exist.\n";
         return;
      }

      Vec3f min_box(0.0f, 0.0f, 0.0f), 
         max_box(1.0f, 1.0f, 1.0f);

      string line;
      vector<Vec3f> vertList;
      vector<Vec3ui> faceList;
      while(!infile.eof()) {
         std::getline(infile, line);
         if(line.substr(0,1) == string("v") && line.substr(1,2) != string("t") && line.substr(1,2) != string("n")) {
            stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vertList.push_back(point);
         }
         else if(line.substr(0,1) == string("f")) {
            stringstream data(line);
            int v0,v1,v2;
            string test = strtok((char*)line.c_str(), "\  ");
            test = strtok(NULL, "\  ");
            v0 = atoi(test.c_str());
            test = strtok(NULL, "\  ");
            v1 = atoi(test.c_str());
            test = strtok(NULL, "\  ");
            v2 = atoi(test.c_str());
            faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
         }
      }
      infile.close();

      float maxDim = max(max(max_box[0]-min_box[0], max_box[1]-min_box[1]), max_box[2]-min_box[2]);
      float dx = maxDim/size;

      cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << size << "." << endl;

      cout << "Computing signed distance field.\n";
      make_level_set3(faceList, vertList, min_box, dx, size, size, size, phi);

      //Flip the sign for solid phi
      if (solid) {
         for (unsigned int i=0; i<phi.a.size(); ++i) {
            phi.a[i] = -phi.a[i];
         }
      }

      //Very hackily strip off file suffix.
      string outname = filename.substr(0, filename.size()-4) + "_" + to_string(size) + std::string(".sdf");
      cout << "Writing results to: " << outname << "\n";
      ofstream outfile( outname.c_str());
      outfile << size << " " << size << " " << size << endl;
      outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;
      outfile << dx << std::endl;
      for(unsigned int i = 0; i < phi.a.size(); ++i) {
         outfile << phi.a[i] << std::endl;
      }
      outfile.close();
      cout << "Processing complete.\n";
   }
}

//Load level set data from sdf file
bool FluidSim::load_levelset(string filename, Array3f& phi) {
   ifstream infile(filename);
   if(!infile) {
      return false;
   }

   cout << "Loading level set data from " << filename << endl;
   string line;
   getline(infile, line);
   getline(infile, line);
   getline(infile, line);
   int size = ni;
   for(int k = 0; k < size; ++k) for(int j = 0; j < size; ++j) for(int i = 0; i < size; ++i) {
      getline(infile, line);
      phi(i,j,k) = atof(line.c_str());
   }
   infile.close();
   cout << "Finished loading" << endl;
   return true;
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float (*phi)(const Vec3f&)) {
   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec3f pos(i*dx,j*dx,k*dx);
      nodal_solid_phi(i,j,k) = phi(pos);
   }

}

//Initialize solid boundary phi from obj file
void FluidSim::set_boundary(string filename) {
   initialize_phi(filename, nodal_solid_phi, true);
}

void FluidSim::set_liquid(float (*phi)(const Vec3f&)) {
   //initialize particles
   int seed = 0;
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vec3f pos(i*dx,j*dx,k*dx);
      float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
      pos += dx * Vec3f(a,b,c);

      if(phi(pos) <= 0) {
         float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
         if(solid_phi >= 0)
            particles.push_back(Particle(pos, particle_mass));
      }
   }
   
   sort_particles();
   rasterize_particle_data();
   compute_particle_vol_dens();
}

//Initialize liquid phi from obj file
void FluidSim::set_liquid(string filename) {
   initialize_phi(filename, liquid_phi, false);

   //initialize particles
   int seed = 0;
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      for (int r=0; r<2; ++r) for (int s=0; s<2; ++s) for (int t=0; t<2; ++t) {
         //add a particle in each of the 8 subgrid cells
         Vec3f pos((i+r/2.0f)*dx,(j+s/2.0f)*dx,(k+t/2.0f)*dx);
         float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
         pos += dx/2.0f * Vec3f(a,b,c);

         if(interpolate_value(pos/dx, liquid_phi) <= 0) {
            float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
            if(solid_phi >= 0) {

               particles.push_back(Particle(pos, particle_mass));
            }
         }
      }
   }

   sort_particles();
   rasterize_particle_data();
   compute_particle_vol_dens();
}

//Compare function for sorting
bool compare_particles(Particle a, Particle b) {
   return a.index < b.index;
}

//Sort particles on their indices in the grid
void FluidSim::sort_particles() {
   particle_indices.set_zero();
   for (int p=0; p<particles.size(); ++p) {
      particles[p].compute_index(dx, ni, nj, nk);
      Particle debug = particles[p];
      //if (p == 0) {
      //   cout << "particle 0:" << endl;
      //   cout << "pos = " << debug.pos << endl;
      //   cout << debug.i << " " << debug.j << " " << debug.k << endl;
      //   cout << "radius = " << debug.radius << endl;
      //}

      //if (debug.i < 0 || debug.i >= ni || debug.j < 0 || debug.j >= nj || debug.k < 0 || debug.k >= nk) {
      //   cout << "particle " << p << ":" << endl;
      //   cout << "pos = " << debug.pos << endl;
      //   cout << debug.i << " " << debug.j << " " << debug.k << endl;
      //   cout << "vel = " << debug.vel << endl;
      //   cout << "vol = " << debug.vol << endl;
      //   cout << "radius = " << debug.radius << endl;
      //   particles.erase(particles.begin() + p);
      //   p--;
      //}
      //else {
         particle_indices(particles[p].i, particles[p].j, particles[p].k)++;
      //}
   }
   Particle debug = particles[0];
   sort(particles.begin(), particles.end(), compare_particles);
   for (int i=1; i<particle_indices.size(); ++i) {
      particle_indices.a[i] += particle_indices.a[i-1];
   }
   debug = particles[0];
   printf("Finished sorting particles\n");
}

void FluidSim::find_particle_indices(int i, int j, int k, int& start, int& end) {
   int ind = i + j * ni + k * ni * nj;
   if (ind == 0) {
      start = 0;
   } else {
      start = particle_indices.a[ind-1];
   }
   end = particle_indices.a[ind];
}

//Rasterize particles' mass and velocity data to grid
void FluidSim::rasterize_particle_data() {
   mass.set_zero();
   for(int k = 0; k < nk; ++k) {
      #pragma omp parallel for
      for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
         for (int kk=max(0,k-2); kk<=min(k+2,nk-1); ++kk) for (int jj=max(0,j-2); jj<=min(j+2,nj-1); ++jj) for (int ii=max(0,i-2); ii<=min(i+2,ni-1); ++ii) {
            int start, end;
            find_particle_indices(ii, jj, kk, start, end);
            for (int p=start; p<end; ++p) {
               float weight = get_weight(dx, i, j, k, particles[p].pos);
               mass(i, j, k) += weight * particles[p].mass;
            }
         }
         if (mass(i, j, k) > EPSILON) {
            for (int kk=max(0,k-2); kk<=min(k+2,nk-1); ++kk) for (int jj=max(0,j-2); jj<=min(j+2,nj-1); ++jj) for (int ii=max(0,i-2); ii<=min(i+2,ni-1); ++ii) {
               int start, end;
               find_particle_indices(ii, jj, kk, start, end);
               for (int p=start; p<end; ++p) {
                  float weight = get_weight(dx, i, j, k, particles[p].pos);
                  u(i, j, k) += particles[p].vel[0] * particles[p].mass * weight / mass(i, j, k);
                  v(i, j, k) += particles[p].vel[1] * particles[p].mass * weight / mass(i, j, k);
                  w(i, j, k) += particles[p].vel[2] * particles[p].mass * weight / mass(i, j, k);
               }
            }
         }
         if (abs(u(i, j, k)) > EPSILON || abs(v(i, j, k)) > EPSILON || abs(w(i, j, k)) > EPSILON) {
            int debug = 1;
         }
      }
   }
   Vec3f debug(u(18, 11, 13), v(18, 11, 13), w(18, 11, 13));
   float debug1 = mass(18, 11, 13);
   printf("Finished rasterizing particle data to the grid\n");
}

//Computer particles' volumes and densities from grid
void FluidSim::compute_particle_vol_dens() {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      Particle pt = particles[p];
      for(int k=max(0,pt.k-2); k<=min(pt.k+2,nk-1); ++k) for(int j=max(0,pt.j-2); j<=min(pt.j+2,nj-1); ++j) for(int i=max(0,pt.i-2); i<=min(pt.i+2,ni-1); ++i) {
         float weight = get_weight(dx, i, j, k, pt.pos);
         particles[p].dens += mass(i, j, k) * weight / pow(dx, 3);
      }
      particles[p].init_vol = particles[p].mass / particles[p].dens;
      particles[p].vol = particles[p].init_vol;
      particles[p].radius = pow(particles[p].vol * 3/4.0f / PI, 1/3.0f);
      if (p == 0) {
         Particle debug = particles[p];
         int debug1 = 1;
      }
   }
   printf("Finished computing particle density and volume\n");
   printf("Particle radius = %f\n", particles[0].radius);
}

//Compute grid forces using deformation gradients
void FluidSim::compute_grid_forces(float dt) {
   for(int k = 0; k < nk; ++k) {
      #pragma omp parallel for
      for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
         Eigen::Vector3f force = Eigen::Vector3f::Zero();

         for (int k1=max(0,k-3); k1<=min(k+2,nk-1); ++k1) for (int j1=max(0,j-3); j1<=min(j+2,nj-1); ++j1) for (int i1=max(0,i-3); i1<=min(i+2,ni-1); ++i1) {
            int start, end;
            find_particle_indices(i1, j1, k1, start, end);
 
            for (int p=start; p<end; ++p) {
               Particle pt = particles[p];
               Eigen::Matrix3f sum = Eigen::Matrix3f::Identity();
               for(int kk=max(0,pt.k-2); kk<=min(pt.k+3,nk-1); ++kk) for(int jj=max(0,pt.j-2); jj<=min(pt.j+3,nj-1); ++jj) for(int ii=max(0,pt.i-2); ii<=min(pt.i+3,ni-1); ++ii) {
                  Eigen::Vector3f vel;
                  vel(0) = u(ii,jj,kk);
                  vel(1) = v(ii,jj,kk);
                  vel(2) = w(ii,jj,kk);
                  Eigen::Vector3f wg = get_weight_gradient(dx, ii, jj, kk, pt.pos);
                  //if (i == 18 && j == 11 && k == 13) {
                  //   cout << "particle index = " << pt.i << " " << pt.j << " " << pt.k << endl;
                  //   cout << "cell index = " << ii << " " << jj << " " << kk << endl;
                  //   cout << "vel = " << vel << endl;
                  //   cout << "wg = " << endl;
                  //   cout << wg;
                  //   cout << endl;
                  //}
                  sum += dt * vel * wg.transpose();
                  //if (i == 18 && j == 11 && k == 13) {
                  //   cout << "sum = " << endl;
                  //   cout << sum << endl;
                  //}
               }

               //compute derivative of psi wrt particle's elastic deformation gradient
               float psi = compute_psi(pt.def_e, pt.def_p, sum);
               //if (i == 18 && j == 11 && k == 13) {
               //   cout << "def_e = " << endl;
               //   cout << pt.def_e << endl;
               //   cout << "def_p = " << endl;
               //   cout << pt.def_p << endl;
               //   cout << "sum = " << endl;
               //   cout << sum << endl;
               //   cout << "psi = " << psi;
               //   cout << endl;
               //}
               Eigen::Matrix3f gradient = Eigen::Matrix3f::Zero();
               float step = dx/2;
               for (int ii=0; ii<3; ++ii) for (int jj=0; jj<3; ++jj) {
                  Eigen::Matrix3f def_e = pt.def_e;
                  def_e(ii, jj) += step;
                  float new_psi = compute_psi(def_e, pt.def_p, sum);
                  gradient(jj, ii) = (new_psi - psi) / step;
               }

               //compute grid force
               force -= pt.init_vol * gradient * get_weight_gradient(dx, i, j, k, pt.pos);
               //if (i == 18 && j == 11 && k == 13) {
               //   cout << pt.init_vol << endl;
               //   cout << gradient << endl;
               //   cout << get_weight_gradient(dx, i, j, k, pt.pos);
               //   cout << endl;
               //}
            }
         }
         fx(i, j, k) = force(0);
         fy(i, j, k) = force(1);
         fz(i, j, k) = force(2);
         if (i == 18 && j == 11 && k == 13) {
            cout << fx(i, j, k) << " " << fy(i, j, k) << " " << fz(i, j, k) << endl;
            cout << mass(i, j, k) << endl;
         }
      }
      printf(" Finished computing forces for layer %d\n", k);
   }
   printf("Finished computing grid forces\n");
}

//Helper function for computing the derivative of psi wrt elastic deformation gradient in compute_grid_forces
float FluidSim::compute_psi(Eigen::Matrix3f def_e, Eigen::Matrix3f def_p, Eigen::Matrix3f sum) {
   Eigen::Matrix3f new_def_e = sum * def_e;

   //SVD and polar decomposition of new_def_e
   Eigen::JacobiSVD<Eigen::Matrix3f> svd(new_def_e, Eigen::ComputeFullU | Eigen::ComputeFullV);
   Eigen::Matrix3f matrixU = svd.matrixU();
   Eigen::Matrix3f matrixV = svd.matrixV();
   Eigen::Matrix3f matrixR = matrixU * matrixV.transpose();

   float squared_norm = (new_def_e - matrixR).squaredNorm();
   float deter_e = new_def_e.determinant();
   float deter_p = def_p.determinant();

   //compute the two functions of plastic deformation gradient
   float mu = young / (2 * (1 + poisson));
   float lambda = young * poisson / ((1 + poisson) * (1 - 2 * poisson));
   float exponent = exp(harden * (1 - deter_p));
   float mu_def_p = mu * exponent;
   float lambda_def_p = lambda * exponent;

   //compute psi
   float psi = mu_def_p * squared_norm + lambda_def_p/2 * pow(deter_e-1, 2);
   return psi;
}

void FluidSim::apply_external_force() {
   //gravity
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      if (mass(i, j, k) > EPSILON) {
         fy(i, j, k) -= 9.81 * mass(i, j, k);
      }
   }
   printf("Finished applying external force\n");
}

//Update temporary velocities
//TODO: use semi-implicit integration
void FluidSim::update_temp_velocities(float dt) {
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      if (mass(i, j, k) > EPSILON) {
         temp_u(i, j, k) = u(i, j, k) + fx(i, j, k) * dt / mass(i, j, k);
         temp_v(i, j, k) = v(i, j, k) + fy(i, j, k) * dt / mass(i, j, k);
         temp_w(i, j, k) = w(i, j, k) + fz(i, j, k) * dt / mass(i, j, k);
      }
   }
   printf("Finished updating temporary velocities\n");
}

//Apply collision to temporary grid velocities
void FluidSim::apply_collision_to_grid() {
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      if (abs(temp_u(i, j, k)) > EPSILON || abs(temp_v(i, j, k)) > EPSILON || abs(temp_w(i, j, k)) > EPSILON) {   //cell has particles
         if (nodal_solid_phi(i, j, k) < 0) {   //collision detected
            Vec3f vel(temp_u(i, j, k), temp_v(i, j, k), temp_w(i, j, k));
            Vec3f normal;
            interpolate_gradient(normal, Vec3f(i, j, k), nodal_solid_phi);
            normalize(normal);
            float vn = dot(vel, normal);
            if (vn < 0) {   //don't apply collision if the bodies are separating
               cout << "collision " << i << " " << j << " " << k << " " << endl;
               Vec3f vt = vel - vn * normal;
               vel = vt + friction * vn * normalized(vt);
               temp_u(i, j, k) = vel[0];
               temp_v(i, j, k) = vel[1];
               temp_w(i, j, k) = vel[2];
            }
         }
      }
   }
   printf("Finished applying collision to grid\n");
}

//Update grid velocities using explicit integration
//TODO: change to semi-implicit integration
void FluidSim::update_grid_velocities() {
   u = temp_u;
   v = temp_v;
   w = temp_w;
   printf("Finished updating grid velocities\n");
}

//Update elastic and plastic deform gradient of particles
void FluidSim::update_deform_gradient(float dt) {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      Particle pt = particles[p];
      Eigen::Matrix3f sum = Eigen::Matrix3f::Identity();
      for(int k=max(0,pt.k-2); k<=min(pt.k+3,nk-1); ++k) for(int j=max(0,pt.j-2); j<=min(pt.j+3,nj-1); ++j) for(int i=max(0,pt.i-2); i<=min(pt.i+3,ni-1); ++i) {
         Eigen::Vector3f v_grid;
         v_grid(0) = temp_u(i, j, k);
         v_grid(1) = temp_v(i, j, k);
         v_grid(2) = temp_w(i, j, k);
         sum += dt * v_grid * get_weight_gradient(dx, i, j, k, pt.pos).transpose();
      }

      //compute new total deformation gradient and particle volume
      Eigen::Matrix3f new_def_e = sum * pt.def_e;
      Eigen::Matrix3f new_def = new_def_e * pt.def_p;
      float deter = max(new_def.determinant(), 0.0f);
      particles[p].vol = particles[p].init_vol * deter;
      particles[p].radius = pow(particles[p].vol * 3/4.0f / PI, 1/3.0f);
      if (p == 4012) {
         cout << "particle " << p << endl;
         cout << "vol = " << particles[p].vol << endl;
         cout << "radius = " << particles[p].radius << endl;
         cout << "old def_e = " << pt.def_e << endl;
         cout << "old def_p = " << pt.def_p << endl;
         cout << "sum = " << sum << endl;
      }

      //recompute new elastic and plastic deformation gradient
      Eigen::JacobiSVD<Eigen::Matrix3f> svd(new_def_e, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::Matrix3f matrixU = svd.matrixU();
      Eigen::Matrix3f matrixV = svd.matrixV();
      Eigen::Vector3f singulars = svd.singularValues();
      Eigen::Matrix3f matrixS = Eigen::Matrix3f::Zero();
      for (int i=0; i<3; ++i) {
         matrixS(i, i) = clamp(singulars(i), 1 - compression, 1 + stretch);
      }
      particles[p].def_e = matrixU * matrixS * matrixV.transpose();
      particles[p].def_p = matrixV * matrixS.inverse() * matrixU.transpose() * new_def;
   }
   printf("Finished updating deformation gradient\n");
}

void FluidSim::update_particle_velocities() {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      Particle pt = particles[p];
      Vec3f v_pic(0);
      Vec3f v_flip = pt.vel;
      for(int k=max(0,pt.k-2); k<=min(pt.k+2,nk-1); ++k) for(int j=max(0,pt.j-2); j<=min(pt.j+2,nj-1); ++j) for(int i=max(0,pt.i-2); i<=min(pt.i+2,ni-1); ++i) {
         float weight = get_weight(dx, i, j, k, pt.pos);
         v_pic[0] += temp_u(i, j, k) * weight;
         v_pic[1] += temp_v(i, j, k) * weight;
         v_pic[2] += temp_w(i, j, k) * weight;
         v_flip[0] += (temp_u(i, j, k) - u(i, j, k)) * weight;
         v_flip[1] += (temp_v(i, j, k) - v(i, j, k)) * weight;
         v_flip[2] += (temp_w(i, j, k) - w(i, j, k)) * weight;

         //if (p == 2555) {
         //   cout << "particle " << p << endl;
         //   cout << "pic = " << v_pic << endl;
         //   cout << "flip = " << v_flip << endl;
         //   cout << i << " " << j << " " << k << endl;
         //}
      } 
      particles[p].vel = (1 - alpha) * v_pic + alpha * v_flip;
      if (p == 0) {
         cout << "particle " << p << endl;
         cout << "vel = " << particles[p].vel << endl;
      }
   }
   printf("Finished updating particle velocities\n");
}

//Apply collision to temporary grid velocities
void FluidSim::apply_collision_to_particles(float dt) {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      //use predicted position instead of current position
      Vec3f pos = particles[p].pos + particles[p].vel * dt;
      float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
      //if (p == 4041) {
      //   cout << "particle " << p << endl;
      //   cout << "predicted pos = " << pos << endl;
      //   cout << "vel = " << particles[p].vel << endl;
      //   cout << "radius = " << particles[p].radius << endl;
      //   Particle debug = particles[p];
      //   cout << "forcey = " << fy(debug.i, debug.j, debug.k) << endl;
      //   cout << "solid phi = " << solid_phi << endl;
      //}
      if (solid_phi < particles[p].radius) {   //collision detected
         //cout << "collision " << p << endl;
         Vec3f vel = particles[p].vel;
         Vec3f normal;
         interpolate_gradient(normal, pos/dx, nodal_solid_phi);
         normalize(normal);
         float vn = dot(vel, normal);
         if (vn < 0) {   //don't apply collision if the bodies are separating
            Vec3f vt = vel - vn * normal;
            vel = vt + friction * vn * normalized(vt);
            //if (p == 2555) {
            //   Particle debug = particles[p];
            //   cout << "particle " << p << endl;
            //   cout << "old vel = " << debug.vel << endl;
            //   cout << "new vel = " << vel << endl;
            //}
            particles[p].vel = vel;
         }

         pos = particles[p].pos + particles[p].vel * dt;
         solid_phi = interpolate_value(pos/dx, nodal_solid_phi);

         //if (p == 4041) {
         //   cout << "particle " << p << endl;
         //   cout << "normal = " << normal << endl;
         //   cout << "vn = " << vn << endl;
         //   cout << "vt = " << vel - vn * normal << endl;
         //   cout << "vel = " << particles[p].vel << endl;
         //   cout << "new pos = " << particles[p].pos + particles[p].vel * dt;
         //   cout << endl;
         //}
      }

      if (solid_phi < particles[p].radius) {
         particles[p].vel = 0;
      }
   }
   printf("Finished applying collision to particles\n");
}

//Update particle positions
void FluidSim::update_particle_positions(float dt) {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      //if (p == 2555) {
      //   Particle debug = particles[p];
      //   cout << "particle " << p << endl;
      //   cout << "pos = " << debug.pos << endl;
      //   cout << "vel = " << debug.vel << endl;
      //}
      particles[p].pos += dt * particles[p].vel;
      //if (p == 2555) {
      //   Particle debug = particles[p];
      //   cout << "pos = " << debug.pos << endl;
      //}
   }
   printf("Finished updating particle positions\n");
}

//Push particles that are inside the solid shape out
void FluidSim::stablize() {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      Particle pt = particles[p];
      float solid_phi = interpolate_value(pt.pos/dx, nodal_solid_phi);
      if (solid_phi < pt.radius) {
         Vec3f normal;
         interpolate_gradient(normal, pt.pos/dx, nodal_solid_phi);
         normalize(normal);
         particles[p].pos += (pt.radius - solid_phi);
      }
   }
}

//The main fluid simulation step
void FluidSim::advance(float dt, bool first_step) {
   //float t = 0;

   //while(t < dt) {
   //   float substep = cfl();   
   //   if(t + substep > dt)
   //      substep = dt - t;
   //   printf("Taking substep of size %f (to %0.3f%% of the frame)\n", substep, 100 * (t+substep)/dt);
   //   t+=substep;
   //}

   if (!first_step) {
      rasterize_particle_data();
   }
   compute_grid_forces(dt);
   ////temporary
   //fx.set_zero();
   //fy.set_zero();
   //fz.set_zero();
   apply_external_force();
   update_temp_velocities(dt);
   apply_collision_to_grid();
   update_deform_gradient(dt);
   update_particle_velocities();
   apply_collision_to_particles(dt);
   update_grid_velocities();
   update_particle_positions(dt);
   sort_particles();
   stablize();
}


float FluidSim::cfl() {

   float maxvel = 0;
   for(unsigned int i = 0; i < u.a.size(); ++i)
      maxvel = max(maxvel, fabs(u.a[i]));
   for(unsigned int i = 0; i < v.a.size(); ++i)
      maxvel = max(maxvel, fabs(v.a[i]));
   for(unsigned int i = 0; i < w.a.size(); ++i)
      maxvel = max(maxvel, fabs(w.a[i]));

   return dx / maxvel;
}

void FluidSim::advect_particles(float dt) { 
   for(unsigned int p = 0; p < particles.size(); ++p) {
      particles[p].pos = trace_rk2(particles[p].pos, dt);

      //check boundaries and project exterior particles back in
      float phi_val = interpolate_value(particles[p].pos/dx, nodal_solid_phi); 
      if(phi_val < 0) {
         Vec3f grad;
         interpolate_gradient(grad, particles[p].pos/dx, nodal_solid_phi);
         if(mag(grad) > 0)
            normalize(grad);
         particles[p].pos -= phi_val * grad;
      }
   }
}

//Apply RK2 to advect a point in the domain.
Vec3f FluidSim::trace_rk2(const Vec3f& position, float dt) {
   Vec3f input = position;
   Vec3f velocity = get_velocity(input);
   velocity = get_velocity(input + 0.5f*dt*velocity);
   input += dt*velocity;
   return input;
}

//Interpolate velocity from the MAC grid.
Vec3f FluidSim::get_velocity(const Vec3f& position) {

   //Interpolate the velocity from the u and v grids
   float u_value = interpolate_value(position / dx - Vec3f(0, 0.5f, 0.5f), u);
   float v_value = interpolate_value(position / dx - Vec3f(0.5f, 0, 0.5f), v);
   float w_value = interpolate_value(position / dx - Vec3f(0.5f, 0.5f, 0), w);

   return Vec3f(u_value, v_value, w_value);
}

//Helper functions for computing interpolating weights
float get_weight_1d(float x) {
   if (abs(x) < 1) {
      return 0.5 * abs(pow(x, 3)) - pow(x, 2) + 2.0/3.0;
   } else if (abs(x) < 2) {
      return -abs(pow(x, 3))/6 + pow(x, 2) - 2 * abs(x) + 4.0/3.0;
   } else {
      return 0;
   }
}

float get_weight(float dx, int i, int j, int k, Vec3f pos) {
   float Nx = get_weight_1d(pos[0]/dx - i);
   float Ny = get_weight_1d(pos[1]/dx - j);
   float Nz = get_weight_1d(pos[2]/dx - k);
   return Nx * Ny * Nz;
}

Eigen::Vector3f get_weight_gradient(float dx, int i, int j, int k, Vec3f pos) {
   float step = dx/2;
   float weight = get_weight(dx, i, j, k, pos);
   float gx = get_weight(dx, i, j, k, pos+Vec3f(step,0,0)) - weight;
   float gy = get_weight(dx, i, j, k, pos+Vec3f(0,step,0)) - weight;
   float gz = get_weight(dx, i, j, k, pos+Vec3f(0,0,step)) - weight;
   Eigen::Vector3f gradient;
   gradient(0) = gx/step;
   gradient(1) = gy/step;
   gradient(2) = gz/step;
   return gradient;
}