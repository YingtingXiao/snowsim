#include "fluidsim.h"

#include "array3_utils.h"
#include "levelset_util.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "makelevelset3.h"
#include "marching_cube.h"

#define EPSILON 0.000001f

float get_weight(float dx, int i, int j, int k, Vec3f pos);
Eigen::Vector3f get_weight_gradient(float dx, int i, int j, int k, Vec3f pos);

//Parameters for computing grid forces
float compression = 0.025;   //critical compression
float stretch = 0.0075;      //critical stretch
float harden = 20;           //hardening coefficient
float young = 140000;        //Young's modulus
float poisson = 0.2;         //Poisson's ratio
float friction = 0;          //dynamic friction
float alpha = 0.95;          //weight parameter used in updating particle velocities
float mu, lambda;            //initial Lame parameters (computed from young and poisson)
float beta = 0.5;            //for semi-implicit update

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
   particle_mass = 0.05;   //initial particle mass

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

   V_s.assign(3*ni*nj*nk, 0);
   V_new.assign(3*ni*nj*nk, 0);
   K.resize(3*ni*nj*nk);

   mu = young / (2 * (1 + poisson));
   lambda = young * poisson / ((1 + poisson) * (1 - 2 * poisson));
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
         //Vec3f pos(i*dx, j*dx, k*dx);
         float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
         pos += dx/2.0f * Vec3f(a,b,c);
         //pos += dx * Vec3f(a,b,c);

         if(interpolate_value(pos/dx, liquid_phi) <= 0) {
            float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
            if(solid_phi >= 0) {
               Particle pt(pos, particle_mass);
               pt.vel = Vec3f(0, -1, 0);
               particles.push_back(pt);
            }
         }
      }
   }

   sort_particles();
   rasterize_particle_data();
   compute_particle_vol_dens();
}

//Grab liquid phi from particles
void FluidSim::compute_phi() {
   liquid_phi.assign(3*dx);
   for(unsigned int p = 0; p < particles.size(); ++p) {
      Particle pt = particles[p];
      for(int k = max(0,pt.k-1); k <= min(pt.k+1,nk-1); ++k) {
         for(int j = max(0,pt.j-1); j <= min(pt.j+1,nj-1); ++j) {
            for(int i = max(0,pt.i-1); i <= min(pt.i+1,ni-1); ++i) {
               Vec3f sample_pos((i+0.5f)*dx, (j+0.5f)*dx, (k+0.5f)*dx);
               float test_val = dist(sample_pos, pt.pos) - 0.5 * dx;
               if(test_val < liquid_phi(i,j,k))
                  liquid_phi(i,j,k) = test_val;
            }
         }
      }
   }
   
   ////extend phi slightly into solids (this is a simple, naive approach, but works reasonably well)
   //Array3f phi_temp = liquid_phi;
   //for(int k = 0; k < nk; ++k) {
   //   for(int j = 0; j < nj; ++j) {
   //      for(int i = 0; i < ni; ++i) {
   //         if(liquid_phi(i,j,k) < 0.5*dx) {
   //            float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
   //               + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
   //            if(solid_phi_val < 0)
   //               phi_temp(i,j,k) = -0.5f*dx;
   //         }
   //      }
   //   }
   //}
   //liquid_phi = phi_temp;
}

void FluidSim::marching_cube(vector<Vec3f>& position, vector<Vec3f>& normal, vector<unsigned int>& indices){
	unsigned int triangleCount(0);
	float iso_value(0.00f);// the value for surface, by default is 0.0f
	// Loop over all cells in the grid.
	// for a current cell ijk:
	for(int k = 0; k < nk; ++k) {
      for(int j = 0; j < nj; ++j) {
         for(int i = 0; i < ni; ++i) {
			{
				Vec3f v_0(i*dx,j*dx,k*dx);
				Vec3f v_1((i+1)*dx,j*dx,k*dx);
				Vec3f v_2((i+1)*dx,(j+1)*dx,k*dx);
				Vec3f v_3(i*dx,(j+1)*dx,k*dx);
				Vec3f v_4(i*dx,j*dx,(k+1)*dx);
				Vec3f v_5((i+1)*dx,j*dx,(k+1)*dx);
				Vec3f v_6((i+1)*dx,(j+1)*dx,(k+1)*dx);
				Vec3f v_7(i*dx,(j+1)*dx,(k+1)*dx);

				std::vector<Vec3f> v;
				v.push_back(v_0);v.push_back(v_1);v.push_back(v_2);v.push_back(v_3);
				v.push_back(v_4);v.push_back(v_5);v.push_back(v_6);v.push_back(v_7);
				//1. sample values of 8 corners for current cell.
				float cornerValue[8];
				cornerValue[0] = interpolate_value(v_0/dx,liquid_phi);
				cornerValue[1] = interpolate_value(v_1/dx,liquid_phi);
				cornerValue[2] = interpolate_value(v_2/dx,liquid_phi);
				cornerValue[3] = interpolate_value(v_3/dx,liquid_phi);
				cornerValue[4] = interpolate_value(v_4/dx,liquid_phi);
				cornerValue[5] = interpolate_value(v_5/dx,liquid_phi);
				cornerValue[6] = interpolate_value(v_6/dx,liquid_phi);
				cornerValue[7] = interpolate_value(v_7/dx,liquid_phi);
				unsigned int valueMask = 0;//8 bit for indicating polarity for 8 corner.

				for(int current_corner  = 0;current_corner<8;current_corner++)
				{
					float corner_value = cornerValue[current_corner];
					
					cornerValue[current_corner] = corner_value - iso_value;

					if(cornerValue[current_corner] <= 0)
						valueMask |= (1 << current_corner);// bit operation to mark corresponding corner.
				}
				//if(valueMask>0)
				//std::cout<<"valueMask:"<<valueMask<<endl;
				//2. get vertex position and normal for intersection point.
				int edgeFlag = aiCubeEdgeFlags[valueMask];// flag for indicating which edges is intersecting
				// not intersecting at all.
				if(edgeFlag == 0)
					continue;
				// 12 correspond to 12 edges. e.g. a intersection vertex on edge 6 would be stored in edgeVertex[6] and edgeNormal[6]
				Vec3f edgeVertex[12];
				Vec3f edgeNormal[12];
				for(int edge_id = 0;edge_id<12;edge_id++)
				{
					// if edge is intersecting with the surface.
					// use bit-wise operation to find out. intersection information is stored in edgeFlag.
					{
						// get the two end point ids for current edge.
						unsigned int v1 = a2iEdgeConnection[edge_id][0];
						unsigned int v2 = a2iEdgeConnection[edge_id][1];
						// fraction is the intersection point between v1 and v2.
						float fraction =  cornerValue[v1]/ (cornerValue[v1] - cornerValue[v2]);
						// edgeDirection is the unit vector of the direction of current edge in world space.
						Vec3f edgeDirection = Vec3f(a2fEdgeDirection[edge_id][0],a2fEdgeDirection[edge_id][1],a2fEdgeDirection[edge_id][2]);
						// get edge direction from a2fEdgeDirection[edge_id]

						Vec3f pos;
						unsigned int i1, j1, k1;
						// i, j, k is the index for the corner with id 0;
						// i1, j1, k1 is the index for the corner with id v1;
						
						// computer i1, j1, k1;
						pos = v[v1];
                        // compute position of corner i1, j1, k1 in world space;
						pos += fraction * dx * edgeDirection;

						edgeVertex[edge_id] = pos;
                        
                  // compute normal for current vertex.
                  // works only for level set field, since gradient of level set is the normal.
                  Vec3f norm;
						interpolate_gradient(norm,pos/dx,liquid_phi);
                  // compute normal using the gradient.
                  // take samples from the gradient.
						if(mag(norm)>0)
							normalize(norm);
                  edgeNormal[edge_id] = norm;
					}
				}

				//3. generate triangles from the vertices and normal.
				// triangle information is stored in a2iTriangleConnectionTable[valueMask]
				int * triangleTablePtr = a2iTriangleConnectionTable[valueMask];
				for(unsigned int num = 0; num < 5; ++num)
				{// up to 5 triangles
					if(*(triangleTablePtr + 3 * num) < 0)
						break;

					for(unsigned int idx = 0; idx < 3; ++idx)
					{
						// vertex id is used for extracting position and normal from edgeVertex and edgeNormal array.
						int vertex_idx = *(triangleTablePtr + 3 * num + idx);

						Vec3f p = edgeVertex[vertex_idx];

						position.push_back(p);

                  Vec3f n = edgeNormal[vertex_idx];
						normal.push_back(n);
						indices.push_back(3 * triangleCount + idx);
					}

					triangleCount++;
				}
			}
		}
	  }
	}
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

      if (debug.i < 0 || debug.i >= ni || debug.j < 0 || debug.j >= nj || debug.k < 0 || debug.k >= nk) {
         cout << "particle " << p << ":" << endl;
         cout << "pos = " << debug.pos << endl;
         cout << "index = " << debug.i << " " << debug.j << " " << debug.k << endl;
         cout << "vel = " << debug.vel << endl;
         cout << "vol = " << debug.vol << endl;
         cout << "radius = " << debug.radius << endl;
         //particles.erase(particles.begin() + p);
         //p--;
      }
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
   u.set_zero();
   v.set_zero();
   w.set_zero();
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
         //if (abs(u(i, j, k)) > EPSILON || abs(v(i, j, k)) > EPSILON || abs(w(i, j, k)) > EPSILON) {
         //   int debug = 1;
         //}
      }
   }
   //Vec3f debug(u(18, 11, 13), v(18, 11, 13), w(18, 11, 13));
   //float debug1 = mass(18, 11, 13);
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

   Particle pt = particles[0];
   cout << "particle 0 density: " << pt.dens << endl;
   cout << "cell (" << pt.i << "," << pt.j << "," << pt.k << ") density = " << mass(pt.i, pt.j, pt.k) / pow(dx,3) << endl;

   printf("Finished computing particle density and volume\n");
   printf("Particle radius = %f\n", particles[0].radius);
}

//Compute grid forces using deformation gradients
void FluidSim::compute_grid_forces(float dt) {
   for(int k = 0; k < nk; ++k) {
      #pragma omp parallel for
      for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
         Eigen::Vector3f force = Eigen::Vector3f::Zero();

         //if (i == 17 && j == 11 && k == 16) {
         //   int debug = 1;
         //}

         for (int k1=max(0,k-2); k1<=min(k+2,nk-1); ++k1) for (int j1=max(0,j-2); j1<=min(j+2,nj-1); ++j1) for (int i1=max(0,i-2); i1<=min(i+2,ni-1); ++i1) {
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
                  sum += 0.5 * dt * vel * wg.transpose();
                  //if (i == 18 && j == 11 && k == 13) {
                  //   cout << "sum = " << endl;
                  //   cout << sum << endl;
                  //}
               }

               //analytically computer derivative of psi wrt particle's elastic deformation gradient
               Eigen::Matrix3f new_def_e = sum * pt.def_e;

               //if (i == 17 && j == 11 && k == 16 && p == 23) {
               //   int debug = 1;
               //}

               //SVD and polar decomposition of new_def_e
               Eigen::JacobiSVD<Eigen::Matrix3f> svd(new_def_e, Eigen::ComputeFullU | Eigen::ComputeFullV);
               Eigen::Matrix3f matrixU = svd.matrixU();
               Eigen::Matrix3f matrixV = svd.matrixV();
               Eigen::Vector3f singulars = svd.singularValues();
               Eigen::Matrix3f sigma = Eigen::Matrix3f::Zero();
               for (int ss=0; ss<3; ++ss) {
                  sigma(ss, ss) = singulars(ss);
               }
               Eigen::Matrix3f matrixR = matrixU * matrixV.transpose();
               Eigen::Matrix3f matrixS = matrixV * sigma * matrixV.transpose();

               //compute the two functions of plastic deformation gradient
               float deter_p = pt.def_p.determinant();
               float exponent = exp(harden * (1 - deter_p));
               float mu_def_p = mu * exponent;
               float lambda_def_p = lambda * exponent;

               float deter_e = new_def_e.determinant();

               Eigen::Matrix3f temp = matrixS - Eigen::Matrix3f::Identity();   //so we don't have to recompute it
               //Eigen::Matrix3f gradient = 2 * mu_def_p * temp + lambda_def_p * temp.trace() * Eigen::Matrix3f::Identity();
               Eigen::Matrix3f gradient = 2 * mu_def_p * (new_def_e - matrixR) + lambda_def_p * (deter_e - 1) * deter_e * new_def_e.transpose().inverse();

               //for new_def_e
               gradient = sum * gradient;

               //compute grid force
               Eigen::Vector3f oldForce = force;
               force -= pt.init_vol * gradient * pt.def_e.transpose() * get_weight_gradient(dx, i, j, k, pt.pos);

               Eigen::Matrix3f trans = pt.def_e.transpose();
               Eigen::Matrix3f temp1 = gradient * trans;
               Eigen::Vector3f grad = get_weight_gradient(dx, i, j, k, pt.pos);
               Eigen::Vector3f f = temp1 * grad;
               //int debug = 1;
               //if (i == 16 && j == 13 && k == 15) {
               //   if (abs(force(1) - oldForce(1)) < 20) {
               //      cout << "vel(16, 13, 15)" << endl;
               //      cout << "particle " << p << endl;
               //      cout << "prev particle position = " << pt.pos - dt * pt.vel << endl;
               //      cout << "deter_def_p = " << deter_p << endl;
               //      cout << "force = " << endl << force << endl;
               //   }
               //}
            }
         }
         fx(i, j, k) = force(0);
         fy(i, j, k) = force(1);
         fz(i, j, k) = force(2);
         //if (i == 16 && j == 13 && k == 15) {
         //   cout << "grid (" << i << ", " << j << ", " << k << ")" << endl;
         //   cout << "vel = (" << u(i, j, k) << ", " << v(i, j, k) << ", " << w(i, j, k) << ")" << endl;
         //   cout << "force = (" << fx(i, j, k) << ", " << fy(i, j, k) << ", " << fz(i, j, k) << ")" << endl;
         //}
      }
      //printf(" Finished computing forces for layer %d\n", k);
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
void FluidSim::update_temp_velocities(float dt) {
   //cout << "vel(16, 13, 15) = (" << u(16,13,15) << ", " << v(16,13,15) << ", " << w(16,13,15) << ")" << endl;
   //cout << "f(16, 13, 15) = (" << fx(16,13,15) << ", " << fy(16,13,15) << ", " << fz(16,13,15) << ")" << endl;
   //cout << "mass(16, 13, 15) = " << mass(16, 13, 15) << endl;
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      if (mass(i, j, k) > EPSILON) {
         temp_u(i, j, k) = u(i, j, k) + fx(i, j, k) * dt / mass(i, j, k);
         temp_v(i, j, k) = v(i, j, k) + fy(i, j, k) * dt / mass(i, j, k);
         temp_w(i, j, k) = w(i, j, k) + fz(i, j, k) * dt / mass(i, j, k);
      }
   }
   //cout << "new_vel(16, 13, 15) = (" << temp_u(16,13,15) << ", " << temp_v(16,13,15) << ", " << temp_w(16,13,15) << ")" << endl;
   //printf("Finished updating temporary velocities\n");
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
               //cout << "collision " << i << " " << j << " " << k << " " << endl;
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
//CRAZY SEMI-IMPLICIT UPDATE
void FluidSim::update_grid_velocities(float dt) {
   //Construct V* and V^n+1 vector
   V_s.assign(V_s.size(), 0);
   V_new.assign(V_new.size(), 0);
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      int index = k * ni * nj + j * ni + i;
      V_s[3*index] = temp_u(i, j, k);
      V_s[3*index+1] = temp_v(i, j, k);
      V_s[3*index+2] = temp_w(i, j, k);
      K.set_element(3*index, 3*index, mass(i, j, k));
      K.set_element(3*index+1, 3*index+1, mass(i, j, k));
      K.set_element(3*index+2, 3*index+2, mass(i, j, k));
   }

   //Construct the crazy 3n*3n matrix
   vector<Particle> Ps;
   vector<Eigen::MatrixXf> Ds;
   vector<Eigen::Vector3f> FWs;

   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      int index = k * ni * nj + j * ni + i;

      Ps.clear();
      Ds.clear();
      FWs.clear();
      
      for (int k1=max(0,k-2); k1<=min(k+2,nk-1); ++k1) for (int j1=max(0,j-2); j1<=min(j+2,nj-1); ++j1) for (int i1=max(0,i-2); i1<=min(i+2,ni-1); ++i1) {
         int start, end;
         find_particle_indices(i1, j1, k1, start, end);
         
         for (int p=start; p<end; ++p) {
            Particle pt = particles[p];
            float deter_p = pt.def_p.determinant();
            float exponent = exp(harden * (1 - deter_p));
            float mu_p = mu * exponent;
            Eigen::MatrixXf D = get_second_derivative(pt.def_e, mu_p);
            Eigen::Vector3f FW = pt.def_e.transpose() * get_weight_gradient(dx, i, j, k, pt.pos);
            
            Ps.push_back(pt);
            Ds.push_back(D);
            FWs.push_back(FW);
         }
      }

      for (int kk=max(0,k-4); kk<=min(k+5,nk-1); ++kk) for (int jj=max(0,j-4); jj<=min(j+5,nj-1); ++jj) for (int ii=max(0,i-4); ii<=min(i+5,ni-1); ++ii) {
         int index1 = kk * ni * nj + jj * ni + ii;
         Eigen::Matrix3f Kij = Eigen::Matrix3f::Zero();

         for (int p=0; p<Ps.size(); ++p) {
            Particle pt = Ps[p];
            Eigen::MatrixXf D = Ds[p];
            Eigen::Vector3f FW = FWs[p];
            Eigen::Vector3f WF = (get_weight_gradient(dx, ii, jj, kk, pt.pos).transpose() * pt.def_e).transpose();

            for (int r=0; r<3; ++r) {
               for (int s=0; s<3; ++s) {
                  for (int t=0; t<3; ++t) {
                     for (int l=0; l<3; ++l) {
                        Kij(r, s) += pt.init_vol * D(3*s+r, 3*l+t) * WF(l) * FW(t);
                     }
                  }
               }
            }
         }

         for (int r=0; r<3; ++r) {
            for (int s=0; s<3; ++s) {
               K.add_to_element(index*3+r, index1*3+s, beta * dt * Kij(r, s));
            }
         }
      }
   }

   printf("Finished constructing semi-implicit update equations\n");

   //Solve for V^n+1 using Conjugate Gradient solver
   double tolerance;
   int iterations;
   solver.set_solver_parameters(1e-18, 1000);
   bool success = solver.solve(K, V_s, V_new, tolerance, iterations);
   printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
   if(!success) {
      printf("WARNING: Pressure solve failed!************************************************\n");
   }

   //Output data in V^n+1 back to temp_u, temp_v, temp_w
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      int index = k * ni * nj + j * ni + i;
      temp_u(i, j, k) = V_new[3*index];
      temp_v(i, j, k) = V_new[3*index+1];
      temp_w(i, j, k) = V_new[3*index+2];
   }

   printf("Finished updating grid velocities\n");
}

//Get second order derivative of psi wrt Fe
Eigen::MatrixXf FluidSim::get_second_derivative(Eigen::Matrix3f F, float mu_p) {
   //Compute SVD
   Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
   Eigen::Matrix3f U = svd.matrixU();
   Eigen::Matrix3f V = svd.matrixV();
   Eigen::Vector3f S = svd.singularValues();

   Eigen::MatrixXf deriv(9, 9);
   for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
         //deriv[i][j] = 2 * mu_p * (dFji + dR/dFji)
         Eigen::Matrix3f dF = Eigen::Matrix3f::Zero();
         //Compute dF/dF
         dF(j, i) = 1;

         //Compute dR/dFji = dU/dF * V^T + U * dV^T/dF
         //Compute dU/dF and dV^T/dF first
         Eigen::Matrix3f L = U.transpose() * dF * V;
         Eigen::MatrixXf C(6, 6);   //coefficients for solving for U^tilde and V^tilde
         C.setZero();
         C(0, 0) = -S(0); C(0, 3) = -S(1);
         C(1, 1) = -S(0); C(1, 4) = -S(2);
         C(2, 0) =  S(1); C(2, 3) =  S(0);
         C(3, 2) = -S(1); C(3, 5) = -S(2);
         C(4, 1) =  S(2); C(4, 4) =  S(0);
         C(5, 2) =  S(2); C(5, 5) =  S(1);
         Eigen::VectorXf D(6);   //constants for solving for U^tilde and V^tilde
         D(0) = L(1, 0); D(1) = L(2, 0); D(2) = L(0, 1);
         D(3) = L(2, 1); D(4) = L(0, 2); D(5) = L(1, 2);
         Eigen::VectorXf UV = C.inverse() * D;   //(u1, u2, u3, v1, v2, v3)

         Eigen::Matrix3f Ut = Eigen::Matrix3f::Zero();   //U^tilde
         Ut(0, 1) =  UV(0); Ut(0, 2) =  UV(1); Ut(1, 2) =  UV(2);
         Ut(1, 0) = -UV(0); Ut(2, 0) = -UV(1); Ut(2, 1) = -UV(2);
         Eigen::Matrix3f Vt = Eigen::Matrix3f::Zero();   //V^tilde
         Vt(0, 1) =  UV(3); Vt(0, 2) =  UV(4); Vt(1, 2) =  UV(5);
         Vt(1, 0) = -UV(3); Vt(2, 0) = -UV(4); Vt(2, 1) = -UV(5);

         Eigen::Matrix3f dU = U.transpose().inverse() * Ut;
         Eigen::Matrix3f dVT = Vt * V.inverse();
         Eigen::Matrix3f dR = dU * V.transpose() + U * dVT;
         Eigen::Matrix3f deriv_ij = 2 * mu_p * (dF + dR);

         for (int i1=0; i1<3; ++i) {
            for (int j1=0; j1<3; ++j) {
               deriv(i*3+i1, j*3+j1) = deriv_ij(i1, j1);
            }
         }
      }
   }
}

//Update elastic and plastic deform gradient of particles
void FluidSim::update_deform_gradient(float dt) {
   #pragma omp parallel for
   for (int p=0; p<particles.size(); ++p) {
      Particle pt = particles[p];
      Eigen::Matrix3f sum = Eigen::Matrix3f::Identity();
      for(int k=max(0,pt.k-2); k<=min(pt.k+2,nk-1); ++k) for(int j=max(0,pt.j-2); j<=min(pt.j+2,nj-1); ++j) for(int i=max(0,pt.i-2); i<=min(pt.i+2,ni-1); ++i) {
         Eigen::Vector3f v_grid;
         v_grid(0) = temp_u(i, j, k);
         v_grid(1) = temp_v(i, j, k);
         v_grid(2) = temp_w(i, j, k);
         Eigen::Vector3f w = get_weight_gradient(dx, i, j, k, pt.pos);
         sum += dt * v_grid * w.transpose();
         //if (p == 21) {
         //   cout << "(" << i << ", " << j << ", " << k << ")" << endl;
         //   cout << "grid velocity = \n" << v_grid << endl;
         //   cout << "weight gradient = \n" << w << endl; 
         //}
      }

      //if (p == 19) {
      //   cout << "particle 74" << endl;
      //   cout << "old def_e" << endl;
      //   cout << pt.def_e << endl;
      //   cout << "old def_p" << endl;
      //   cout << pt.def_p << endl;
      //   cout << sum << endl;
      //}

      //compute new total deformation gradient and particle volume
      Eigen::Matrix3f new_def_e = sum * pt.def_e;
      Eigen::Matrix3f new_def = new_def_e * pt.def_p;
      float deter = new_def.determinant();
      particles[p].vol = particles[p].init_vol * deter;
      particles[p].radius = pow(particles[p].vol * 3/4.0f / PI, 1/3.0f);
      //if (p == 4012) {
      //   cout << "particle " << p << endl;
      //   cout << "vol = " << particles[p].vol << endl;
      //   cout << "radius = " << particles[p].radius << endl;
      //   cout << "old def_e = " << pt.def_e << endl;
      //   cout << "old def_p = " << pt.def_p << endl;
      //   cout << "sum = " << sum << endl;
      //}

      //recompute new elastic and plastic deformation gradient
      Eigen::JacobiSVD<Eigen::Matrix3f> svd(new_def_e, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::Matrix3f matrixU = svd.matrixU();
      Eigen::Matrix3f matrixV = svd.matrixV();
      Eigen::Vector3f singulars = svd.singularValues();
      Eigen::Matrix3f matrixS = Eigen::Matrix3f::Zero();
      for (int s=0; s<3; ++s) {
         matrixS(s, s) = clamp(singulars(s), 1 - compression, 1 + stretch);
      }
      particles[p].def_e = matrixU * matrixS * matrixV.transpose();
      particles[p].def_p = matrixV * matrixS.inverse() * matrixU.transpose() * new_def;

      //if (p == 21) {
      //   cout << "particle " << p << endl;
      //   cout << "old def_p determinant = " << pt.def_p.determinant() << endl;
      //   cout << "def_p determinant = " << particles[p].def_p.determinant() << endl;
      //   cout << "sum = " << endl << sum << endl;
      //}

      //if (p == 74) {
      //   cout << "particle 74" << endl;
      //   cout << "new def_p" << endl;
      //   cout << particles[p].def_p << endl;
      //}
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
      //if (p == 0) {
      //   cout << "particle " << p << endl;
      //   cout << "vel = " << particles[p].vel << endl;
      //}
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
         particles[p].pos += (pt.radius - solid_phi) * normal;
      }
   }
   printf("Finished stablization\n");
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

   //   if (!first_step) {
   //      rasterize_particle_data();
   //   }
   //   compute_grid_forces(substep);
   //   ////temporary
   //   //fx.set_zero();
   //   //fy.set_zero();
   //   //fz.set_zero();
   //   apply_external_force();
   //   update_temp_velocities(substep);
   //   apply_collision_to_grid();
   //   update_deform_gradient(substep);
   //   update_particle_velocities();
   //   apply_collision_to_particles(substep);
   //   update_grid_velocities();
   //   update_particle_positions(substep);
   //   sort_particles();
   //   stablize();

   //   first_step = false;
   //}

   compute_grid_forces(dt);
   ////temporary
   //fx.set_zero();
   //fy.set_zero();
   //fz.set_zero();
   apply_external_force();
   update_temp_velocities(dt);
   apply_collision_to_grid();
   //update_grid_velocities(dt);
   update_deform_gradient(dt);
   update_particle_velocities();
   apply_collision_to_particles(dt);
   update_particle_positions(dt);
   sort_particles();
   //stablize();
   rasterize_particle_data();
   compute_phi();
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

float get_derivative_of_N(float x) {
   if (x >= 0 && x < 1) {
      return 1.5 * pow(x, 2) - 2 * x;
   } else if (x >= 1 && x < 2) {
      return -0.5 * pow(x, 2) + 2 * x - 2;
   } else if (x < 0 && x > -1) {
      return -1.5 * pow(x, 2) - 2 * x;
   } else if (x <= -1 && x > -2) {
      return 0.5 * pow(x, 2) + 2 * x + 2;
   } else {
      return 0;
   }
}

Eigen::Vector3f get_weight_gradient(float dx, int i, int j, int k, Vec3f pos) {
   float px = pos[0]/dx - i;
   float py = pos[1]/dx - j;
   float pz = pos[2]/dx - k;

   float Nx = get_weight_1d(px);
   float Ny = get_weight_1d(py);
   float Nz = get_weight_1d(pz);

   float gx = get_derivative_of_N(px);
   float gy = get_derivative_of_N(py);
   float gz = get_derivative_of_N(pz);

   Eigen::Vector3f gradient;
   gradient(0) = 1/dx * Ny * Nz * gx;
   gradient(1) = 1/dx * Nx * Nz * gy;
   gradient(2) = 1/dx * Nx * Ny * gz;
   return gradient;
}