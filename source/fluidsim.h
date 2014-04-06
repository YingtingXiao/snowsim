#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array3.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "particle.h"

#include <vector>

using namespace std;

class FluidSim {

public:
	void initialize(float width, int ni_, int nj_, int nk_);
	void set_boundary(float (*phi)(const Vec3f&));
	void set_boundary(string filename);
	void set_liquid(float (*phi)(const Vec3f&));
	void set_liquid(string filename);

	void advance(float dt, bool first_step);

	//Grid dimensions
	int ni,nj,nk;
	float dx;

	//Fluid velocity
	Array3f u, v, w;
	Array3f temp_u, temp_v, temp_w;

   //Grid force
   Array3f fx, fy, fz;

	//Static geometry representation
	Array3f nodal_solid_phi;

   Array3i particle_indices;   //the indices of the first particles in each grid
   std::vector<Particle> particles;
   float particle_mass;   //initial particle mass

	Array3f liquid_phi;

   //Grid mass and rasterized from particles
   Array3f mass;

	//Solver data
	PCGSolver<double> solver;
	SparseMatrixd matrix;
	std::vector<double> rhs;
	std::vector<double> pressure;

	Vec3f get_velocity(const Vec3f& position);

private:
	void initialize_phi(string filename, Array3f& phi, bool solid);
   bool load_levelset(string filename, Array3f& phi);
   void sort_particles();
   void find_particle_indices(int i, int j, int k, int& start, int& end);

   void rasterize_particle_data();
   void compute_particle_vol_dens();
   void compute_grid_forces(float dt);
   void apply_external_force();
   float compute_psi(Eigen::Matrix3f def_e, Eigen::Matrix3f def_p, Eigen::Matrix3f sum);
   void update_temp_velocities(float dt);
   void update_grid_velocities();
   void apply_collision_to_grid();
   void update_deform_gradient(float dt);
   void update_particle_velocities();
   void apply_collision_to_particles(float dt);
   void update_particle_positions(float dt);

	Vec3f trace_rk2(const Vec3f& position, float dt);

	float cfl();

	void advect_particles(float dt);
	void add_force(float dt);
};


#endif