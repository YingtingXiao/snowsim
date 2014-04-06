#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec.h"
#include "Eigen/Dense"

#define PI 3.14159265359

struct Particle {
   Vec3f pos;   //position
   float mass;  //mass
   Vec3f vel;   //velocity
   float dens;  //density
   float vol;   //volume
   float init_vol;
   float radius;
   int i;
   int j;
   int k;
   int index;
   Eigen::Matrix3f def_e;   //elastic deformation gradient
   Eigen::Matrix3f def_p;   //plastic deformation gradient

   Particle() {}
   Particle(Vec3f p, float m) : pos(p), mass(m), vel(0), dens(0), vol(0), init_vol(0), radius(0) {
      def_e = Eigen::Matrix3f::Identity();
      def_p = Eigen::Matrix3f::Identity();
   }

   void compute_index(float dx, int ni, int nj, int nk) {
      i = floor(pos[0]/dx);
      j = floor(pos[1]/dx);
      k = floor(pos[2]/dx);
      index = i + j * ni + k * ni * nj;
   }
};

#endif