#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec.h"

struct Particle {
   Vec3f pos;   //position
   Vec3f vel;   //velocity
   float dens;  //density
   float vol;   //volume

   Particle() {}
   Particle(Vec3f p) : pos(p) {}
};

#endif