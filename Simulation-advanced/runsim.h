// this file is distributed under
// GPL license
#ifndef ______RUNADVANCEDSIMULATION_H_____
#define ______RUNADVANCEDSIMULATION_H_____
#include <string>
#include <functional>
#include <random>
#include <math_h/vectors.h>
#include <Kinematics/particles.h>
struct particle_sim{Particle type;MathTemplates::Vector<3> P;};
typedef std::function<std::list<particle_sim>()> EventGenerator;
void Simulate(const std::string&filename,const EventGenerator gen,const size_t millions=10);
#endif

