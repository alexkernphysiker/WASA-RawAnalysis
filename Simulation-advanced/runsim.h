// this file is distributed under
// GPL license
#ifndef ______RUNADVANCEDSIMULATION_H_____
#define ______RUNADVANCEDSIMULATION_H_____
#include <string>
#include <functional>
#include <random>
#include <math_h/vectors.h>
#include <Kinematics/particles.h>
struct particle_sim{Particle type;MathTemplates::Vector3<double> P;};
typedef std::function<std::list<particle_sim>(std::mt19937&)> EventGenerator;
void Simulate(const std::string&filename,const EventGenerator gen);
#endif

