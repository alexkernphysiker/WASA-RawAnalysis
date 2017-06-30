// this file is distributed under
// GPL license
#ifndef ______RUNADVANCEDSIMULATION_H_____
#define ______RUNADVANCEDSIMULATION_H_____
#include <string>
#include <functional>
#include <Kinematics/particles.h>
typedef std::function<std::list<particle_kinematics>(std::mt19937&)> EventGenerator;
void Simulate(const std::string&filename,const EventGenerator gen);
#endif

