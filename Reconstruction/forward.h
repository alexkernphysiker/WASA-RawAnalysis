// this file is distributed under 
// GPL license
#ifndef ______RECONSTRUCTION____FORWARD_H_______________
#define ______RECONSTRUCTION____FORWARD_H_______________
#include <Kinematics/reactions.h>
#include <analysis.h>
#include <trackprocessing.h>
shared_ptr<TrackAnalyse::AbstractChain> ForwardHe3Reconstruction(const Analysis&data,particle_kinematics&kin_rec);
#endif 
