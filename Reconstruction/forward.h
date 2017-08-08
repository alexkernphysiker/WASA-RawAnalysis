// this file is distributed under 
// GPL license
#ifndef ______RECONSTRUCTION____FORWARD_H_______________
#define ______RECONSTRUCTION____FORWARD_H_______________
#include <Kinematics/reactions.h>
#include <analysis.h>
#include <trackprocessing.h>
particle_kine&___test_mode___();
shared_ptr<TrackAnalyse::AbstractChain> ForwardHe3Reconstruction(const string&prefix,const Analysis&data,particle_kine&kin_rec=___test_mode___());
shared_ptr<TrackAnalyse::AbstractChain> ForwardPReconstruction(const Analysis&data,particle_kine&kin_rec=___test_mode___());
shared_ptr<TrackAnalyse::AbstractChain> ForwardDReconstruction(const Analysis&data,particle_kine&kin_rec=___test_mode___());
#endif 
