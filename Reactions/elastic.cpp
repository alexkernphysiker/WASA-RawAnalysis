// this file is distributed under
// GPL license
#include <Wasa.hh>
#include <TCutG.h>
#include <math_h/functions.h>
#include <math_h/tabledata.h>
#include <math_h/error.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include <trackprocessing.h>
#include <detectors.h>
#include <reconstruction.h>
#include <data.h>
#include <Reconstruction/forward.h>
#include "elastic.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

void ReconstructP(Analysis&res){
    const Axis Ep([&res]()->double{
        for(const auto&P:res.Vertex(0))
            if(P.particle==Particle::p())return P.E;
        return INFINITY;
    },0.0,2.5,250);
    const Axis Tp([&res]()->double{
        for(const Analysis::Kinematic&P:res.Vertex(0))
            if(P.particle==Particle::p())return P.Th*180.0/PI();
        return INFINITY;
    },0.0,20.0,200);
    const Axis Ed([&res]()->double{
        for(const auto&P:res.Vertex(0))
            if(P.particle==Particle::d())return P.E;
        return INFINITY;
    },0.0,2.5,250);
    const Axis Td([&res]()->double{
        for(const Analysis::Kinematic&P:res.Vertex(0))
            if(P.particle==Particle::d())return P.Th*180.0/PI();
        return INFINITY;
    },0.0,20.0,200);
    res.Trigger(0).pre()
	<<make_shared<Hist2D>("D","vertex-E-vs-theta",Ed,Td)
	<<make_shared<Hist2D>("P","vertex-E-vs-theta",Ep,Tp)
    ;

    res.Trigger(trigger_he3_forward.number).pre()
    ;
    static particle_kinematics p,d;
    res.Trigger(trigger_he3_forward.number).per_track()
	<<ForwardDReconstruction(res,d)
    ;
}
