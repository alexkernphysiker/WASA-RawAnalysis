// this file is distributed under 
// GPL license
#include <Wasa.hh>
#include <TCutG.h>
#include <math_h/functions.h>
#include <math_h/tabledata.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include <trackprocessing.h>
#include <detectors.h>
#include <reconstruction.h>
#include <data.h>
#include <Reconstruction/forward.h>
#include "analyses.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
Axis Q_axis_over(const Analysis&res){
	return Axis([&res](){return 1000.0*He3eta.P2Q(res.PBeam());},-70.0,30.0,40);
}
const string dir_v_name="He3Forward_Vertices";
const string dir_r_name="He3Forward_Reconstruction";

void He3_X_analyse(Analysis&res){
    const Axis MM_vertex([&res]()->double{
        for(const auto&P:res.Vertex(0))if(P.particle==Particle::he3())
            return He3eta.MissingMass({{.index=0,.E=P.E,.theta=P.Th,.phi=P.Phi}},res.PBeam());
        return INFINITY;
    },0.0,0.8,800);
    const Axis Ev([&res]()->double{
	for(const auto&P:res.Vertex(0))
	    if(P.particle==Particle::he3())
		return P.E;
	return INFINITY;
    },0.1,0.6,500);
    const Axis Tv([&res]()->double{
	for(const Analysis::Kinematic&P:res.Vertex(0))
	    if(P.particle==Particle::he3())
		return P.Th*180.0/PI();
	return INFINITY;
    },3.5,9.0,550);
    res.Trigger(0).pre()<<make_shared<SetOfHists1D>(dir_v_name,"MissingMass-vertex",Q_axis_over(res),MM_vertex);
    res.Trigger(0).pre()<<make_shared<SetOfHists2D>(dir_v_name,"Kinematic-vertex",Q_axis_over(res),Ev,Tv);
		
    static long trackcount,he3count;
    res.Trigger(trigger_he3_forward.number).pre()<<(make_shared<ChainCheck>()
	<<make_shared<Hist1D>(dir_r_name,"0-Reference",Q_axis_over(res))
	<<[](){trackcount=0;he3count=0;return true;}
    );
    static particle_kine he3;
    res.Trigger(trigger_he3_forward.number).per_track()<<(make_shared<ChainOr>()
	<<(make_shared<ChainCheck>()
	    <<[](WTrack&T){return T.Type()==kFDC;}
	    <<[](){trackcount++;return true;}
	)
	<<(make_shared<ChainCheck>()
	    <<ForwardHe3Reconstruction("He3Forward",res,he3)
	    <<[](){he3count++;return true;}
	    <<make_shared<SetOfHists1D>(dir_r_name,"MissingMass",
		Q_axis_over(res),
		Axis([&res](){
    	            return He3eta.MissingMass({
        	        {.index=0,.E=he3.E,.theta=he3.theta,.phi=he3.phi}
            	    },res.PBeam());
		},0.0,0.8,800)
	    )
	    <<make_shared<SetOfHists2D>(dir_r_name,"Kinematic-reconstructed",
		Q_axis_over(res),
		Axis([](){return he3.E;},0.1,0.6,500),
		Axis([](){return (he3.theta*180.0)/PI();},3.5,9.0,550)
	    )
	)
    );
    res.Trigger(trigger_he3_forward.number).post()<<(make_shared<ChainCheck>()
	<<make_shared<Hist1D>(dir_r_name,"Forward_charged_tracks",Axis([]()->double{return trackcount;},-0.5,9.5,10))
	<<make_shared<Hist1D>(dir_r_name,"Forward_he3_tracks",Axis([]()->double{return he3count;},-0.5,9.5,10))
    );
}
