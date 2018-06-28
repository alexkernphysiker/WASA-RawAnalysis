// this file is distributed under
// GPL license
#include <iostream>
#include <Wasa.hh>
#include <TCutG.h>
#include <math_h/functions.h>
#include <math_h/tabledata.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include <Parameters/parameters.h>
#include <trackprocessing.h>
#include <detectors.h>
#include <reconstruction.h>
#include <data.h>
#include "forward.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;
particle_kine testmode;
particle_kine&___test_mode___(){return testmode;}
shared_ptr<AbstractChain> ForwardHe3Reconstruction(const string&prefix,const Analysis&data,particle_kine&kin_rec){
    const string dir_v_name=prefix+"_Vertices";
    const string dir_r_name=prefix+"_Reconstruction";
    const string dir_dbg_name=prefix+"_Debug";
    static const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
    const Axis Q_axis_over([&data](){return 1000.0*He3eta.P2Q(data.PBeam());},0.0,30.0,12);
    const Axis Phi_deg([&kin_rec](){return kin_rec.phi*180./PI();},-180.0,180.0,360);
    const function<double(WTrack&)> vertex_energy=[&data](WTrack&)->double{
        for(const auto&P:data.Vertex(0)){
            if(P.particle==Particle::he3())return P.E;
        }
        return INFINITY;
    };
    auto corr_hists=make_shared<ChainOr>()<<[](){return true;};
    if(&kin_rec==&testmode){
	const Axis Theta_lr([&kin_rec](){return kin_rec.theta*180./PI();},0.0,20.0,20);
	const Axis Edep([](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},0.0,0.4,200);
	const Axis Erec([&kin_rec](){return kin_rec.E;},0.0,0.6,400);
	const Axis Ever(vertex_energy,0.0,0.6,400);
	corr_hists 
	    <<make_shared<SetOfHists2D>("dir_dbg_name","2-Edep-vs-Vertex",Theta_lr,Edep,Ever)
	    <<make_shared<SetOfHists2D>("dir_dbg_name","2-Erec-vs-Vertex",Theta_lr,Erec,Ever)
	;
    }
    return make_shared<ChainCheck>()
	<<[](WTrack&T){return T.Type()==kFDC;}
        <<Forward::Get().CreateMarker(dir_r_name,"0-AllTracks")
        <<make_shared<Hist1D>(dir_r_name,"0-AllTracks",Q_axis_over)
        <<[&kin_rec](WTrack&T){
		kin_rec.particle=Particle::he3();
		return true;
	}
        <<Forward::Get().CreateMarker(dir_r_name,"1-PhiIsFinite")
        <<make_shared<SetOfHists1D>(dir_dbg_name,"1-PhiDistribution",Q_axis_over,Phi_deg)
        <<make_shared<Hist1D>(dir_r_name,"1-PhiIsFinite",Q_axis_over)
	<<[](WTrack&T){return (T.Theta()!=0.125);}
	<<[](WTrack&T){return (T.Theta()>(PI()*(4.5/180.0)));}//because events amount with theta<4 deg is not in agreement with MC
        <<[&kin_rec](WTrack&T){
	    kin_rec.theta=T.Theta();
	    return true;
	}
        <<make_shared<SetOfHists1D>(dir_dbg_name,"2-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"2-ThetaIsAccepted")
	<<make_shared<Hist1D>(dir_r_name,"2-ThetaIsAccepted",Q_axis_over)
        <<[](WTrack&T){return Forward::Get().StoppingLayer(T)==kFRH1;}
        <<make_shared<SetOfHists1D>(dir_dbg_name,"3-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"3-StopInFRH1")
	<<make_shared<Hist1D>(dir_r_name,"3-StopInFRH1",Q_axis_over)
        <<[](WTrack&T){
            const double x=Forward::Get()[kFRH1].Edep(T);
            const double y=Forward::Get()[kFTH1].Edep(T);
	    const double h=getParameter(he3_cut_h);
	    return (y>(h+((0.1-x)*0.05)))&&(y>(h-((x-0.1)*0.02)));
        }
        <<make_shared<SetOfHists1D>(dir_dbg_name,"4-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"4-GeomCut")
	<<make_shared<Hist1D>(dir_r_name,"4-GeomCut",Q_axis_over)
        <<[&data,&kin_rec,vertex_energy](WTrack&track){
            static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
                "He3.E.FRH1",
		{
		    [](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},
		    [](WTrack&track){return track.Theta();}
		},
		vertex_energy
            );
            kin_rec.E=energy(track);
            double phi=track.Phi();
	    if(dynamic_cast<const RealData*>(&data)){
		const double e=kin_rec.E,m=Particle::he3().mass(),p=sqrt(pow(m+e,2)-pow(m,2));
                phi+=data.FinderFD().GetPhiCorrection(p,phi,2,10.); //MF (kG),                                                      
	    }
            kin_rec.phi=phi;
	    return true;
        }
        <<make_shared<SetOfHists1D>(dir_dbg_name,"5-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"5-Reconstructed")
        <<make_shared<Hist1D>(dir_r_name,"5-Reconstructed",Q_axis_over)
        <<corr_hists
    ;
}
