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
#include "forward.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

shared_ptr<AbstractChain> ForwardHe3Reconstruction(const Analysis&data,particle_kinematics&kin_rec){
    static const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
    static const string dir_v_name="He3Forward_Vertices";
    static const string dir_r_name="He3Forward_Reconstruction";
    static const string dir_dbg_name="He3Forward_Debug";
    const Axis Q_axis_over([&data](){return 1000.0*He3eta.P2Q(data.PBeam());},0.0,30.0,12);
    const Axis Phi_deg([&kin_rec](){return kin_rec.phi;},0.0,360.0,360);
    return make_shared<ChainCheck>()

	<<[](WTrack&T){return T.Type()==kFDC;}
        <<Forward::Get().CreateMarker(dir_r_name,"0-AllTracks")
        <<make_shared<Hist1D>(dir_r_name,"0-AllTracks",Q_axis_over)

        <<[&kin_rec](WTrack&T){
		kin_rec.particle=Particle::he3();
		kin_rec.phi=T.Phi();
		return true;
	}
        <<Forward::Get().CreateMarker(dir_r_name,"1-PhiIsFinite")
        <<make_shared<SetOfHists1D>(dir_dbg_name,"1-PhiDistribution",Q_axis_over,Phi_deg)
        <<make_shared<Hist1D>(dir_r_name,"1-PhiIsFinite",Q_axis_over)

	<<[](WTrack&T){return (T.Theta()!=0.125);}
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
            static TCutG *cut=nullptr;
            if(cut==nullptr){
                cut=new TCutG("FRH1_cut",16);
                cut->SetVarX("FRH1");
                cut->SetVarY("FTH1");
                cut->SetPoint(16,0.018,0.025);
    	        cut->SetPoint(15,0.019,0.030);
                cut->SetPoint(14,0.069,0.021);
                cut->SetPoint(13,0.121,0.018);
                cut->SetPoint(12,0.162,0.016);
                cut->SetPoint(11,0.206,0.015);
                cut->SetPoint(10,0.304,0.014);
                cut->SetPoint(9,0.35,0.014);
                cut->SetPoint(8,0.35,0.006);
                cut->SetPoint(7,0.298,0.006);
                cut->SetPoint(6,0.201,0.007);
                cut->SetPoint(5,0.141,0.009);
                cut->SetPoint(4,0.105,0.011);
                cut->SetPoint(3,0.061,0.014);
                cut->SetPoint(2,0.027,0.019);
                cut->SetPoint(1,0.018,0.025);
            }
            double x=Forward::Get()[kFRH1].Edep(T);
            double y=Forward::Get()[kFTH1].Edep(T);
            return cut->IsInside(x,y);
        }
        <<make_shared<SetOfHists1D>(dir_dbg_name,"4-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"4-GeomCut")
	<<make_shared<Hist1D>(dir_r_name,"4-GeomCut",Q_axis_over)

        <<[&data,&kin_rec](WTrack&track){
            static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
                "He3.E.FRH1",
		{
		    [](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},
		    [](WTrack&track){return track.Theta();}
		},
                [&data](WTrack&)->double{
                    for(const auto&P:data.Vertex(0)){
                	if(P.particle==Particle::he3())return P.E;
                    }
                    return INFINITY;
                }
            );
            kin_rec.E=energy(track)+(dynamic_cast<const RealData*>(&data)?0.007:0.0);
	    return true;
        }
        <<make_shared<SetOfHists1D>(dir_dbg_name,"5-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"5-Reconstructed")
        <<make_shared<Hist1D>(dir_r_name,"5-Reconstructed",Q_axis_over)
    ;
}
