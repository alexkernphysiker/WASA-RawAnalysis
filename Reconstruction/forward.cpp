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
particle_kinematics testmode;
particle_kinematics&___test_mode___(){return testmode;}
shared_ptr<AbstractChain> ForwardHe3Reconstruction(const Analysis&data,particle_kinematics&kin_rec){
    static const string dir_v_name="He3Forward_Vertices";
    static const string dir_r_name="He3Forward_Reconstruction";
    static const string dir_dbg_name="He3Forward_Debug";
    static const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
    const Axis Q_axis_over([&data](){return 1000.0*He3eta.P2Q(data.PBeam());},0.0,30.0,12);
    const Axis Phi_deg([&kin_rec](){return kin_rec.phi*180./PI();},-180.0,180.0,360);
    const function<double(WTrack&)> vertex_energy=[&data](WTrack&)->double{
        for(const auto&P:data.Vertex(0)){
            if(P.particle==Particle::he3())return P.E;
        }
        return INFINITY;
    };
    auto corr_hists=make_shared<ChainOr>();
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
                cut=new TCutG("He3_cut",16);
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

        <<[&data,&kin_rec,vertex_energy](WTrack&track){
            static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
                "He3.E.FRH1",
		{
		    [](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},
		    [](WTrack&track){return track.Theta();}
		},
		vertex_energy
            );
            kin_rec.E=energy(track)+(dynamic_cast<const RealData*>(&data)?0.007:0.0);
	    return true;
        }
        <<make_shared<SetOfHists1D>(dir_dbg_name,"5-PhiDistribution",Q_axis_over,Phi_deg)
        <<Forward::Get().CreateMarker(dir_r_name,"5-Reconstructed")
        <<make_shared<Hist1D>(dir_r_name,"5-Reconstructed",Q_axis_over)
        <<corr_hists
    ;
}

shared_ptr<AbstractChain> ForwardDReconstruction(const Analysis&data,particle_kinematics&kin_rec){
    const Axis Ed([&data,&kin_rec]()->double{
        for(const auto&P:data.Vertex(0)){
	    if(P.particle==Particle::d()){
    	        if(
    	    	    (pow(P.Th-kin_rec.theta,2)<0.0003)&&
    	    	    (pow(P.Phi-kin_rec.phi,2)<0.001)
    	    	)return P.E;
            }
        }
        return INFINITY;
    },0.0,2.5,250);
    const Axis Td([&data,&kin_rec]()->double{
        for(const auto&P:data.Vertex(0)){
            if(P.particle==Particle::d()){
                if(
                    (pow(P.Th-kin_rec.theta,2)<0.0003)&&
                    (pow(P.Phi-kin_rec.phi,2)<0.001)
                )return P.Th*180.0/PI();
            }
        }
	return INFINITY;
    },0.0,20.0,200);
                                                                                                                                                
    const Axis Phi_deg([&kin_rec](){return kin_rec.phi*180./PI();},-180.0,180.0,360);
    const Axis Theta_deg([&kin_rec](){return kin_rec.theta*180./PI();},0.0,20.0,200);
    const Axis Theta_lr([&kin_rec](){return kin_rec.theta*180./PI();},0.0,20.0,20);
    const Axis Edep([](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},0.01,0.06,500);
    return make_shared<ChainCheck>()
        <<[](WTrack&T)->bool{return T.Type()==kFDC;}
        <<[&kin_rec](WTrack&T){
	    kin_rec.phi=T.Phi();
	    kin_rec.theta=T.Theta();
	    return true;
	}
	<<[](WTrack&T){return (T.Theta()!=0.125);}
        <<Forward::Get().CreateMarker("D","1-ChargedTracks")
	<<make_shared<Hist1D>("D","1-Phi",Phi_deg)
	<<make_shared<Hist2D>("D","1-Edep-vs-Theta",Edep,Theta_deg)

        <<[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)<0.050;}
        <<[](WTrack&T){
            static TCutG *cut=nullptr;
            if(cut==nullptr){
                cut=new TCutG("D_cut",16);
                cut->SetVarX("FRH2");
                cut->SetVarY("Theta");
                cut->SetPoint(16,0.042,17.94);
                cut->SetPoint(15,0.043,16.54);
                cut->SetPoint(14,0.041,13.70);
                cut->SetPoint(13,0.039,10.53);
                cut->SetPoint(12,0.037,6.37);
                cut->SetPoint(11,0.036,4.39);
                cut->SetPoint(10,0.034,3.81);
                cut->SetPoint(9,0.032,6.04);
                cut->SetPoint(8,0.032,8.93);
                cut->SetPoint(7,0.033,12.55);
                cut->SetPoint(6,0.033,15.14);
                cut->SetPoint(5,0.034,16.91);
                cut->SetPoint(4,0.035,17.74);
                cut->SetPoint(3,0.038,18.19);
                cut->SetPoint(2,0.040,18.19);
                cut->SetPoint(1,0.042,17.94);
            }
            const double x=Forward::Get()[kFRH2].Edep(T);
            const double y=T.Theta()*180./PI();
            return cut->IsInside(x,y);
        }
        <<Forward::Get().CreateMarker("D","2-GeomCut")
        <<make_shared<Hist1D>("D","2-Phi",Phi_deg)
        <<make_shared<Hist2D>("D","2-Edep-vs-Theta",Edep,Theta_deg)
	<<make_shared<Hist2D>("D","2-Edep-vs-Theta-true",Ed,Td)
	
        <<make_shared<Hist2D>("D","2-ThetaRec-vs-Vertex",Theta_deg,Td)
        <<make_shared<SetOfHists2D>("D","2-Edep-vs-Vertex",Theta_lr,Edep,Ed)
    ;
}
shared_ptr<AbstractChain> ForwardPReconstruction(const Analysis&data,particle_kinematics&kin_rec){
    const Axis Ep([&data,&kin_rec]()->double{
        for(const auto&P:data.Vertex(0)){
            if(P.particle==Particle::p()){
        	if(
		    (pow(P.Th-kin_rec.theta,2)<0.0003)&&
        	    (pow(P.Phi-kin_rec.phi,2)<0.001)
        	)return P.E;
    	    }
    	}
	return INFINITY;
    },0.0,2.5,250);
    const Axis Tp([&data,&kin_rec]()->double{
        for(const auto&P:data.Vertex(0)){
	    if(P.particle==Particle::p()){
		if(
		    (pow(P.Th-kin_rec.theta,2)<0.0003)&&
		    (pow(P.Phi-kin_rec.phi,2)<0.001)
		)return P.Th*180.0/PI();
	    }
	}
        return INFINITY;
    },0.0,20.0,200);
                                                                        
    const Axis Phi_deg([&kin_rec](){return kin_rec.phi*180./PI();},-180.0,180.0,360);
    const Axis Theta_deg([&kin_rec](){return kin_rec.theta*180./PI();},0.0,20.0,200);
    const Axis Theta_lr([&kin_rec](){return kin_rec.theta*180./PI();},0.0,20.0,20);
    const Axis Edep([](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},0.01,0.06,500);
    return make_shared<ChainCheck>()
        <<[](WTrack&T)->bool{return T.Type()==kFDC;}
        <<Forward::Get().CreateMarker("P","0-AllTracks")
        <<make_shared<Hist1D>("P","0-Phi",Phi_deg)
        <<make_shared<Hist2D>("P","0-Edep-vs-Theta",Edep,Theta_deg)

        <<[&kin_rec](WTrack&T){
	    kin_rec.phi=T.Phi();
	    kin_rec.theta=T.Theta();
	    return true;
	}
	<<[](WTrack&T){return (T.Theta()!=0.125);}
        <<Forward::Get().CreateMarker("P","1-Reconstructable")
	<<make_shared<Hist1D>("P","1-Phi",Phi_deg)
	<<make_shared<Hist2D>("P","1-Edep-vs-Theta",Edep,Theta_deg)

        <<[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)<0.035;}
        <<[](WTrack&T){
    	    static TCutG *cut=nullptr;
            if(cut==nullptr){
                cut=new TCutG("P_cut",16);
                cut->SetVarX("FRH2");
                cut->SetVarY("Theta");
                cut->SetPoint(16,0.030,17.50);
                cut->SetPoint(15,0.029,13.08);
                cut->SetPoint(14,0.029,10.00);
                cut->SetPoint(13,0.027,4.47);
                cut->SetPoint(12,0.026,3.44);
                cut->SetPoint(11,0.025,2.95);
                cut->SetPoint(10,0.024,3.32);
                cut->SetPoint(9,0.023,3.01);
                cut->SetPoint(8,0.022,4.76);
                cut->SetPoint(7,0.021,6.33);
                cut->SetPoint(6,0.021,9.74);
                cut->SetPoint(5,0.021,12.85);
                cut->SetPoint(4,0.022,16.14);
                cut->SetPoint(3,0.023,17.41);
                cut->SetPoint(2,0.025,18.14);
                cut->SetPoint(1,0.030,17.50);
            }
            const double x=Forward::Get()[kFRH2].Edep(T);
            const double y=T.Theta()*180./PI();
            return cut->IsInside(x,y);
        }
        <<Forward::Get().CreateMarker("P","2-GeomCut")
        <<make_shared<Hist1D>("P","2-Phi",Phi_deg)
        <<make_shared<Hist2D>("P","2-Edep-vs-Theta",Edep,Theta_deg)
        <<make_shared<Hist2D>("P","2-Edep-vs-Theta-true",Ep,Tp)

        <<make_shared<Hist2D>("P","2-ThetaRec-vs-Vertex",Theta_deg,Tp)
        <<make_shared<SetOfHists2D>("P","2-Edep-vs-Vertex",Theta_lr,Edep,Ep)
    ;
}
