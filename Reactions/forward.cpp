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

const double e_corr=0.007;

	using namespace std;
	using namespace MathTemplates;
	using namespace TrackAnalyse;
	const Reaction& He3eta(){
	    static Reaction res(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	    return res;
	}
	Axis Q_axis_over(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta().P2Q(res.PBeam());},0.0,30.0,12);}

	inline bool IsIn(double value,const pair<double,double>&&border){
		return (value>=border.first)&&(value<=border.second);
	}
	string dir_v_name(){return "He3Forward_Vertices";};
	string dir_r_name(){return "He3Forward_Reconstruction";};
	string dir_dbg_name(){return "He3Forward_Debug";};
	
	
	Axis MM_vertex(const Analysis&res){
		return Axis([&res]()->double{
			for(const auto&P:res.Vertex(0))if(P.particle==Particle::he3())
				return He3eta().MissingMass({{.index=0,.E=P.E,.theta=P.Th,.phi=P.Phi}},res.PBeam());
			return INFINITY;
		},0.0,0.8,800);
	}
	Axis Th_deg([](const vector<double>&P)->double{return P[0];},3.5,9.0,550);
	Axis Phi_deg([](const vector<double>&P)->double{return P[1];},0.0,360.0,360);
	Axis Ek_GeV([](const vector<double>&P)->double{return P[2];},0.1,0.6,500);
	Axis MM_GeV([](const vector<double>&P)->double{return P[3];},0.0,0.8,800);
	Axis Ev(const Analysis&data){return Axis([&data]()->double{
		for(const auto&P:data.Vertex(0))if(P.particle==Particle::he3())
			return P.E;
		return INFINITY;
	},Ek_GeV);}
	Axis Tv(const Analysis&data){return Axis([&data]()->double{
		for(const Analysis::Kinematic&P:data.Vertex(0))if(P.particle==Particle::he3())
			return P.Th*180.0/PI();
		return INFINITY;
	},Th_deg);}

	shared_ptr<AbstractChain> ForwardReconstructionProcess(const Analysis&data){
		return make_shared<ChainCheck>()
		<<Forward::Get().CreateMarker(dir_r_name(),"1-AllTracks")
		<<make_shared<Hist1D>(dir_r_name(),"1-AllTracks",Q_axis_over(data))
		<<[](WTrack&T)->bool{return (T.Theta()!=0.125);}
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Theta()*180.0/PI();})
		<<make_shared<Parameter>([](WTrack&T)->double{return NormPhi(T.Phi())*180.0/PI();})
		<<make_shared<SetOfHists1D>(dir_dbg_name(),"2-PhiDistribution",Q_axis_over(data),Phi_deg)
		<<Forward::Get().CreateMarker(dir_r_name(),"2-FPC")<<make_shared<Hist1D>(dir_r_name(),"2-FPC",Q_axis_over(data))
		<<(make_shared<ChainOr>()
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH1;}
				<<[](WTrack&T)->bool{
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
				<<Forward::Get().CreateMarker(dir_r_name(),"2.5-FRH1")
				<<make_shared<Parameter>([&data](WTrack&track)->double{
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
						"He3.E.FRH1",{[](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},[](WTrack&track){return track.Theta();}},
						[&data](WTrack&)->double{
							for(const auto&P:data.Vertex(0))if(P.particle==Particle::he3())
								return P.E;
							return INFINITY;
						}
					);
					return energy(track)+(dynamic_cast<const RealData*>(&data)?e_corr:0.0);
				})
			)
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T)->bool{return false;}//if I need FRH2 I'll delete this line
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH2;}
				<<[](WTrack&T)->bool{
					if(Forward::Get()[kFRH2].Edep(T)>0.22)return false;
					double locusline=0.3-0.417*Forward::Get()[kFRH2].Edep(T);
					return IsIn(Forward::Get()[kFRH1].Edep(T),make_pair(locusline-0.05,locusline+0.05));
				}
				<<Forward::Get().CreateMarker(dir_r_name(),"2.5-FRH2")
				<<make_shared<Parameter>([&data](WTrack&track)->double{
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH2,WTrack&> energy(
						"He3.E.FRH2",{[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)+Forward::Get()[kFRH2].Edep(T);},[](WTrack&T){return T.Theta();}},
						[&data](WTrack&)->double{
							for(const auto&P:data.Vertex(0))if(P.particle==Particle::he3())
								return P.E;
							return INFINITY;
						}
					);
					return energy(track)+(dynamic_cast<const RealData*>(&data)?e_corr:0.0);
				})	
			)
		)
		<<make_shared<SetOfHists1D>(dir_dbg_name(),"3-PhiDistribution",Q_axis_over(data),Phi_deg)
		<<Forward::Get().CreateMarker(dir_r_name(),"3-AllCuts")<<make_shared<Hist1D>(dir_r_name(),"3-AllCuts",Q_axis_over(data))
		<<[](const vector<double>&P)->bool{return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);}
		<<make_shared<SetOfHists1D>(dir_dbg_name(),"4-PhiDistribution",Q_axis_over(data),Phi_deg)
		<<Forward::Get().CreateMarker(dir_r_name(),"4-Reconstructed")
		<<make_shared<Hist1D>(dir_r_name(),"4-Reconstructed",Q_axis_over(data));
	}
	shared_ptr<AbstractChain> He3MissingMass(const Analysis&data){
		return make_shared<Chain>()
		<<make_shared<Parameter>([&data](WTrack&T,const vector<double>&P)->double{
			return He3eta().MissingMass({{.index=0,.E=Ek_GeV(T,P),.theta=T.Theta(),.phi=T.Phi()}},data.PBeam());
		})
		<<make_shared<SetOfHists1D>(dir_r_name(),"MissingMass",Q_axis_over(data),MM_GeV);
	}
	shared_ptr<AbstractChain> He3KinematicHe3Test(const Analysis&data){
		return make_shared<Chain>()
			<<make_shared<SetOfHists2D>(dir_r_name(),"Kinematic-reconstructed",Q_axis_over(data),Ek_GeV,Th_deg);
	}
	
	
	void He3_X_analyse(Analysis&res){
		res.Trigger(0).pre()<<make_shared<SetOfHists1D>(dir_v_name(),"MissingMass-vertex",Q_axis_over(res),MM_vertex(res));
		res.Trigger(0).pre()<<make_shared<SetOfHists2D>(dir_v_name(),"Kinematic-vertex",Q_axis_over(res),Ev(res),Tv(res));
		
		auto trackcount=make_shared<long>(0);
		res.Trigger(trigger_he3_forward.number).pre()<<make_shared<Hist1D>(dir_r_name(),"0-Reference",Q_axis_over(res))
			<<[trackcount](){(*trackcount)=0;return true;};
		res.Trigger(trigger_he3_forward.number).per_track()<<(make_shared<ChainCheck>()
			<<[](WTrack&T)->bool{return T.Type()==kFDC;}
			<<[trackcount](){(*trackcount)++;return true;}
			<<ForwardReconstructionProcess(res)
			<<[](WTrack&T,const vector<double>&P)->bool{
				return (Ek_GeV(T,P)<0.45)&&(Ek_GeV(T,P)>0.15)&&(Th_deg(T,P)<9.0);
			}
			<<He3MissingMass(res)
			<<He3KinematicHe3Test(res)
		);
		res.Trigger(trigger_he3_forward.number).post()
			<<make_shared<Hist1D>(
				dir_dbg_name(),
				"Forward_charged_tracks",
				Axis([trackcount]()->double{return *trackcount;},-0.5,9.5,10)
			);
	}
	void He3_X_reconstruction(Analysis&res){
		res.Trigger(trigger_he3_forward.number).pre()<<make_shared<Hist1D>(dir_r_name(),"0-Reference",Q_axis_over(res));
		res.Trigger(trigger_he3_forward.number).per_track()<<(make_shared<ChainCheck>()
			<<[](WTrack&T)->bool{return T.Type()==kFDC;}
			<<ForwardReconstructionProcess(res)
			<<He3KinematicHe3Test(res)
		);
	}
