// this file is distributed under 
// MIT license
#include <Wasa.hh>
#include <TCutG.h>
#include "math_h/functions.h"
#include "math_h/tabledata.h"
#include "ReconstructionFit/reconstruction_types.h"
#include "Experiment/experiment_conv.h"
#include "Kinematics/reactions.h"
#include "trackprocessing.h"
#include "detectors.h"
#include "reconstruction.h"
#include "data.h"
#include "montecarlo.h"
#include "he3_X.h"
namespace ReactionSetup{
	using namespace std;
	using namespace MathTemplates;
	using namespace TrackAnalyse;
	inline bool IsIn(double value,const pair<double,double>&&border){
		return (value>=border.first)&&(value<=border.second);
	}
	string dir_v_name(){return "He3Forward_Vertices";};
	string dir_r_name(){return "He3Forward_Reconstruction";};
	string dir_dbg_name(){return "He3Forward_Debug";};
	Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	Axis Q_axis(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta.P2Q(res.PBeam());},0.0,30.0,12);}
	Axis MM_vertex(const Analysis&res){
		return Axis([&res]()->double{
			Analysis::Kinematic p=res.FromFirstVertex(kHe3);
			return He3eta.MissingMass({{.index=0,.E=p.E,.theta=p.Th,.phi=p.Phi}},res.PBeam());
		},0.0,0.8,800);
	}
	Axis Th_deg([](const vector<double>&P)->double{return P[0];},3.5,9.0,550);
	Axis Phi_deg([](const vector<double>&P)->double{return P[1];},0.0,360.0,360);
	Axis Ek_GeV([](const vector<double>&P)->double{return P[2];},0.1,0.6,500);
	Axis MM_GeV([](const vector<double>&P)->double{return P[3];},0.0,0.8,800);
	Axis Ev(const Analysis&data){return Axis([&data]()->double{return data.FromFirstVertex(kHe3).E;},Ek_GeV);}
	Axis Tv(const Analysis&data){return Axis([&data]()->double{return data.FromFirstVertex(kHe3).Th*180.0/PI();},Th_deg);}
	Analysis*Prepare(He3Modification mode){
		Analysis* res=nullptr;
		switch(mode){
			case forData:
				return new RealData();
			case forEta:
				res=new MonteCarlo();
				res->AddParticleToFirstVertex(kHe3,Particle::he3().mass());
				res->AddParticleToFirstVertex(kEta,Particle::eta().mass());
				break;
			case forPi0:
				res=new MonteCarlo();
				res->AddParticleToFirstVertex(kHe3,Particle::he3().mass());
				res->AddParticleToFirstVertex(kPi0,Particle::pi0().mass());
				break;
		};
		res->Trigger(0).pre()<<make_shared<SetOfHists1D>(dir_v_name(),"MissingMass-vertex",Q_axis(*res),MM_vertex(*res));
		res->Trigger(0).pre()<<make_shared<SetOfHists2D>(dir_v_name(),"Kinematic-vertex",Q_axis(*res),Ev(*res),Tv(*res));
		return res;
	}

	shared_ptr<AbstractChain> ForwardReconstructionProcess(const Analysis&data){
		return make_shared<ChainCheck>()
		<<Forward::Get().CreateMarker(dir_r_name(),"1-AllTracks")
		<<make_shared<Hist1D>(dir_r_name(),"1-AllTracks",Q_axis(data))
		<<[](WTrack&T)->bool{return (T.Theta()!=0.125);}
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Theta()*180.0/PI();})
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Phi()*180.0/PI();})
		<<Forward::Get().CreateMarker(dir_r_name(),"2-FPC")<<make_shared<Hist1D>(dir_r_name(),"2-FPC",Q_axis(data))
		<<(make_shared<ChainOr>()
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH1;}
				<<[](WTrack&T)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("FRH1_cut",15);
						cut->SetVarX("FRH1");
						cut->SetVarY("FTH1");
						cut->SetPoint(15,0.019,0.025);
						cut->SetPoint(14,0.069,0.019);
						cut->SetPoint(13,0.121,0.016);
						cut->SetPoint(12,0.162,0.014);
						cut->SetPoint(11,0.206,0.013);
						cut->SetPoint(10,0.304,0.012);
						cut->SetPoint(9,0.4,0.012);
						cut->SetPoint(8,0.4,0.008);
						cut->SetPoint(7,0.298,0.008);
						cut->SetPoint(6,0.201,0.010);
						cut->SetPoint(5,0.141,0.012);
						cut->SetPoint(4,0.105,0.013);
						cut->SetPoint(3,0.061,0.017);
						cut->SetPoint(2,0.027,0.021);
						cut->SetPoint(1,0.019,0.025);
					}
					double x=Forward::Get()[kFRH1].Edep(T);
					double y=Forward::Get()[kFTH1].Edep(T);
					return cut->IsInside(x,y);
				}
				<<Forward::Get().CreateMarker(dir_r_name(),"2.5-FRH1")
				<<make_shared<Parameter>([&data](WTrack&track)->double{
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
						"He3.E.FRH1",{[](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},[](WTrack&track){return track.Theta();}},
						[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
					);
					return energy(track);
				})
			)
			<<(make_shared<ChainCheck>()
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
						[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
					);
					return energy(track);
				})	
			)
		)
		<<Forward::Get().CreateMarker(dir_r_name(),"3-AllCuts")<<make_shared<Hist1D>(dir_r_name(),"3-AllCuts",Q_axis(data))
		<<[](const vector<double>&P)->bool{return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);}
		<<Forward::Get().CreateMarker(dir_r_name(),"4-Reconstructed")
		<<make_shared<Hist1D>(dir_r_name(),"4-Reconstructed",Q_axis(data));
	}
	shared_ptr<AbstractChain> He3MissingMass(const Analysis&data){
		return make_shared<Chain>()
		<<make_shared<Parameter>([&data](WTrack&T,const vector<double>&P)->double{
			return He3eta.MissingMass({{.index=0,.E=Ek_GeV(T,P),.theta=T.Theta(),.phi=T.Phi()}},data.PBeam());
		})
		<<make_shared<SetOfHists1D>(dir_r_name(),"MissingMass",Q_axis(data),MM_GeV);
	}
	shared_ptr<AbstractChain> He3KinematicHe3Test(const Analysis&data){
		return make_shared<Chain>()<<make_shared<SetOfHists2D>(dir_r_name(),"Kinematic-reconstructed",Q_axis(data),Ek_GeV,Th_deg);
	}
	
	void SearchGammaTracks(Analysis* res){
		auto data=make_shared<vector<particle_kinematics>>();
		res->Trigger(17).pre()<<[data](){data->clear(); return true;};
		res->Trigger(17).per_track()<<(make_shared<ChainCheck>()
			<<[](WTrack&T)->bool{return T.Type()==kCDN;}
			<<[data](WTrack&T)->bool{
				data->push_back({.particle=Particle::gamma(),.E=T.Ek(),.theta=T.Theta(),.phi=T.Phi()});
				return true;
			}
		);
		auto im_val=make_shared<double>(INFINITY);
		res->Trigger(17).post()
			<<[data,im_val](){
				(*im_val)=INFINITY;
				SortedPoints<double> table;
				for(size_t i=0;i<data->size();i++)
					for(size_t j=i+1;j<data->size();j++){
						double im=InvariantMass({data->operator[](i),data->operator[](j)});
						table<<point<double>(pow(im-Particle::eta().mass(),2),im);
					}
				if(table.size()>0)
					(*im_val)=table[0].Y();
				return true;
			}
			<<make_shared<Hist1D>("CentralGammas","inv_mass_2gamma",Axis([im_val]()->double{return *im_val;},0.0,0.8,800));
	}
	
	Analysis* He3_X_analyse(He3Modification mode){
		auto res=Prepare(mode);
		auto trackcount=make_shared<long>(0);
		res->Trigger(trigger_he3_forward.number).pre()<<make_shared<Hist1D>(dir_r_name(),"0-Reference",Q_axis(*res))
			<<[trackcount](){(*trackcount)=0;return true;};
		res->Trigger(trigger_he3_forward.number).per_track()<<(make_shared<ChainCheck>()
			<<[](WTrack&T)->bool{return T.Type()==kFDC;}
			<<[trackcount](){(*trackcount)++;return true;}
			<<ForwardReconstructionProcess(*res)
			<<[](WTrack&T,const vector<double>&P)->bool{
				return (Ek_GeV(T,P)<0.45)&&(Ek_GeV(T,P)>0.15)&&(Th_deg(T,P)<9.0);
			}
			<<He3MissingMass(*res)
			<<He3KinematicHe3Test(*res)
		);
		res->Trigger(trigger_he3_forward.number).post()
			<<make_shared<Hist1D>(
				dir_dbg_name(),
				"Forward_charged_tracks",
				Axis([trackcount]()->double{return *trackcount;},-0.5,9.5,10)
			);
		SearchGammaTracks(res);
		return res;
	}
	Analysis* He3_X_reconstruction(He3Modification mode){
		auto res=Prepare(mode);
		res->Trigger(trigger_he3_forward.number).pre()<<make_shared<Hist1D>(dir_r_name(),"0-Reference",Q_axis(*res));
		res->Trigger(trigger_he3_forward.number).per_track()<<(make_shared<ChainCheck>()
			<<[](WTrack&T)->bool{return T.Type()==kFDC;}
			<<ForwardReconstructionProcess(*res)
			<<He3KinematicHe3Test(*res)
		);
		return res;
	}
}
