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
Axis Q_axis(const Analysis&res){
	return Axis([&res](){return 1000.0*He3eta.P2Q(res.PBeam());},-70.0,30.0,40);
}
void p_or_d_analyse(Analysis&res){
	static SortedPoints<double,particle_kinematics> tracks;
	static size_t FC,CC;
	res.Trigger(21).pre()
		<<[](){FC=0;CC=0;return true;}
		<<[](){tracks.clear();return true;}
		<<make_shared<Hist1D>("elastic","0-Reference",Q_axis(res))
	;
	res.Trigger(21).per_track()<<(make_shared<ChainOr>()
		<<(make_shared<ChainCheck>()
			<<[](WTrack&T){return T.Type()==kFDC;}
			<<[](WTrack&T){return T.Theta()!=0.125;}
			<<[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)<0.035;}
			<<[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)>0.020;}
			<<[](WTrack&T){return Forward::Get()[kFRH2].Edep(T)<0.045;}
			<<[](WTrack&T){return Forward::Get()[kFRH2].Edep(T)>0.020;}
			<<[](WTrack&T){return T.Theta()<0.30;}
			<<[](){FC++;return true;}
			<<[](WTrack&T){
				tracks<<point<double,particle_kinematics>(
					T.Theta(),
					{.particle=Particle::p(),.E=Forward::Get()[kFRH1].Edep(T)+Forward::Get()[kFRH2].Edep(T),.theta=T.Theta(),.phi=T.Phi()}
				);
				return true;
			}
		)
		<<(make_shared<ChainCheck>()
			<<[](WTrack&T){return T.Type()==kCDC;}
			<<[](WTrack&T){return T.Theta()>0.40;}
			<<[](){CC++;return true;}
			<<[](WTrack&T){
				tracks<<point<double,particle_kinematics>(
					T.Theta(),
					{.particle=Particle::p(),.E=T.Edep(),.theta=T.Theta(),.phi=T.Phi()}
				);
				return true;
			}
		)
	);
	static SortedPoints<double,pair<particle_kinematics,particle_kinematics>> trackpairs;
	static auto ct_axis=make_pair(
		Axis([](){return trackpairs[0].Y().first.theta*180./PI();},0,90,90),
		Axis([](){return trackpairs[0].Y().second.theta*180./PI();},0,180,180)
	);
	static auto ct_sum=Axis([](){return (trackpairs[0].Y().first.theta+trackpairs[0].Y().second.theta*0.5)*180./PI();},0,180,180);
	static const auto ed_axis=make_pair(
		Axis([](){return trackpairs[0].Y().first.E;} ,0.0,0.5,50),
		Axis([](){return trackpairs[0].Y().second.E;},0.0,0.5,50)
	);
	res.Trigger(21).post()<<(make_shared<ChainCheck>()
		<<make_shared<Hist1D>("elastic","Forward_charged_tracks",Axis([]()->double{return FC;},-0.5,9.5,10))
		<<make_shared<Hist1D>("elastic","Central_charged_tracks",Axis([]()->double{return CC;},-0.5,9.5,10))
		<<make_shared<Hist1D>("elastic","Vectors",Axis([]()->double{return tracks.size();},-0.5,9.5,10))
		<<[](){
			trackpairs.clear();
			for(size_t i=0;i<tracks.size();i++)for(size_t j=(i+1);j<tracks.size();j++){
				trackpairs<<point<double,pair<particle_kinematics,particle_kinematics>>(
					sqrt(pow(tracks[i].Y().phi-tracks[j].Y().phi-PI(),2))*180/PI(),
					make_pair(tracks[i].Y(),tracks[j].Y())
				);
			}
			return trackpairs.size()>0;
		}
		<<make_shared<Hist1D>("elastic","pair_phi_diff_0",Axis([](){return trackpairs[0].X();},0.0,90.0,90))
		<<make_shared<Hist1D>("elastic","count_0",Q_axis(res))
		<<[](){return trackpairs[0].X()<15.0;}
		<<make_shared<Hist1D>("elastic","pair_phi_diff_1",Axis([](){return trackpairs[0].X();},0.0,90.0,90))
		<<make_shared<Hist2D>("elastic","t_vs_e_1",ct_axis.first,ed_axis.first)
		<<make_shared<Hist2D>("elastic","t_vs_t_1",ct_axis.first,ct_axis.second)
		<<make_shared<Hist2D>("elastic","e_vs_e_1",ed_axis.first,ed_axis.second)
		<<make_shared<Hist1D>("elastic","count_1",Q_axis(res))
		<<(make_shared<ChainOr>()
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T){return((ct_axis.first(T)>20.0)&&(ct_axis.first(T)<40.0));}
				<<make_shared<Hist1D>("elastic","pair_phi_diff_21",Axis([](){return trackpairs[0].X();},0.0,90.0,90))
				<<make_shared<Hist2D>("elastic","t_vs_e_21",ct_axis.first,ed_axis.first)
				<<make_shared<Hist2D>("elastic","t_vs_t_21",ct_axis.first,ct_axis.second)
				<<make_shared<Hist2D>("elastic","e_vs_e_21",ed_axis.first,ed_axis.second)
				<<make_shared<Hist1D>("elastic","theta_sum_21",ct_sum)
				<<make_shared<Hist1D>("elastic","count_21",Q_axis(res))
			)
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T){return(ct_axis.first(T)<20.0);}
				<<make_shared<Hist1D>("elastic","pair_phi_diff_22",Axis([](){return trackpairs[0].X();},0.0,90.0,90))
				<<make_shared<Hist2D>("elastic","t_vs_e_22",ct_axis.first,ed_axis.first)
				<<make_shared<Hist2D>("elastic","t_vs_t_22",ct_axis.first,ct_axis.second)
				<<make_shared<Hist2D>("elastic","e_vs_e_22",ed_axis.first,ed_axis.second)
				<<make_shared<Hist1D>("elastic","theta_sum_22",ct_sum)
				<<make_shared<Hist1D>("elastic","count_22",Q_axis(res))
			)
		)
	);
}
