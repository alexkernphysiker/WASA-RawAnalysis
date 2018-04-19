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
#include "montecarlo.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
Axis Q_axis(const Analysis&res){
	return Axis([&res](){return 1000.0*He3eta.P2Q(res.PBeam());},-70.0,30.0,40);
}
struct track_info{LorentzVector<>L;double t;};
void qe_central_analysis(Analysis&res){
	static SortedPoints<double,track_info> tracks;
	static auto Ptotal=LorentzVector<>::zero();
	res.Trigger(trigger_elastic1.number).pre()
		<<[&res](){
	                const auto p=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass());
        	        const auto d=lorentz_byPM(Zero(),Particle::d().mass());
                	Ptotal=p+d;
			tracks.clear();
			return true;
		}
		<<make_shared<Hist1D>("quasielastic","0-Reference",Q_axis(res))
	;
	res.Trigger(trigger_elastic1.number).per_track()<<(make_shared<ChainOr>()
		<<(make_shared<ChainCheck>()
			<<[](WTrack&T){return T.Type()==kCDC;}
			<<[](WTrack&T){return T.Theta()*180./PI()>23.;}
			<<[](WTrack&T){return T.Edep()>0.03;}
			<<[](WTrack&T){
				double phi=T.Phi();
				tracks<<point<double,track_info>(
					T.Theta(),
					{.L=lorentz_byEkM(double(T.Edep()),Particle::p().mass(),direction(phi,T.Theta())),.t=T.Time()}
				);
				return true;
			}
		)
	);
	static SortedPoints<double,pair<track_info,track_info>> trackpairs;
	static auto ct_axis=make_pair(
		Axis([](){return direction(trackpairs[0].Y().first.L.P()).th()*180./PI();},0.,180.,180),
		Axis([](){return direction(trackpairs[0].Y().second.L.P()).th()*180./PI();},0.,180.,180)
	);
	static auto time_axis=Axis([](){return trackpairs[0].Y().second.t-trackpairs[0].Y().first.t;},-50.,50.,100);
	static auto coplanarity=Axis([](){
		auto dphi=
			direction(trackpairs[0].Y().first.L.P()).phi()
			-direction(trackpairs[0].Y().second.L.P()).phi();
		dphi*=180./PI();
		while(dphi<0.)dphi+=360;
		while(dphi>360.)dphi-=360;
		return dphi;
	},0.0,360.0,360);
	static const auto ed_axis=make_pair(
		Axis([](){return trackpairs[0].Y().first.L.Ekin();} ,0.0,0.5,50),
		Axis([](){return trackpairs[0].Y().second.L.Ekin();},0.0,0.5,50)
	);
	static const Axis mm_axis([](){
		const auto& p1=trackpairs[0].Y().first;
		const auto& p2=trackpairs[0].Y().second;
		return (Ptotal-p1.L-p2.L).M();		
	},0.0,2.0,2000);
	res.Trigger(trigger_elastic1.number).post()<<(make_shared<ChainCheck>()
		<<make_shared<Hist1D>("quasielastic","Vectors",Axis([]()->double{return tracks.size();},-0.5,9.5,10))
		<<[](){
			trackpairs.clear();
			for(size_t i=0;i<tracks.size();i++)for(size_t j=(i+1);j<tracks.size();j++){
                                double dphi=direction(tracks[i].Y().L.P()).phi()
					-direction(tracks[j].Y().L.P()).phi();
                                dphi*=180./PI();
                                while(dphi<0.)dphi+=360;while(dphi>360.)dphi-=360;
				if(abs(dphi-180)<90){
					const bool C=(direction(tracks[i].Y().L.P()).th()<direction(tracks[j].Y().L.P()).th());
					const auto&T1=C?tracks[i]:tracks[j];
					const auto&T2=C?tracks[j]:tracks[i];
					trackpairs<<make_point(abs(dphi-180),make_pair(T1.Y(),T2.Y()));
				}
			}
			return trackpairs.size()>0;
		}
		<<make_shared<SetOfHists1D>("quasielastic","pair_phi_diff_0",Q_axis(res),coplanarity)
		<<make_shared<SetOfHists1D>("quasielastic","pair_time_diff_0",Q_axis(res),time_axis)
		<<make_shared<Hist2D>("quasielastic","t_vs_e_0",ct_axis.first,ed_axis.first)
		<<make_shared<Hist2D>("quasielastic","t_vs_t_0",ct_axis.first,ct_axis.second)
		<<make_shared<Hist2D>("quasielastic","e_vs_e_0",ed_axis.first,ed_axis.second)
		<<make_shared<Hist1D>("quasielastic","pp_mm_0",mm_axis)

		<<[](WTrack&T){return (ct_axis.first(T)<40.)&&(ct_axis.second(T)>45.);}
		<<make_shared<SetOfHists1D>("quasielastic","pair_phi_diff_1",Q_axis(res),coplanarity)
		<<make_shared<SetOfHists1D>("quasielastic","pair_time_diff_1",Q_axis(res),time_axis)
		<<make_shared<Hist2D>("quasielastic","t_vs_e_1",ct_axis.first,ed_axis.first)
		<<make_shared<Hist2D>("quasielastic","t_vs_t_1",ct_axis.first,ct_axis.second)
		<<make_shared<Hist2D>("quasielastic","e_vs_e_1",ed_axis.first,ed_axis.second)
		<<make_shared<Hist1D>("quasielastic","pp_mm_1",mm_axis)

		<<[](WTrack&T){return mm_axis(T)>1.1;}
                <<make_shared<SetOfHists1D>("quasielastic","pair_phi_diff_2",Q_axis(res),coplanarity)
                <<make_shared<SetOfHists1D>("quasielastic","pair_time_diff_2",Q_axis(res),time_axis)
                <<make_shared<Hist2D>("quasielastic","t_vs_e_2",ct_axis.first,ed_axis.first)
                <<make_shared<Hist2D>("quasielastic","t_vs_t_2",ct_axis.first,ct_axis.second)
                <<make_shared<Hist2D>("quasielastic","e_vs_e_2",ed_axis.first,ed_axis.second)
                <<make_shared<Hist1D>("quasielastic","pp_mm_2",mm_axis)
 
		<<[&res](WTrack&T){
			if(dynamic_cast<const MonteCarlo*>(&res))return true;
			return (time_axis(T)>-20.)&&(time_axis(T)<-5.);
		}
		<<make_shared<SetOfHists1D>("quasielastic","pair_phi_diff_3",Q_axis(res),coplanarity)
		<<make_shared<SetOfHists1D>("quasielastic","pair_time_diff_3",Q_axis(res),time_axis)
		<<make_shared<Hist2D>("quasielastic","t_vs_e_3",ct_axis.first,ed_axis.first)
		<<make_shared<Hist2D>("quasielastic","t_vs_t_3",ct_axis.first,ct_axis.second)
		<<make_shared<Hist2D>("quasielastic","e_vs_e_3",ed_axis.first,ed_axis.second)
		<<make_shared<Hist1D>("quasielastic","pp_mm_3",mm_axis)

                <<[](WTrack&T){return (ct_axis.first(T)<35);}
                <<make_shared<SetOfHists1D>("quasielastic","pair_phi_diff_4",Q_axis(res),coplanarity)
                <<make_shared<SetOfHists1D>("quasielastic","pair_time_diff_4",Q_axis(res),time_axis)
                <<make_shared<Hist2D>("quasielastic","t_vs_e_4",ct_axis.first,ed_axis.first)
                <<make_shared<Hist2D>("quasielastic","t_vs_t_4",ct_axis.first,ct_axis.second)
                <<make_shared<Hist2D>("quasielastic","e_vs_e_4",ed_axis.first,ed_axis.second)
                <<make_shared<Hist1D>("quasielastic","pp_mm_4",mm_axis)
                <<[](WTrack&T){return (ct_axis.first(T)<30);}
                <<make_shared<SetOfHists1D>("quasielastic","pair_phi_diff_5",Q_axis(res),coplanarity)
                <<make_shared<SetOfHists1D>("quasielastic","pair_time_diff_5",Q_axis(res),time_axis)
                <<make_shared<Hist2D>("quasielastic","t_vs_e_5",ct_axis.first,ed_axis.first)
                <<make_shared<Hist2D>("quasielastic","t_vs_t_5",ct_axis.first,ct_axis.second)
                <<make_shared<Hist2D>("quasielastic","e_vs_e_5",ed_axis.first,ed_axis.second)
                <<make_shared<Hist1D>("quasielastic","pp_mm_5",mm_axis)
 
	);
}
