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
#include "forward.h"
#include "central.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

Axis Q_axis_full(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta().P2Q(res.PBeam());},-70.0,30.0,40);}

void Search6Gamma(Analysis&res){
        auto data=make_shared<vector<particle_kinematics>>();
        const auto&tn=trigger_gammas_central.number;
	auto im_val=make_shared<point<double>>(INFINITY,INFINITY);
	res.Trigger(tn).pre()
		<< make_shared<Hist1D>("CentralGammas","0-Reference",Q_axis_full(res))
		<<[data](){data->clear();return true;}
	;
	res.Trigger(tn).per_track()
		<<(make_shared<ChainCheck>()
			<< [](WTrack&T)->bool{return T.Type()==kCDN;}
		        << [data](WTrack&T)->bool{
                		data->push_back({.particle=Particle::gamma(),.E=T.Edep(),.theta=T.Theta(),.phi=T.Phi()});
				return true;
			}
		)
	;
	res.Trigger(tn).post()
		<< make_shared<Hist1D>("CentralGammas","neutral_tracks_count",Axis([data]()->double{return data->size();},-0.5,9.5,10))
		<<(make_shared<ChainCheck>()
			<< [data](){return data->size()>=6;}
			<< [data,im_val]()->bool{
				vector<size_t> inc_track;
				for(size_t i=0;i<data->size();i++)inc_track.push_back(false);
				function<point<double>(const size_t)> look_for=[&inc_track,data,&look_for](const size_t i)->point<double>{
					if(i>data->size())throw Exception<vector<particle_kinematics>>("error in 6-gamma search");
					if(i<data->size()){
						SortedPoints<double> res;
						inc_track[i]=true;
						res<<look_for(i+1);
						inc_track[i]=false;
						res<<look_for(i+1);
						return res[0];
					}else{
						vector<particle_kinematics> used,notused;
						for(size_t i=0;i<data->size();i++){
							if(inc_track[i])used.push_back(data->operator[](i));
							else notused.push_back(data->operator[](i));
						}
						double rest=0;
						for(const auto&p:notused)rest+=p.E;
						return point<double>(rest,InvariantMass(used));
					}
				};
				(*im_val)=look_for(0);
				return true;
			}
			<< make_shared<SetOfHists1D>("CentralGammas","InvMass6Gamma",Q_axis_full(res),Axis([im_val]()->double{return im_val->Y();},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("CentralGammas","RestAfter6Gamma",Q_axis_full(res),Axis([im_val]()->double{return im_val->X();},0.0,1.0,1000))
		)
	;
}

void Search2Gamma(Analysis&res){
	auto data=make_shared<vector<particle_kinematics>>();
	auto im_val=make_shared<point<double>>(INFINITY,INFINITY);
	const auto&tn=trigger_gammas_central.number;
	res.Trigger(tn).pre()
	    << [data](){data->clear(); return true;}
	;
	res.Trigger(tn).per_track()
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << [data](WTrack&T)->bool{
			data->push_back({.particle=Particle::gamma(),.E=T.Edep(),.theta=T.Theta(),.phi=T.Phi()});
			return true;
		    }
		)
	;
	res.Trigger(tn).post() 
		<< (make_shared<ChainCheck>()
		    << [data](){return data->size()>=2;}
		    << [data,im_val](){
			SortedPoints<double> table;
			for(size_t i=0;i<data->size();i++)for(size_t j=i+1;j<data->size();j++){
			    double im=InvariantMass({data->operator[](i),data->operator[](j)});
			    double rest=0;
			    for(size_t k=0;k<data->size();k++)if((k!=i)&&(k!=j))
				rest+=data->operator[](k).E;
			    table<<point<double>(rest,im);
			}
			if(table.size()>0)if(table[0].X()<0.050){
			    (*im_val)=table[0];
			    return true;
			}
			return false;
		    }
                    << make_shared<SetOfHists1D>("CentralGammas","InvMass2Gamma",Q_axis_full(res),Axis([im_val]()->double{return im_val->Y();},0.0,1.0,1000))
                    << make_shared<SetOfHists1D>("CentralGammas","RestAfter2Gamma",Q_axis_full(res),Axis([im_val]()->double{return im_val->X();},0.0,1.0,1000))
		)
	;
}
