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
#include "central.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

const Reaction& He3eta(){
    static Reaction res(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
    return res;
}
Axis Q_axis_full(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta().P2Q(res.PBeam());},-70.0,30.0,40);}


void SearchGammaTracks(Analysis&res){
	auto data=make_shared<vector<particle_kinematics>>();
	//ToDo: move this trigger number to experiment_conv.h
	res.Trigger(0).pre()
	    << make_shared<Hist1D>("CentralGammas","0-Reference",Q_axis_full(res))
	    << [data](){data->clear(); return true;};
	res.Trigger(0).per_track()<<(make_shared<ChainCheck>()
	    << [](WTrack&T)->bool{
		return T.Type()==kCDN;
	    }
	    << [data](WTrack&T)->bool{
		data->push_back({.particle=Particle::gamma(),.E=T.Edep(),.theta=T.Theta(),.phi=T.Phi()});
		return true;
	    }
	);
	auto im_val=make_shared<double>(INFINITY);
	res.Trigger(0).post() << (make_shared<ChainCheck>()
	    << [data,im_val](){
		(*im_val)=INFINITY;
		SortedPoints<double> table;
		for(size_t i=0;i<data->size();i++)for(size_t j=i+1;j<data->size();j++){
		    double im=InvariantMass({data->operator[](i),data->operator[](j)});
		    double rest=0;
		    for(size_t k=0;k<data->size();k++)if((k!=i)&&(k!=j))
			rest+=data->operator[](k).E;
		    table<<point<double>(rest,im);
		}
		if(table.size()>0)if(table[0].X()<0.050){
		    (*im_val)=table[0].Y();
		    return true;
		}
		return false;
	    }
	    << make_shared<SetOfHists1D>("CentralGammas","InvMass2Gamma",Q_axis_full(res),Axis([im_val]()->double{return *im_val;},0.0,0.8,800))
	)
	<< make_shared<Hist1D>("CentralGammas","neutral_tracks_count",Axis([data]()->double{return data->size();},-0.5,9.5,10));
    }
