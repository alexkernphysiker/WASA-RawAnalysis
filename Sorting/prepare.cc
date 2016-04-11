// this file is distributed under 
// MIT license
#include "prepare.h"
#include "data.h"
#include "montecarlo.h"
namespace ReactionSetup{
	using namespace TrackAnalyse;
	const Reaction& He3eta(){
		static Reaction res(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
		return res;
	}
	Axis Q_axis(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta().P2Q(res.PBeam());},0.0,30.0,12);}
	Analysis*Prepare(const AnalysisModification mode){
		switch(mode){
			case forData:
				return new RealData();
			case forMC:
				return new MonteCarlo();
		};
		throw;
	}
}