// this file is distributed under 
// GPL license
#include <list>
#include <string>
#include <sstream>
#include <random>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include <Experiment/str_get.h>
#include <Experiment/experiment_conv.h>
#include <config.h>
#include "runpluto.h"
using namespace std;
void SimulatePluto(const string&reaction,const string&filename,const double&from,const double&to){
    PUSHD();
    CD(PLUTO);
    std::mt19937 gen;
    std::uniform_int_distribution<int> d(1,254);
    for(size_t runindex=1;runindex<=40;runindex++){
	TF1 *mf=new TF1(CSTR("Uniform"),CSTR("1"),from,to);
	PBeamSmearing *smear = new PBeamSmearing(CSTR("beam_smear"),CSTR("Beam smearing"));
	smear->SetReaction(CSTR("p + d"));
	smear->SetMomentumFunction(mf);
	makeDistributionManager()->Add(smear);
	PUtils::SetSeed(d(gen));
	PReaction my_reaction(CSTR(to_string(to)),CSTR("p"),CSTR("d"),CSTR(reaction),
	    CSTR(filename+"-"+to_string(runindex)),1,0,0,0
	);
	my_reaction.Print();
	my_reaction.Loop(1000000);
    }
    POPD();
}
