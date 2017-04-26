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
#include <Kinematics/reactions.h>
#include "config.h"
using namespace std;
void SimulatePluto(const string&reaction,const string&filename,const double&from,const double&to){
    std::mt19937 gen;
    std::uniform_int_distribution<int> d(1,254);
    for(ALLMC){
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
}
int main(){
    PUSHD();
    CD(PLUTO);
    SimulatePluto("He3 eta","He3eta",1.573,p_beam_hi);
    SimulatePluto("He3 pi0","He3pi0",p_beam_low,p_beam_hi);
    SimulatePluto("He3 pi0 pi0","He3pi0pi0",p_beam_low,p_beam_hi);
    SimulatePluto("He3 pi0 pi0 pi0","He3pi0pi0pi0",p_beam_low,p_beam_hi);
    POPD();
    return 0;
}