// this file is distributed under 
// GPL license
#include <list>
#include <string>
#include <sstream>
#include <random>
#include <Experiment/str_get.h>
#include <Experiment/experiment_conv.h>
#include <config.h>
#include <TFile.h>
#include "runsim.h"
using namespace std;
void Simulate(const std::string&filename,const EventGenerator gen){
    PUSHD();
    CD(PLUTO);
    std::mt19937 gen;
    std::uniform_int_distribution<int> d(1,254);
    for(ALLMC){
	TFile* f=new TFile(CSTR(filename+"-"+to_string(runindex)),"RECREATE");
    }
    POPD();
}

