// this file is distributed under 
// GPL license
#include <list>
#include <iostream>
#include <string>
#include <sstream>
#include <random>
#include <math_h/error.h>
#include <Experiment/str_get.h>
#include <Experiment/experiment_conv.h>
#include <config.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <PParticle.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
const string ParticleName(const Particle&p){
    if(p==Particle::gamma())return "g";
    if(p==Particle::n())return "n";
    if(p==Particle::p())return "p";
    if(p==Particle::d())return "d";
    if(p==Particle::he3())return "He3";
    if(p==Particle::he4())return "He4";
    if(p==Particle::eta())return "eta";
    if(p==Particle::pi0())return "pi0";
    if(p==Particle::pi_plus())return "pi+";
    if(p==Particle::pi_minus())return "pi-";
    throw Exception<Particle>("Unknown particle type: cannot find it's name");
}
void Simulate(const std::string&filename,const EventGenerator Generator){
    PUSHD();
    CD(PLUTO);
    for(size_t runindex=1;runindex<=4;runindex++){
	cerr<<"Running simulation number "<<runindex<<" started"<<endl;
	TFile*f=new TFile(CSTR(filename+"-"+to_string(runindex)+".root"),"RECREATE");
	auto Result=Generator();
	Int_t Npart=Result.size();
	Float_t Impact=0;
	Float_t Phi=0;
	TClonesArray*Particles=new TClonesArray("PParticle",Result.size());
	TTree*T=new TTree("data","");
	T->Branch("Npart",&Npart,"Npart/I");
	T->Branch("Impact",&Impact,"Impact/F");
	T->Branch("Phi",&Phi,"Phi/F");
	T->Branch("Particles",&Particles);
	for(size_t event=1;event<=10000000;event++){
	    if((event%10000)==0)cerr<<event<<" events"<<endl;
	    while((Result=Generator()).size()==0);
	    Npart=Result.size();
	    Impact=0;
	    Phi=0;
	    Particles->Clear();
	    size_t index=0;
	    for(const auto&p:Result){
		new ((*Particles)[index]) PParticle(ParticleName(p.type).c_str(),p.P.x(),p.P.y(),p.P.z(),p.type.mass());
		index++;
	    }
	    T->Fill();
	}
	T->Write();
	f->Close();
	cerr<<"done."<<endl;
    }
    POPD();
}
