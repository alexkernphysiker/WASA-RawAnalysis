// this file is distributed under 
// GPL license
#include <list>
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
#include <PPatricte.h>
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
    throw Exception<Particle>("");
}
void Simulate(const std::string&filename,const EventGenerator Generator){
    PUSHD();
    CD(PLUTO);
    std::mt19937 gen;
    std::uniform_int_distribution<int> d(1,254);
    for(ALLMC){
	TFile*f=new TFile(CSTR(filename+"-"+to_string(runindex)),"RECREATE");
	Int_t Npart;
	Float_t Impact;
	Float_t Phi;
	TClonesArray*Particles;
	TTree*T=new TTree("data","");
	T->Branch("Npart",&Npart,"Npart/I");
	T->Branch("Impact",&Impact,"Impact/F");
	T->Branch("Phi",&Phi,"Phi/F");
	T->Branch("Particles",&Particles);
	for(size_t event=0;event<1000000;event++){
	    const auto Result=Generator(gen);
	    Particles=new TClonesArray("PParticle",Result.size());
	    Npart=Result.size();Impact=0;Phi=0;
	    size_t index=0;
	    for(const auto&p:Result){
		new ((*Particles)[0]) PParticle(ParticleName(p.type),p.P.x(),p.P.y(),p.P.z(),p.type.mass());
		index++;
	    }
	    T->Fill();
	    Particles->Clear();
	    delete Particles;
	}
	T->Write();
	f->Close();
	delete T;
	delete f;
    }
    POPD();
}

