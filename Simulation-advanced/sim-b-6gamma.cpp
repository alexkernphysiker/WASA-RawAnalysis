// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
#include "bound.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const EventGenerator BoundSimulation6Gamma(RANDOM&RG,const RandomValueGenerator<>&Pb_distr,const RandomValueGenerator<>&Pf_distr){
	return [&RG,&Pb_distr,&Pf_distr]()->list<particle_sim>{
		static PlotDistr1D<> pbplot("he3eta","P_{beam,lab}, GeV/c",BinsByCount(40,p_beam_low,p_beam_hi));
		static PlotDistr1D<> pbplot2("he36gamma","P_{beam,lab},Gev/c",BinsByCount(40,p_beam_low,p_beam_hi));
		while(true){
			const auto C=Compound(RG,Pb_distr,Pf_distr,pow(3.0*Particle::pi0().mass(),2));
			if(C.eta_.M()<3.0*Particle::pi0().mass())continue;
			auto output=ThreePi0Decay(RG,C.eta_);
			output.push_back({.type=Particle::he3(),.P=C.he3.P()});
			pbplot.Fill((C.he3+C.eta_).P().M());
			auto P=LorentzVector<>::zero();
			for(const auto&p:output){
				P+=lorentz_byPM(p.P,p.type.mass());
			}
			pbplot2.Fill(P.P().M());
			return output;
		};
	};
}
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-6gamma");
	const RandomValueTableDistr<>P=Points<>{{p_beam_low,2.0},{p_beam_hi,1.0}}; 
	const auto
	pf1=ReadPfFromFile("distributions/he3eta-pf-75-20.txt"),
	pf2=ReadPfFromFile("distributions/he3eta-pf-80-20.txt"),
	pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
	Plot().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3")<<"set title 'read from file'";
	Simulate("bound1-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf1)),10);
	Simulate("bound2-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf2)),10);
	Simulate("bound3-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf3)),10);
	return 0;
}
