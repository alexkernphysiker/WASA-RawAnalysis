// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
#include "bound.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const EventGenerator BoundSimulation2Gamma(RANDOM&RG,const RandomValueGenerator<>&Pb_distr,const RandomValueGenerator<>&Pf_distr){
	return [&RG,&Pb_distr,&Pf_distr]()->list<particle_sim>{
		const auto C=Compound(RG,Pb_distr,Pf_distr);
		const auto&etaPlab=C.eta_;
		const auto&he3Plab=C.he3;
		
		static PlotDistr1D<> pbplot("","P_{beam,lab}, GeV/c",BinsByCount(40,p_beam_low,p_beam_hi));
		static PlotDistr1D<> pfplot("","P_{eta,lab}, GeV/c",BinsByCount(1000,0.0,1.0));
		static PlotDistr1D<> mplot("","m_{eta}, GeV",BinsByCount(200,0.4,0.6));
		pbplot.Fill((etaPlab+he3Plab).space_component().mag());
		pfplot.Fill(etaPlab.space_component().mag());
		mplot.Fill(etaPlab.length4());
		
		const auto gammas_cme=binaryDecay(etaPlab.length4(),0.0,0.0,RandomIsotropicDirection3<>(RG));
		const auto g1Plab=gammas_cme.first.Lorentz(-etaPlab.Beta());
		const auto g2Plab=gammas_cme.second.Lorentz(-etaPlab.Beta());
		return {
			{.type=Particle::gamma(),.P=g1Plab.space_component()},
			{.type=Particle::gamma(),.P=g2Plab.space_component()},
			{.type=Particle::he3() ,.P=he3Plab.space_component()}
		};
	};
}
int main(){
	RANDOM RG;
	Plotter<>::Instance().SetOutput(".","sim-2g");
	const RandomUniform<> P(p_beam_low,p_beam_hi); 
	const auto
	pf1=ReadPfFromFile("distributions/he3eta-pf-75-20.txt"),
	pf2=ReadPfFromFile("distributions/he3eta-pf-80-20.txt"),
	pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
	Plot<>().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3")<<"set title 'read from file'";
	Simulate("bound1-2g",BoundSimulation2Gamma(RG,P,RandomValueTableDistr<>(pf1)));
	Simulate("bound2-2g",BoundSimulation2Gamma(RG,P,RandomValueTableDistr<>(pf2)));
	Simulate("bound3-2g",BoundSimulation2Gamma(RG,P,RandomValueTableDistr<>(pf3)));
	return 0;
}
