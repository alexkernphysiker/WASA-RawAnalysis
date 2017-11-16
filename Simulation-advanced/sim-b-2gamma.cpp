// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
#include "bound.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const EventGenerator BoundSimulation2Gamma(RANDOM&RG,const RandomValueGenerator<>&Pb_distr,const RandomValueGenerator<>&Pf_distr,vector<PlotDistr1D<>>&plots){
	return [&RG,&Pb_distr,&Pf_distr,&plots]()->list<particle_sim>{
		const auto C=Compound(RG,Pb_distr,Pf_distr);
		const auto&etaPlab=C.eta_;
		const auto&he3Plab=C.he3;
		plots[0].Fill((etaPlab+he3Plab).P().M());
		plots[1].Fill(etaPlab.P().M());
		plots[2].Fill(etaPlab.M());
		const auto gammas_cme=binaryDecay(etaPlab.M(),0.0,0.0,randomIsotropic<3>(RG));
		const auto g1Plab=gammas_cme.first.Transform(-etaPlab.Beta());
		const auto g2Plab=gammas_cme.second.Transform(-etaPlab.Beta());
		return {
			{.type=Particle::gamma(),.P=g1Plab.P()},
			{.type=Particle::gamma(),.P=g2Plab.P()},
			{.type=Particle::he3() ,.P=he3Plab.P()}
		};
	};
}
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-2g");
	const RandomUniform<> P(p_beam_low,p_beam_hi); 
	const auto
	pf1=ReadPfFromFile("distributions/he3eta-pf-75-20.txt"),
	pf2=ReadPfFromFile("distributions/he3eta-pf-80-20.txt"),
	pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
	Plot().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3")<<"set title 'read from file'";
	vector<PlotDistr1D<>>
	plots1{
		PlotDistr1D<>("1","P_{beam,lab},GeV/C",BinsByCount(40,p_beam_low,p_beam_hi)),
		PlotDistr1D<>("1","P_{eta,lab},GeV/C",BinsByCount(1000,0.0,1.0)),
		PlotDistr1D<>("1","m_{eta_b},GeV",BinsByCount(200,0.4,0.6))
	},
	plots2{
		PlotDistr1D<>("2","P_{beam,lab},GeV/C",BinsByCount(40,p_beam_low,p_beam_hi)),
		PlotDistr1D<>("2","P_{eta,lab},GeV/C",BinsByCount(1000,0.0,1.0)),
		PlotDistr1D<>("2","m_{eta_b},GeV",BinsByCount(200,0.4,0.6))
	},
	plots3{
		PlotDistr1D<>("3","P_{beam,lab},GeV/C",BinsByCount(40,p_beam_low,p_beam_hi)),
		PlotDistr1D<>("3","P_{eta,lab},GeV/C",BinsByCount(1000,0.0,1.0)),
		PlotDistr1D<>("3","m_{eta_b},GeV",BinsByCount(200,0.4,0.6))
	};
	Simulate("bound1-2g",BoundSimulation2Gamma(RG,P,RandomValueTableDistr<>(pf1),plots1),10);
	Simulate("bound2-2g",BoundSimulation2Gamma(RG,P,RandomValueTableDistr<>(pf2),plots2),10);
	Simulate("bound3-2g",BoundSimulation2Gamma(RG,P,RandomValueTableDistr<>(pf3),plots3),10);
	return 0;
}
