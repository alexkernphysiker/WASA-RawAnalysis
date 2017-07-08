// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
#include "bound.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const EventGenerator BoundSimulation6Gamma(
	RANDOM&RG,const IFunction<double,RANDOM&>&Pb_distr,
	const IFunction<double,RANDOM&>&Pf_distr
){
	return [&RG,&Pb_distr,&Pf_distr]()->list<particle_sim>{
		const auto C=Compound(RG,Pb_distr,Pf_distr,pow(3.0*Particle::pi0().mass(),2));
		const auto&etaPlab=C.first;
		const auto&he3Plab=C.second;
		static PlotDistr1D<> mplot("6g","m_{eta'}, GeV",BinsByCount(1000,0.0,1.0));
		mplot.Fill(etaPlab.length4());
		const static RandomUniform<> im_pair(
			2.0*Particle::pi0().mass(),
			etaPlab.length4()-Particle::pi0().mass()
		);
		double im12=0,im23=0;
		return{
			{.type=Particle::he3() ,.P=he3Plab.space_component()}
		};
	};
}
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-6gamma");
	const RandomUniform<> P(p_beam_low,p_beam_hi); 
	const auto
	pf1=ReadPfFromFile("distributions/he3eta-pf-75-20.txt"),
	pf2=ReadPfFromFile("distributions/he3eta-pf-80-20.txt"),
	pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
	Plot<>().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3")<<"set title 'read from file'";
	Simulate("bound1-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf1)));
	return 0;
}
