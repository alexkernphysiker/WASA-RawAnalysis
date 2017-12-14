// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
#include "bound.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-he3eta-6g");
	const RandomUniform<>Pb_distr(1.573,p_beam_hi);
	const RandomValueTableDistr<> THETA=LinearInterpolation<>([](double t)->double{
		return sin(t)*(3.+cos(t));
	},ChainWithStep(0.0,0.001,PI()));
	Simulate("He3eta-6g",[&RG,&Pb_distr,&THETA]()->list<particle_sim>{
		const auto pb=Pb_distr(RG);
		const auto p0=lorentz_byPM(Z()*pb,Particle::p().mass());
    		const auto d0=lorentz_byPM(Zero(),Particle::d().mass());
           	const auto total=p0+d0;
	   	const static RandomUniform<> PHI(-PI(),PI());
            	const auto V0=binaryDecay(total.M(),Particle::he3().mass(),Particle::eta().mass(),direction(PHI(RG),THETA(RG)));
            	const auto he3=V0.first.Transform(-total.Beta());
            	const auto eta=V0.second.Transform(-total.Beta());
		list<particle_sim> result=ThreePi0Decay(RG,eta);
		result.push_back({.type=Particle::he3(),.P=he3.P()});
		return result;
        },10);
	return 0;
}
