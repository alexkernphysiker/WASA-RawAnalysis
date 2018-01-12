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
            const auto V0=Direct_eta_production(RG,Pb_distr);
            list<particle_sim> result=ThreePi0Decay(RG,V0.eta_);
            result.push_back({.type=Particle::he3(),.P=V0.he3.P()});
            return result;
        },10);
	return 0;
}
