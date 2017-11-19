// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-he3eta-gg");
	const RandomUniform<>Pb_distr(1.573,p_beam_hi);
	const RandomValueTableDistr<> THETA=LinearInterpolation<>([](double t)->double{
		return sin(t)*(3.+cos(t));
	},ChainWithStep(0.0,0.001,PI()));
	Simulate("He3eta_gg",[&RG,&Pb_distr,&THETA]()->list<particle_sim>{
	    const auto pb=Pb_distr(RG);
            const auto p0=lorentz_byPM(Z()*pb,Particle::p().mass());
            const auto d0=lorentz_byPM(Zero(),Particle::d().mass());
            const auto total=p0+d0;
	    const static RandomUniform<> PHI(-PI(),PI());
            const auto V0=binaryDecay(total.M(),Particle::he3().mass(),Particle::eta().mass(),direction(PHI(RG),THETA(RG)));
            const auto he3=V0.first.Transform(-total.Beta());
            const auto eta=V0.second.Transform(-total.Beta());
	    const auto V1=binaryDecay(eta.M(),0.,0.,randomIsotropic<3>(RG));
	    const auto g1=V1.first.Transform(-eta.Beta());
	    const auto g2=V1.second.Transform(-eta.Beta());
            return {
                   {.type=Particle::he3(),.P=he3.P()},
                   {.type=Particle::gamma(),.P=g1.P()},
                   {.type=Particle::gamma(),.P=g2.P()}
            };
        },10);
	return 0;
}
