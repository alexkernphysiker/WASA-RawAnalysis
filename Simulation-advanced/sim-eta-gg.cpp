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
	Plotter::Instance().SetOutput(".","sim-he3eta-gg");
	const RandomUniform<>Pb_distr(p_beam_low,p_beam_hi);
	Simulate("He3eta-gg",[&RG,&Pb_distr]()->list<particle_sim>{
            const auto V0=Direct_eta_production(RG,Pb_distr);
	    const auto V1=binaryDecay(V0.eta_.M(),0.,0.,randomIsotropic<3>(RG));
	    const auto g1=V1.first.Transform(-V0.eta_.Beta());
	    const auto g2=V1.second.Transform(-V0.eta_.Beta());
            return {
                   {.type=Particle::he3(),.P=V0.he3.P()},
                   {.type=Particle::gamma(),.P=g1.P()},
                   {.type=Particle::gamma(),.P=g2.P()}
            };
        },10);
	return 0;
}
