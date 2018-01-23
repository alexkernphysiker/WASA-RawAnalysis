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
	Plotter::Instance().SetOutput(".","sim-he3eta-6g");
	const RandomUniform<>Pb_distr(p_beam_low,p_beam_hi);
	Simulate("He3eta-6g",[&Pb_distr]()->list<particle_sim>{
            const auto V0=Direct_eta_production(Pb_distr);
            list<particle_sim> result=ThreePi0Decay(V0.eta_);
            result.push_back({.type=Particle::he3(),.P=V0.he3.P()});
            return result;
        },10);
	return 0;
}
