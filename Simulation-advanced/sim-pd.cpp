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
	Plotter::Instance().SetOutput(".","sim-pd");
	const RandomUniform<>Pb_distr(p_beam_low,p_beam_hi);
	const RandomValueTableDistr<> T=LinearInterpolation<>([](double t)->double{
		return exp(12.6932-33.2654*t+66.2139*t*t-68.9458*t*t*t);
	},ChainWithStep(0.0,0.001,0.45));
	vector<LinearInterpolation<>> theta_cm_table;
	for(double pb=p_beam_low;pb<=p_beam_hi;pb+=0.001){
        	SortedPoints<> out;
	        for(double theta_cm=0;theta_cm<PI();theta_cm+=0.001){
        	    const auto p0=lorentz_byPM(x()*pb,Particle::p().mass());
	            const auto d0=lorentz_byPM(zero(),Particle::d().mass());
        	    const auto total=p0+d0;
	            const auto final_cm=binaryDecay(total.M(),Particle::p().mass(),Particle::d().mass(),direction(theta_cm));
        	    const auto beta=-total.Beta();
	            const auto p1=final_cm.first.Transform(beta);
        	    const auto d1=final_cm.second.Transform(beta);
	            out<<make_point(d1.P().M(),theta_cm);
		}
		theta_cm_table.push_back(out);
        }
	const auto THETA=[&theta_cm_table](const double&Pb,const double&t)->double{
		int index=int((Pb-p_beam_low)/0.001);
		if(index<0)
			index=0;
		if(size_t(index)>=theta_cm_table.size())
			index=theta_cm_table.size()-1;
		return theta_cm_table[index](t);
	};
	Simulate("pd_",[&RG,&Pb_distr,&T,&THETA]()->list<particle_sim>{
	    const auto pb=Pb_distr(RG);
            const auto p0=lorentz_byPM(Z()*pb,Particle::p().mass());
            const auto d0=lorentz_byPM(Zero(),Particle::d().mass());
            const auto total=p0+d0;
	    const static RandomUniform<> PHI(-PI(),PI());
            const auto final_cm=binaryDecay(total.M(),Particle::p().mass(),Particle::d().mass(),direction(PHI(RG),THETA(pb,T(RG))));
            const auto beta=-total.Beta();
            const auto p1=final_cm.first.Transform(beta);
            const auto d1=final_cm.second.Transform(beta);
            return {
                        {.type=Particle::p(),.P=p1.P()},
                        {.type=Particle::d(),.P=d1.P()}
            };
        },10);
	return 0;
}
