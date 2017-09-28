// this file is distributed under
// GPL license
#include <fstream>
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "bound.h"
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const LinearInterpolation<> ReadPfFromFile(const string&name){
	SortedPoints<> data;
	ifstream file(name);
	double x,y;
	while(file>>x>>y)data<<make_point(x/1000.,y);
	return data;
}
const etamesic Compound(
	RANDOM&RG,
	const RandomValueGenerator<>&Pb_distr,
	const RandomValueGenerator<>&Pf_distr,
	const double&s_thr
){
	const auto P=Pb_distr(RG);
	const auto TotalP=lorentz_byPM(Z<>()*P,Particle::p().mass())+lorentz_byPM(Zero<>(),Particle::d().mass());
	while(true){
		const auto he3Pcm=lorentz_byPM(randomIsotropic<3>(RG)*Pf_distr(RG),Particle::he3().mass());
		const auto etaPcm=lorentz_byPM(Zero<>(),TotalP.M())-he3Pcm;
		if(etaPcm.M_sqr()>s_thr){
			return {.he3=he3Pcm.Transform(-TotalP.Beta()),.eta_=etaPcm.Transform(-TotalP.Beta())};
		}
	}
}
