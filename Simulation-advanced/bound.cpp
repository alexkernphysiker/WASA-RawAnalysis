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
	while(file>>x>>y)data<<point<>(x/1000.,y);
	return data;
}
const etamesic Compound(
	RANDOM&RG,
	const RandomValueGenerator<>&Pb_distr,
	const RandomValueGenerator<>&Pf_distr,
	const double&s_thr
){
	static PlotDistr1D<>
                PP("bound","P, GeV/c",BinsByCount(40,p_beam_low,p_beam_hi));
	const auto P=Pb_distr(RG);
	PP.Fill(P);
	const auto TotalP=
		Vector4<>::bySpaceC_and_Length4(Vector3<>::basis_z()*P,Particle::p().mass())
		+Particle::d().mass();
	while(true){
		const auto he3Pcm=Vector4<>::bySpaceC_and_Length4(
			Vector3<>::RandomIsotropicDirection(RG)*Pf_distr(RG),
			Particle::he3().mass()
		);
		const auto etaPcm=Vector4<>(TotalP.length4())-he3Pcm;
		if(etaPcm.Sqr4()>s_thr){
			return {.he3=he3Pcm.Lorentz(-TotalP.Beta()),.eta_=etaPcm.Lorentz(-TotalP.Beta())};
		}
	}
}
