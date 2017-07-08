// this file is distributed under
// GPL license
#include <fstream>
#include "bound.h"
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
const SortedPoints<> ReadPfFromFile(const string&name){
	SortedPoints<> data;
	ifstream file(name);
	double x,y;
	while(file>>x>>y)data<<point<>(x/1000.,y);
	return data;
}
const pair<Vector4<>,Vector4<>> Compound(
	RANDOM&RG,
	const IFunction<double,RANDOM&>&Pb_distr,
	const IFunction<double,RANDOM&>&Pf_distr,
	const double&s_thr
){
	const auto TotalP=
		Vector4<>::bySpaceC_and_Length4(Vector3<>::basis_z()*Pb_distr(RG),Particle::p().mass())
		+Particle::d().mass();
	while(true){
		const auto he3Pcm=Vector4<>::bySpaceC_and_Length4(
			Vector3<>::RandomIsotropicDirection(RG)*Pf_distr(RG),
			Particle::he3().mass()
		);
		const auto etaPcm=Vector4<>(TotalP.length4())-he3Pcm;
		if(etaPcm.Sqr4()>s_thr){
			return make_pair(etaPcm.Lorentz(-TotalP.Beta()),he3Pcm.Lorentz(-TotalP.Beta()));
		}
	}
}
