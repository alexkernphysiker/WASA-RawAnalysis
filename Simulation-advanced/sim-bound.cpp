// this file is distributed under
// GPL license
#include <math_h/vectors.h>
#include <math_h/randomfunc.h>
#include <math_h/functions.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
const EventGenerator BoundSimulation2Gamma(mt19937&RG,const pair<double,double>&Bs){
	const RandomValueTableDistr<double> BW_distr(
		[Bs](const double&x){return BreitWigner(x,Bs.first,Bs.second);},
		ChainWithCount(1000,0.0,(Particle::he3().mass()+Particle::eta().mass())*2.0)
	);
	const double mp=Particle::p().mass(),md=Particle::d().mass();
	return [&RG,BW_distr,mp,md](){
		const double InvMass=BW_distr(RG);
		const double Pb=sqrt(((pow(InvMass,2)-pow(mp,2)-pow(md,2))/(2.0*md))-pow(mp,2));
		const auto TotalP=Vector4<double>::SpaceLength4(Vector3<double>::basis_z()*Pb,mp)+
			Vector4<double>::SpaceLength4(Vector3<double>::zero(),md);
		//ToDo: implement Fermi momentum distribution
		const auto etaPcm=Vector4<double>::SpaceLength4(Vector3<double>::RandomIsotropicDirection(RG)*0.1,Particle::eta().mass());
		
		list<particle_sim> res;
		return res;
	};
}
int main(){
	mt19937 RG;
	Simulate("bound-2g-04",BoundSimulation2Gamma(RG,make_pair(0.00402,0.0156/2.0)));
	return 0;
}
