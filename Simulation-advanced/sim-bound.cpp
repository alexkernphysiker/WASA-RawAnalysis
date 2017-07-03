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
	const double mp=Particle::p().mass(),md=Particle::d().mass(),
		mhe=Particle::he3().mass(),meta=Particle::eta().mass();
	const RandomValueTableDistr<double> BW_distr(
		[Bs,meta,mhe](const double&x){return BreitWigner(x,mhe+meta-Bs.first,Bs.second);},
		ChainWithCount(1000,0.0,(mhe+meta)*2.0)
	);
	return [&RG,BW_distr,mp,md,mhe]()->list<particle_sim>{
		const double InvMass=BW_distr(RG);
		const double Pb=sqrt(((pow(InvMass,2)-pow(mp,2)-pow(md,2))/(2.0*md))-pow(mp,2));
		const auto TotalP=Vector4<double>::SpaceLength4(Vector3<double>::basis_z()*Pb,mp)+
			Vector4<double>::SpaceLength4(Vector3<double>::zero(),md);
		//ToDo: implement Fermi momentum distribution
		const double pfe=0.1;
		const auto he3Pcm=Vector4<double>::SpaceLength4(Vector3<double>::RandomIsotropicDirection(RG)*pfe,mhe);
		const double meta_=sqrt(pow(InvMass,2)+pow(mhe,2)-2.0*InvMass*sqrt(pow(mhe,2)+pow(pfe,2)));
		const auto etaPcm=Vector4<double>::SpaceLength4(-he3Pcm.space_component(),meta_);
		const auto g1Pcm2=Vector4<double>::SpaceLength4(Vector3<double>::RandomIsotropicDirection(RG)*meta_/2.0,0.0);
		const auto g2Pcm2=Vector4<double>::SpaceLength4(-g1Pcm2.space_component(),0.0);
		const auto he3Plab=he3Pcm.Lorentz(-TotalP.Beta());
		const auto g1Plab=g1Pcm2.Lorentz(-etaPcm.Beta()).Lorentz(-TotalP.Beta());
		const auto g2Plab=g2Pcm2.Lorentz(-etaPcm.Beta()).Lorentz(-TotalP.Beta());
		return {
			{.type=Particle::he3(),.P=he3Plab.space_component()},
			{.type=Particle::gamma(),.P=g1Plab.space_component()},
			{.type=Particle::gamma(),.P=g2Plab.space_component()}
		};
	};
}
int main(){
	mt19937 RG;
	Simulate("bound2g-04",BoundSimulation2Gamma(RG,make_pair(0.00402,0.0156/2.0)));
	return 0;
}
