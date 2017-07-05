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
	const double mthr=Particle::he3().mass()+Particle::eta().mass();
	const RandomValueTableDistr<double> BW_distr(
		[Bs,mthr](const double&x){return BreitWigner(x,mthr-Bs.first,Bs.second);},
		ChainWithCount(1000,mthr-Bs.second*4.0,mthr)
	);
	return [&RG,BW_distr]()->list<particle_sim>{
		const double mB=BW_distr(RG);
		const double Pb=sqrt(
			pow(
				(
					pow(mB,2)-pow(Particle::p().mass(),2)
					-pow(Particle::d().mass(),2) 
				)
				/(2.0*Particle::d().mass())
			,2)
			-pow(Particle::p().mass(),2)
		);
		const auto TotalP=Vector4<>::bySpaceC_and_Length4(Vector3<>::basis_z()*Pb,Particle::p().mass())+Particle::d().mass();
		//ToDo: Implement fermi P distribution
		const double pfe=0;
		const double m_eta_=sqrt(
			pow(mB,2)+pow(Particle::he3().mass(),2)
			-2.0*mB*sqrt(pow(Particle::he3().mass(),2)+pow(pfe,2))
		);
		const auto etaPcm=Vector4<>::bySpaceC_and_Length4(Vector3<>::RandomIsotropicDirection(RG)*pfe,m_eta_);
		const auto he3Pcm=Vector4<>::bySpaceC_and_Length4(-etaPcm.space_component(),Particle::he3().mass());

		const auto he3Plab=he3Pcm.Lorentz(-TotalP.Beta());
		const auto etaPlab=etaPcm.Lorentz(-TotalP.Beta());

		const auto g1Pcme=Vector4<>::bySpaceC_and_Length4(Vector3<double>::RandomIsotropicDirection(RG)*(etaPlab.length4()/2.0),0.0);
		const auto g2Pcme=Vector4<>::bySpaceC_and_Length4(-g1Pcme.space_component(),0.0);

		const auto g1Plab=g1Pcme.Lorentz(-etaPlab.Beta());
		const auto g2Plab=g2Pcme.Lorentz(-etaPlab.Beta());
		return {
			{.type=Particle::he3() ,.P=he3Plab.space_component()},
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
