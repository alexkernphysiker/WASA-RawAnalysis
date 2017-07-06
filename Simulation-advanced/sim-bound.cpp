// this file is distributed under
// GPL license
#include <fstream>
#include <math_h/vectors.h>
#include <math_h/randomfunc.h>
#include <math_h/functions.h>
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const SortedPoints<> ReadFromFile(const string&name){
	SortedPoints<> data;
	ifstream file(name);
	double x,y;
	while(file>>x>>y)data<<point<>(x/1000.,y);
	return data;
}
const pair<Vector4<>,Vector4<>> Compound(RANDOM&RG,const RandomValueTableDistr<>&BW_distr,const RandomValueTableDistr<>&PF_distr){
	const double mB=BW_distr(RG);
	const double Pb=sqrt(pow((pow(mB,2)-pow(Particle::p().mass(),2)-pow(Particle::d().mass(),2))/(2.0*Particle::d().mass()),2)-pow(Particle::p().mass(),2));
	const auto TotalP=Vector4<>::bySpaceC_and_Length4(Vector3<>::basis_z()*Pb,Particle::p().mass())+Particle::d().mass();
	const double pfe=PF_distr(RG);
	const double m_eta_=sqrt(pow(mB,2)+pow(Particle::he3().mass(),2)-2.0*mB*sqrt(pow(Particle::he3().mass(),2)+pow(pfe,2)));
	const auto etaPcm=Vector4<>::bySpaceC_and_Length4(Vector3<>::RandomIsotropicDirection(RG)*pfe,m_eta_);
	const auto he3Pcm=Vector4<>::bySpaceC_and_Length4(-etaPcm.space_component(),Particle::he3().mass());
	return make_pair(etaPcm.Lorentz(-TotalP.Beta()),he3Pcm.Lorentz(-TotalP.Beta()));
}
const EventGenerator BoundSimulation2Gamma(RANDOM&RG,const RandomValueTableDistr<>&BW_distr,const RandomValueTableDistr<>&PF_distr){
	return [&RG,BW_distr,PF_distr]()->list<particle_sim>{
		const auto C=Compound(RG,BW_distr,PF_distr);
		const auto&etaPlab=C.first;
		const auto&he3Plab=C.second;
		const auto g1Pcme=Vector4<>::bySpaceC_and_Length4(Vector3<double>::RandomIsotropicDirection(RG)*(etaPlab.length4()/2.0),0.0);
		const auto g2Pcme=Vector4<>::bySpaceC_and_Length4(-g1Pcme.space_component(),0.0);
		const auto g1Plab=g1Pcme.Lorentz(-etaPlab.Beta());
		const auto g2Plab=g2Pcme.Lorentz(-etaPlab.Beta());
		return {
			{.type=Particle::gamma(),.P=g1Plab.space_component()},
			{.type=Particle::gamma(),.P=g2Plab.space_component()},
			{.type=Particle::he3() ,.P=he3Plab.space_component()}
		};
	};
}
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-bound");
	const double mthr=Particle::he3().mass()+Particle::eta().mass();
	const SortedPoints<> 
	bw1([mthr](const double&x){return BreitWigner(x,mthr-0.00402,0.01560/2.0);},
                ChainWithCount(1000,mthr-0.07,mthr+0.03)),
	bw2([mthr](const double&x){return BreitWigner(x,mthr-0.00619,0.01739/2.0);},
		ChainWithCount(1000,mthr-0.07,mthr+0.03)),
	bw3([mthr](const double&x){return BreitWigner(x,mthr-0.01110,0.02059/2.0);},
		ChainWithCount(1000,mthr-0.07,mthr+0.03));
	Plot<>().Line(bw1,"1").Line(bw2,"2").Line(bw3,"3");

	const auto
	pf1=ReadFromFile("distributions/he3eta-pf-75-20.txt").XRange(0.0,0.5),
	pf2=ReadFromFile("distributions/he3eta-pf-80-20.txt").XRange(0.0,0.5),
	pf3=ReadFromFile("distributions/he3eta-pf-90-20.txt").XRange(0.0,0.5);
	Plot<>().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3");
	
	Simulate("bound1-2g",BoundSimulation2Gamma(RG,bw1,pf1));
	Simulate("bound2-2g",BoundSimulation2Gamma(RG,bw2,pf2));
	Simulate("bound3-2g",BoundSimulation2Gamma(RG,bw3,pf3));
	return 0;
}
