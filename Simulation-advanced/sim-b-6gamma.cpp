// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
#include "bound.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const double lambda(const double&x,const double&y,const double&z){
	return x*x+y*y+z*z-2.*x*y-2.*y*z-2.*z*x;
}
const EventGenerator BoundSimulation6Gamma(
	RANDOM&RG,const IFunction<double,RANDOM&>&Pb_distr,
	const IFunction<double,RANDOM&>&Pf_distr
){
	return [&RG,&Pb_distr,&Pf_distr]()->list<particle_sim>{
		const auto C=Compound(RG,Pb_distr,Pf_distr,pow(3.0*Particle::pi0().mass(),2));
		const auto&etaPlab=C.first;
		const auto&he3Plab=C.second;
		static PlotDistr1D<> 
		mplot("6g","m_{eta}, GeV",BinsByCount(200,0.4,0.6)),
		s12plot("6g","s12, GeV^2",BinsByCount(200,0.0,0.2)),
		s13plot("6g","s13, GeV^2",BinsByCount(200,0.0,0.2)),
		s23plot("6g","s23, GeV^2",BinsByCount(200,0.0,0.2)),
		p1plot("6g","p1, GeV/c",BinsByCount(200,0.0,0.2)),
		p2plot("6g","p2, GeV/c",BinsByCount(200,0.0,0.2)),
		p3plot("6g","p3, GeV/c",BinsByCount(200,0.0,0.2));
		mplot.Fill(etaPlab.length4());
		const double M=etaPlab.length4(),s=M*M,m=Particle::pi0().mass();
		const double smin=pow(2.0*m,2),smax=pow(M-m,2);
		const RandomUniform<> s_distr(smin,smax);
		double s12,s13,s23;
		for(bool done=false;!done;){
			s12=s_distr(RG);
			s13=s_distr(RG);
			if(isfinite(s12)&&isfinite(s13)){
				if((s12>=smin)&&(s12<=smax)&&(s13>=smin)&&(s13<=smax)){
					s23=s+3.*m*m-s12-s13;
					done=(s23<=smax)&&(s23>=smin);
				}
			}
		}
		s12plot.Fill(s12);s13plot.Fill(s13);s23plot.Fill(s23);
		const double E2=(M*M+m*m-s13)/(2.*M),E3=(M*M+m*m-s12)/(2.*M),E1=(M*M+m*m-s23)/(2.*M);
		if((E1<m)||(E2<m)||(E3<m))throw Exception<EventGenerator,1>("Energies error");
		const double p1=sqrt(E1*E1-m*m),p2=sqrt(E2*E2-m*m),p3=sqrt(E3*E3-m*m);
		if((p1<0)||(p2<0)||(p3<0))throw Exception<EventGenerator,2>("Momenta error");
		const double p2x=(E3*E3-E2*E2-p1*p1)/(2.*p1);
		const double p2y=sqrt(p2*p2-p2x*p2x);
		const double p3x=-(p1+p2x);
		const LorentzVector<Vector2<>> P1(E1,{p1,0.}),P2(E2,{p2x,p2y}),P3(E3,{p3x,-p2y});
		p1plot.Fill(p1);p2plot.Fill(p2);p3plot.Fill(p3);
		static RandomUniform<> decayorientation(0,2.0*PI());
		const auto DecayPlane=Plane3D<>::ByNormalVectorAndTheta(
			Vector3<>::RandomIsotropicDirection(RG),
			decayorientation(RG)
		);
		const list<LorentzVector<>> pizeros={
			DecayPlane(P1).Lorentz(-etaPlab.Beta()),
			DecayPlane(P2).Lorentz(-etaPlab.Beta()),
			DecayPlane(P3).Lorentz(-etaPlab.Beta())
		};
		list<particle_sim> output;
		for(const auto PiP:pizeros){
			const auto g1Pcm=LorentzVector<>(m/2.,Vector3<>::RandomIsotropicDirection(RG)*m/2.);
			const auto g2Pcm=LorentzVector<>(m/2.,-g1Pcm.space_component());
			output.push_back({.type=Particle::gamma(),.P=g1Pcm.Lorentz(-PiP.Beta()).space_component()});
			output.push_back({.type=Particle::gamma(),.P=g2Pcm.Lorentz(-PiP.Beta()).space_component()});
		}
		output.push_back({.type=Particle::he3() ,.P=he3Plab.space_component()});
		return output;
	};
}
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-6gamma");
	const RandomUniform<> P(p_beam_low,p_beam_hi); 
	const auto
	pf1=ReadPfFromFile("distributions/he3eta-pf-75-20.txt"),
	pf2=ReadPfFromFile("distributions/he3eta-pf-80-20.txt"),
	pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
	Plot<>().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3")<<"set title 'read from file'";
	Simulate("bound1-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf1)));
	Simulate("bound2-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf2)));
	Simulate("bound3-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf3)));
	return 0;
}
