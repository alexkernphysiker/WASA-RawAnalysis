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
const EventGenerator BoundSimulation6Gamma(RANDOM&RG,const RandomValueGenerator<>&Pb_distr,const RandomValueGenerator<>&Pf_distr){
	const auto generator=[&RG](const etamesic&C)->list<particle_sim>{
		const auto&etaPlab=C.eta_;
		const auto&he3Plab=C.he3;
		static PlotDistr1D<> 
		mplot("6g","m_{eta}, GeV",BinsByCount(100,0.3,0.7)),
		im1plot("6g","IM before, GeV",BinsByCount(100,0.3,0.7)),
		im2plot("6g","IM after, GeV",BinsByCount(100,0.3,0.7)),
		s1plot("6g","s1, GeV^2",BinsByCount(100,0.0,0.4)),
		s2plot("6g","s2, GeV^2",BinsByCount(100,0.0,0.4)),
		s3plot("6g","s3, GeV^2",BinsByCount(100,0.0,0.4)),
		p1plot("6g","p1, GeV/c",BinsByCount(100,0.0,0.4)),
		p2plot("6g","p2, GeV/c",BinsByCount(100,0.0,0.4)),
		p3plot("6g","p3, GeV/c",BinsByCount(100,0.0,0.4)),
		plplot("pl","",BinsByCount(100,-1.0,1.0));
		static PlotDistr2D<>
		s12plot("s1 vs s2",BinsByCount(100,0.0,0.4),BinsByCount(100,0.0,0.4)),
		s13plot("s1 vs s3",BinsByCount(100,0.0,0.4),BinsByCount(100,0.0,0.4)),
		s23plot("s2 vs s3",BinsByCount(100,0.0,0.4),BinsByCount(100,0.0,0.4));
		mplot.Fill(etaPlab.M());
		const double M=etaPlab.M(),s=M*M,m=Particle::pi0().mass();
		const double smin=pow(m,2),smax=pow(M-m,2);
		const RandomUniform<> s_distr(smin,smax);
		const double s1=s_distr(RG);
		const double s2=s_distr(RG);
		const double s3=s+3.*m*m-s1-s2;
		if(!isfinite(s3))return {};
		if((s3<=smin)||(s3>=smax))return {};
		const double E2=(M*M+m*m-s2)/(2.*M),E3=(M*M+m*m-s3)/(2.*M),E1=(M*M+m*m-s1)/(2.*M);
		if((E1<m)||(E2<m)||(E3<m))return {};
		const double p1=sqrt(E1*E1-m*m),p2=sqrt(E2*E2-m*m),p3=sqrt(E3*E3-m*m);
		if((p1<0)||(p2<0)||(p3<0))throw Exception<EventGenerator,2>("Momenta error");
		const double p2x=(E3*E3-E2*E2-p1*p1)/(2.*p1);
		const double p2y=sqrt(p2*p2-p2x*p2x);
		const double p3x=-(p1+p2x);
		const auto 
		P1=lorentz_byPM(desCartes(p1,0.,0.),m),
		P2=lorentz_byPM(desCartes(p2x,p2y,0.),m),
		P3=lorentz_byPM(desCartes(p3x,-p2y,0.),m);
		if(!isfinite(P1.M()))return {};
		if(!isfinite(P2.M()))return {};
		if(!isfinite(P3.M()))return {};
		s1plot.Fill(s1);s2plot.Fill(s2);s3plot.Fill(s3);
		s12plot.Fill(s1,s2);s13plot.Fill(s1,s3);s23plot.Fill(s2,s3);
		p1plot.Fill(p1);p2plot.Fill(p2);p3plot.Fill(p3);
		const RandomUniform<> TH(-PI<>(),PI<>());
		const auto 
		RR=randomIsotropic<3>(RG).Rotations()*Rotation(direction(P1.S()),TH(RG));
		plplot.Fill(((RR*P1.S())^(RR*P2.S()))*(RR*P3.S()));
		const vector<LorentzVector<>> pizeros={
			lorentz_byPM(RR*P1.S(),P1.M()).Transform(-etaPlab.Beta()),
			lorentz_byPM(RR*P2.S(),P2.M()).Transform(-etaPlab.Beta()),
			lorentz_byPM(RR*P3.S(),P3.M()).Transform(-etaPlab.Beta()),
		};
		im1plot.Fill((pizeros[0]+pizeros[1]+pizeros[2]).M());
		auto G=LorentzVector<>::zero();
		list<particle_sim> output;
		for(const auto PiP:pizeros){
			const auto g1Pcm=lorentz_byPM(randomIsotropic<3>(RG)*m/2.,0.);
			const auto g2Pcm=lorentz_byPM(-g1Pcm.S(),0.);
			const auto piframe=PiP.Beta();
			if(!isfinite(piframe.M()))throw Exception<EventGenerator,106>("Invalid pi0 frame");
			output.push_back({.type=Particle::gamma(),.P=g1Pcm.Transform(-piframe).S()});
			output.push_back({.type=Particle::gamma(),.P=g2Pcm.Transform(-piframe).S()});
			G+=g1Pcm.Transform(-piframe)+g2Pcm.Transform(-piframe);
		}
		im2plot.Fill(G.M());
		output.push_back({.type=Particle::he3() ,.P=he3Plab.S()});
		return output;
	};
	return [generator,&RG,&Pb_distr,&Pf_distr]()->list<particle_sim>{
		static PlotDistr1D<> pbplot("he3eta","P_{beam,lab}, GeV/c",BinsByCount(40,p_beam_low,p_beam_hi));
		static PlotDistr1D<> pbplot2("he36gamma","P_{beam,lab},Gev/c",BinsByCount(40,p_beam_low,p_beam_hi));
		while(true){
			const auto C=Compound(RG,Pb_distr,Pf_distr,pow(3.0*Particle::pi0().mass(),2));
			const auto res=generator(C);
			if(res.size()>0){
				pbplot.Fill((C.he3+C.eta_).S().M());
				auto P=LorentzVector<>::zero();
				for(const auto&p:res){
					P+=lorentz_byPM(p.P,p.type.mass());
				}
				pbplot2.Fill(P.S().M());
				return res;
			}
		};
	};
}
int main(){
	RANDOM RG;
	Plotter::Instance().SetOutput(".","sim-6gamma");
	const RandomValueTableDistr<>P=Points<>{{p_beam_low,3.0},{p_beam_hi,1.0}}; 
	const auto
	pf1=ReadPfFromFile("distributions/he3eta-pf-75-20.txt"),
	pf2=ReadPfFromFile("distributions/he3eta-pf-80-20.txt"),
	pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
	Plot().Line(pf1,"1").Line(pf2,"2").Line(pf3,"3")<<"set title 'read from file'";
	Simulate("bound1-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf1)));
	Simulate("bound2-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf2)));
	Simulate("bound3-6g",BoundSimulation6Gamma(RG,P,RandomValueTableDistr<>(pf3)));
	return 0;
}
