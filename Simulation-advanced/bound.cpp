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
LinearInterpolation<> ReadPfFromFile(const string&name){
	SortedPoints<> data;
	ifstream file(name);
	double x,y;
	while(file>>x>>y)data<<make_point(x/1000.,y);
	return data;
}
etamesic Compound(
	const RandomValueGenerator<>&Pb_distr,
	const RandomValueGenerator<>&Pf_distr,
	const double&s_thr
){
	const auto P=Pb_distr();
	const auto TotalP=lorentz_byPM(Z<>()*P,Particle::p().mass())+lorentz_byPM(Zero<>(),Particle::d().mass());
	while(true){
		const auto he3Pcm=lorentz_byPM(randomIsotropic<3>()*Pf_distr(),Particle::he3().mass());
		const auto etaPcm=lorentz_byPM(Zero<>(),TotalP.M())-he3Pcm;
		if(etaPcm.M_sqr()>s_thr){
			return {.he3=he3Pcm.Transform(-TotalP.Beta()),.eta_=etaPcm.Transform(-TotalP.Beta())};
		}
	}
}
etamesic Direct_eta_production(
	const RandomValueGenerator<>&Pb_distr
){
	while(true){
		const auto P=Pb_distr();
		const auto TotalP=lorentz_byPM(Z<>()*P,Particle::p().mass())+lorentz_byPM(Zero<>(),Particle::d().mass());
		if(TotalP.M()>(Particle::he3().mass()+Particle::eta().mass())){
		    static auto filltable=[](){
			static auto f=[](const double&costh,const double&a,const double&b,const double&g){
				return 1.0 + a*costh + b*costh*costh + g*costh*costh*costh;
			};
			static const LinearInterpolation<>alpha=Points<>{
                            {-0.001,0},{0.04,0},
                            {0.112,0.52},{0.130,0.62},{0.145,0.74},{0.170,1.2}
                    	};
                    	static const LinearInterpolation<>beta=Points<>{
                            {-0.001,-0.32},{0.04,-0.32},
                            {0.112,-0.33},{0.130,-0.34},{0.145,-0.30},{0.170,-0.34}
        		};
			static const LinearInterpolation<>gamma=Points<>{
                            {-0.001,0},{0.04,0},
                            {0.112,-0.10},{0.130,-0.18},{0.145,-0.21},{0.170,-0.29}
			};
			vector<RandomValueTableDistr<>> res;
			for(size_t i=0;i<170;i++){
				const double p=0.001*i;
				res.push_back(RandomValueTableDistr<>(
					[&p](const double&theta){
						const auto ct=cos(theta);
						const auto rho=f(ct,alpha(p),beta(p),gamma(p))*sin(theta);
						return (rho>0)?rho:0.0;
					},ChainWithStep(0.0,0.001,PI())
				));
			}
			return res;
		    };
		    const static auto THETA=filltable();
		    const static RandomUniform<> PHI(-PI(),PI());
                    const auto V0=binaryDecay(TotalP.M(),Particle::eta().mass(),Particle::he3().mass(),direction(0.0,0.0));
		    const double theta=THETA[size_t(V0.first.P().M()/0.001)]();
		    const double phi=PHI();
                    const auto eta= V0.first.Rotate(direction(Y()),theta).Rotate(direction(Z()),phi).Transform(-TotalP.Beta());
                    const auto he3=V0.second.Rotate(direction(Y()),theta).Rotate(direction(Z()),phi).Transform(-TotalP.Beta());
                    return {.he3=he3,.eta_=eta};
		}
	}
}
double lambda(const double&x,const double&y,const double&z){
	return x*x+y*y+z*z-2.*x*y-2.*y*z-2.*z*x;
}
list<particle_sim> ThreePi0Decay(
	const LorentzVector<>&eta
){
	while(true){
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
		mplot.Fill(eta.M());
		const double M=eta.M(),s=M*M,m=Particle::pi0().mass();
		const double smin=pow(m,2),smax=pow(M-m,2);
		const RandomUniform<> s_distr(smin,smax);
		const double s1=s_distr();
		const double s2=s_distr();
		const double s3=s+3.*m*m-s1-s2;
		if(!isfinite(s3))continue;
		if((s3<=smin)||(s3>=smax))continue;
		const double E2=(M*M+m*m-s2)/(2.*M),E3=(M*M+m*m-s3)/(2.*M),E1=(M*M+m*m-s1)/(2.*M);
		if((E1<m)||(E2<m)||(E3<m))continue;
		const double p1=sqrt(E1*E1-m*m),p2=sqrt(E2*E2-m*m),p3=sqrt(E3*E3-m*m);
		if((p1<0)||(p2<0)||(p3<0))throw Exception<EventGenerator,2>("Momenta error");
		const double p2x=(E3*E3-E2*E2-p1*p1)/(2.*p1);
		const double p2y=sqrt(p2*p2-p2x*p2x);
		const double p3x=-(p1+p2x);
		const auto 
		P1=lorentz_byPM(vec(p1,0.,0.),m),
		P2=lorentz_byPM(vec(p2x,p2y,0.),m),
		P3=lorentz_byPM(vec(p3x,-p2y,0.),m);
		if(!isfinite(P1.M()))continue;
		if(!isfinite(P2.M()))continue;
		if(!isfinite(P3.M()))continue;
		s1plot.Fill(s1);s2plot.Fill(s2);s3plot.Fill(s3);
		s12plot.Fill(s1,s2);s13plot.Fill(s1,s3);s23plot.Fill(s2,s3);
		p1plot.Fill(p1);p2plot.Fill(p2);p3plot.Fill(p3);
		const RandomUniform<> TH(-PI<>(),PI<>());
		const auto 
		RR=randomIsotropic<3>().Rotations()*Rotation(direction(P1.P()),TH());
		plplot.Fill(((RR*P1.P())^(RR*P2.P()))*(RR*P3.P()));
		const vector<LorentzVector<>> pizeros={
			lorentz_byPM(RR*P1.P(),P1.M()).Transform(-eta.Beta()),
			lorentz_byPM(RR*P2.P(),P2.M()).Transform(-eta.Beta()),
			lorentz_byPM(RR*P3.P(),P3.M()).Transform(-eta.Beta()),
		};
		im1plot.Fill((pizeros[0]+pizeros[1]+pizeros[2]).M());
		auto G=LorentzVector<>::zero();
		list<particle_sim> output;
		for(const auto PiP:pizeros){
			const auto g1Pcm=lorentz_byPM(randomIsotropic<3>()*m/2.,0.);
			const auto g2Pcm=lorentz_byPM(-g1Pcm.P(),0.);
			const auto piframe=PiP.Beta();
			if(!isfinite(piframe.M()))throw Exception<EventGenerator,106>("Invalid pi0 frame");
			output.push_back({.type=Particle::gamma(),.P=g1Pcm.Transform(-piframe).P()});
			output.push_back({.type=Particle::gamma(),.P=g2Pcm.Transform(-piframe).P()});
			G+=g1Pcm.Transform(-piframe)+g2Pcm.Transform(-piframe);
		}
		im2plot.Fill(G.M());
		return output;
	}
}
