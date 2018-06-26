// this file is distributed under 
// GPL license
#include <iostream>
#include <Wasa.hh>
#include <TCutG.h>
#include <math_h/functions.h>
#include <math_h/tabledata.h>
#include <math_h/error.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include <Parameters/parameters.h>
#include <Reconstruction/forward.h>
#include <trackprocessing.h>
#include <detectors.h>
#include <reconstruction.h>
#include <data.h>
#include <montecarlo.h>
#include "analyses.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;
const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
Axis Q_axis_full(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta.P2Q(res.PBeam());},-70.0,30.0,40);}
struct track_info{LorentzVector<>L;double t;};
struct eta_decay_gg{
	track_info A;
	track_info B;
	double expected;
	inline eta_decay_gg&operator=(const eta_decay_gg&src){A=src.A;B=src.B;return *this;}
	inline const double dt()const{return abs(A.t-B.t);}
	inline const double t()const{return (A.t<B.t)?A.t:B.t;}
	inline const LorentzVector<> L()const{return A.L+B.L;}
	inline const double IM()const{return L().M();}
	inline const double diff()const{return abs(IM()-expected);}
	inline const bool operator<(const eta_decay_gg&other)const{return diff()<other.diff();}
	inline const bool operator>(const eta_decay_gg&other)const{return diff()>other.diff();}
};
struct pi0_decay{
        track_info A;
        track_info B;
	inline pi0_decay&operator=(const pi0_decay&src){A=src.A;B=src.B;return *this;}
        inline const double dt()const{return abs(A.t-B.t);}
        inline const double t()const{return (A.t<B.t)?A.t:B.t;}
	inline const LorentzVector<> L()const{return A.L+B.L;}
        inline const double IM()const{return L().M();}
        inline const double diff()const{return pow(IM()-Particle::pi0().mass(),2);}
        inline const bool operator<(const pi0_decay&other)const{return diff()<other.diff();}
        inline const bool operator>(const pi0_decay&other)const{return diff()>other.diff();}
};
struct eta_decay_ppp{
	pi0_decay I;
	pi0_decay J;
	pi0_decay K;
	inline eta_decay_ppp&operator=(const eta_decay_ppp&src){I=src.I;J=src.J;K=src.K;return *this;}
        inline const double t()const{
                const auto t1=(I.t()<J.t())?I.t():J.t();
                return (t1<K.t())?t1:K.t();
        }
	inline const double t2()const{
		const auto t1=(I.t()>J.t())?I.t():J.t();
		return (t1>K.t())?t1:K.t();
	}
	inline const double dt()const{return t2()-t();}
        inline const LorentzVector<> L()const{return I.L()+J.L()+K.L();}
        inline const double IM()const{return L().M();}
        inline const double diff()const{return sqrt(I.diff()+J.diff()+K.diff());}
        inline const bool operator<(const eta_decay_ppp&other)const{return diff()<other.diff();}
        inline const bool operator>(const eta_decay_ppp&other)const{return diff()>other.diff();}
};
void Search3He6Gamma(Analysis&res){
	const auto&tn=trigger_he3_forward.number;
	static auto Ptotal=LorentzVector<>::zero();
	static vector<track_info> gammas;
	static particle_kine he3;
	static track_info He3{.L=LorentzVector<>::zero(),.t=INFINITY};
	res.Trigger(tn).pre()<<(make_shared<ChainOr>()
	    <<make_shared<Hist1D>("He3nCentralGammas6","0-Reference",Q_axis_full(res))
	    << [&res](){
		const auto p=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass());
		const auto d=lorentz_byPM(Zero(),Particle::d().mass());
		Ptotal=p+d;
		gammas.clear();
		He3.t=INFINITY;
		return true;
	    }
	);
	res.Trigger(tn).per_track()<<(make_shared<ChainOr>()
		<<(make_shared<ChainCheck>()
			<<[](WTrack&T){return T.Type()==kFDC;}
			<<ForwardHe3Reconstruction("CentralGammas6andHe3",res,he3)
			<<[](WTrack&T){
				He3={.L=lorentz_byEkM(he3.E,Particle::he3().mass(),direction(he3.phi,he3.theta)),.t=T.Time(kFTH1)};
				return true;
			}
		)
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << make_shared<Hist1D>("He3nCentralGammas6","GammaEnergy",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=getParameter(gamma_E_thr);}
		    << [](WTrack&T)->bool{
			gammas.push_back({.L=lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.),.t=T.Time()});
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas6","GammaEnergyCut",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static eta_decay_ppp six_gamma{.I={.A=He3,.B=He3},.J={.A=He3,.B=He3},.K={.A=He3,.B=He3}};
	Axis 
	he3mmc([](){return (Ptotal-He3.L).M()-Ptotal.M()+Particle::he3().mass()+Particle::eta().mass();},0.4,0.6,200),
	he3me([](){return (Ptotal.E()-He3.L.E());},0.0,1.0,1000),
	ggggggdiff([]()->double{return six_gamma.diff();},0.0,0.2,200),
        he3ggggggimdiff([](){return (He3.L+six_gamma.L()).M()-Ptotal.M();},-0.5,0.5,500),
	ggggggimc([](){return six_gamma.IM()-Ptotal.M()+Particle::he3().mass()+Particle::eta().mass();},0.0,1.0,1000),
        ggggggmm([](){return (Ptotal-six_gamma.L()).M();},0.0,4.0,4000),
        ggggggt([](){return He3.t-six_gamma.t();},-50,100,150),
        ggggggdt([](){return six_gamma.dt();},0,50,50),
	measured_eta_angle([](){return direction(six_gamma.L().P()).th()*180/PI();},0,180,180),
        gamma_gamma_cosi([]()->double{
              const auto d1=direction(six_gamma.I.A.L.P());
              const auto d2=direction(six_gamma.I.B.L.P());
              return ((d1*1.0)*(d2*1.0));
        },-1.0,1.0,100),
        gamma_gamma_cosj([]()->double{
              const auto d1=direction(six_gamma.J.A.L.P());
              const auto d2=direction(six_gamma.J.B.L.P());
              return ((d1*1.0)*(d2*1.0));
        },-1.0,1.0,100),
        gamma_gamma_cosk([]()->double{
              const auto d1=direction(six_gamma.K.A.L.P());
              const auto d2=direction(six_gamma.K.B.L.P());
              return ((d1*1.0)*(d2*1.0));
        },-1.0,1.0,100);
	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas6","Events0",Q_axis_full(res))
	    << make_shared<Hist1D>("He3nCentralGammas6","He3MM0",he3mmc)
	    <<[he3mmc](WTrack&T){return (he3mmc(T)>getParameter(he3mm_cut));}
	    << make_shared<Hist1D>("He3nCentralGammas6","Events1",Q_axis_full(res))
	    << make_shared<Hist1D>("He3nCentralGammas6","He3MM1",he3mmc)
	    << make_shared<Hist1D>("He3nCentralGammas6","GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))

	    <<(make_shared<ChainOr>()		
	        << ( make_shared<ChainCheck>()
			<<[]()->bool{return gammas.size()>=6;}
	                <<[&res]()->bool{
	                        SortedChain<eta_decay_ppp> combinations;
        	                for(size_t i=0;i<(gammas.size()-5);++i)for(size_t j=i+1;j<gammas.size();++j){
	                            for(size_t k=i+1;k<(gammas.size()-3);(++k)+=((k==j)?1:0))for(size_t l=k+1;l<gammas.size();(++l)+=((l==j)?1:0)){
                                        for(size_t o=k+1;o<(gammas.size()-1);(++o)+=(((o==j)||(o==l))?1:0))for(size_t p=o+1;p<gammas.size();(++p)+=(((p==j)||(p==l))?1:0)){
	 				        const auto candidate=eta_decay_ppp{
							.I={.A=gammas[i],.B=gammas[j]},
							.J={.A=gammas[k],.B=gammas[l]},
							.K={.A=gammas[o],.B=gammas[p]}
						};
						vector<track_info> check{
							gammas[i],gammas[j],gammas[k],
							gammas[l],gammas[o],gammas[p]
						};
						//bool passed=true;
						//for(size_t r=0;(r<check.size())&&passed;r++)for(size_t q=r+1;(q<check.size())&&passed;q++){
						//	const auto d1=direction(check[r].L.P());
						//	const auto d2=direction(check[q].L.P());
						//	const auto cosa=(d1*1.0)*(d2*1.0);
						//	if(cosa>0.99)passed=false;
						//}
						if(
						//	passed &&
							(dynamic_cast<const MonteCarlo*>(&res)||(
                                                                (candidate.dt()<getParameter(time_dt))&&
                                                                ((He3.t-candidate.t())>getParameter(time_t1))&&((He3.t-candidate.t())<getParameter(time_t2))
                                                        ))
						)
							combinations<<candidate;
                	        	}
				    }	
	                        }
				if(combinations.size()==0)return false;
	                        six_gamma=combinations[0];
	                        return true;
	                }
			<< make_shared<Hist1D>("He3nCentralGammas6","cosi2",gamma_gamma_cosi)
			<< make_shared<Hist1D>("He3nCentralGammas6","cosj2",gamma_gamma_cosj)
			<< make_shared<Hist1D>("He3nCentralGammas6","cosk2",gamma_gamma_cosk)
			<< make_shared<Hist1D>("He3nCentralGammas6","t2",ggggggt)
			<< make_shared<Hist1D>("He3nCentralGammas6","dt2",ggggggdt)
                        << make_shared<Hist1D>("He3nCentralGammas6","He3MM2",he3mmc)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMM2",ggggggmm)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff2",ggggggdiff)
                        << make_shared<Hist1D>("He3nCentralGammas6","GIM2",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET2",measured_eta_angle)
			<< make_shared<Hist1D>("He3nCentralGammas6","TIM2",he3ggggggimdiff)

			<<[gamma_gamma_cosi,gamma_gamma_cosj,gamma_gamma_cosk](WTrack&T)->bool{
				const double cut=1.0;//I decided to remove this condition
				return (gamma_gamma_cosi(T)<cut)&&(gamma_gamma_cosj(T)<cut)&&(gamma_gamma_cosk(T)<cut);
			}
                        << make_shared<Hist1D>("He3nCentralGammas6","cosi3",gamma_gamma_cosi)
                        << make_shared<Hist1D>("He3nCentralGammas6","cosj3",gamma_gamma_cosj)
                        << make_shared<Hist1D>("He3nCentralGammas6","cosk3",gamma_gamma_cosk)
                        << make_shared<Hist1D>("He3nCentralGammas6","t3",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt3",ggggggdt)
                        << make_shared<Hist1D>("He3nCentralGammas6","He3MM3",he3mmc)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMM3",ggggggmm)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff3",ggggggdiff)
                        << make_shared<Hist1D>("He3nCentralGammas6","GIM3",ggggggimc)
                        << make_shared<Hist1D>("He3nCentralGammas6","ET3",measured_eta_angle)
                        << make_shared<Hist1D>("He3nCentralGammas6","TIM3",he3ggggggimdiff)

			<<[measured_eta_angle](WTrack&T)->bool{return measured_eta_angle(T)<getParameter(eta_theta_thr);}
                        << make_shared<Hist1D>("He3nCentralGammas6","cosi4",gamma_gamma_cosi)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosj4",gamma_gamma_cosj)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosk4",gamma_gamma_cosk)
                        << make_shared<Hist1D>("He3nCentralGammas6","t4",ggggggt)                                                                          
                        << make_shared<Hist1D>("He3nCentralGammas6","dt4",ggggggdt)                                                                        
                        << make_shared<Hist1D>("He3nCentralGammas6","He3MM4",he3mmc)                                                                       
                        << make_shared<Hist1D>("He3nCentralGammas6","GMM4",ggggggmm)                                                                       
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff4",ggggggdiff)                                                                
                        << make_shared<Hist1D>("He3nCentralGammas6","GIM4",ggggggimc)                                                                       
                        << make_shared<Hist1D>("He3nCentralGammas6","ET4",measured_eta_angle)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","TIM4",he3ggggggimdiff)

                        <<[ggggggdiff](WTrack&T)->bool{return ggggggdiff(T)<getParameter(three_pi0);}
                        << make_shared<Hist1D>("He3nCentralGammas6","cosi5",gamma_gamma_cosi)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosj5",gamma_gamma_cosj)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosk5",gamma_gamma_cosk)
                        << make_shared<Hist1D>("He3nCentralGammas6","t5",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt5",ggggggdt)
			<< make_shared<Hist1D>("He3nCentralGammas6","He3MM5",he3mmc)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMM5",ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff5",ggggggdiff)
			<< make_shared<Hist1D>("He3nCentralGammas6","GIM5",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET5",measured_eta_angle)
			<< make_shared<Hist1D>("He3nCentralGammas6","TIM5",he3ggggggimdiff)

                        <<[ggggggmm](WTrack&T)->bool{return (ggggggmm(T)>getParameter(gamma_mm_lo))&&(ggggggmm(T)<getParameter(gamma_mm_hi));}
                        << make_shared<Hist1D>("He3nCentralGammas6","cosi6",gamma_gamma_cosi)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosj6",gamma_gamma_cosj)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosk6",gamma_gamma_cosk)
                        << make_shared<Hist1D>("He3nCentralGammas6","t6",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt6",ggggggdt)
			<< make_shared<Hist1D>("He3nCentralGammas6","He3MM6",he3mmc)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMM6",ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff6",ggggggdiff)
			<< make_shared<Hist1D>("He3nCentralGammas6","GIM6",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET6",measured_eta_angle)
			<< make_shared<Hist1D>("He3nCentralGammas6","TIM6",he3ggggggimdiff)

                        <<[ggggggimc](WTrack&T)->bool{return (ggggggimc(T)>getParameter(gamma_im_lo6))&&(ggggggimc(T)<getParameter(gamma_im_hi));}
                        << make_shared<Hist1D>("He3nCentralGammas6","cosi7",gamma_gamma_cosi)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosj7",gamma_gamma_cosj)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","cosk7",gamma_gamma_cosk)
                        << make_shared<Hist1D>("He3nCentralGammas6","Events7",Q_axis_full(res))
                        << make_shared<Hist1D>("He3nCentralGammas6","t7",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt7",ggggggdt)
			<< make_shared<Hist1D>("He3nCentralGammas6","He3MM7",he3mmc)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMM7",ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff7",ggggggdiff)
			<< make_shared<Hist1D>("He3nCentralGammas6","GIM7",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET7",measured_eta_angle)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM7",Q_axis_full(res),he3ggggggimdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3ME7",Q_axis_full(res),he3me)
	        )
	    )
	);
}
void Search3He2Gamma(Analysis&res){
	const auto&tn=trigger_he3_forward.number;
	static auto Ptotal=LorentzVector<>::zero();
	static vector<track_info> gammas;
	static particle_kine he3;
	static track_info He3{.L=LorentzVector<>::zero(),.t=INFINITY};
	res.Trigger(tn).pre()<<(make_shared<ChainOr>()
	    <<make_shared<Hist1D>("He3nCentralGammas2","0-Reference",Q_axis_full(res))
	    << [&res](){
		const auto p=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass());
		const auto d=lorentz_byPM(Zero(),Particle::d().mass());
		Ptotal=p+d;
		gammas.clear();
		He3.t=INFINITY;
		return true;
	    }
	);
	res.Trigger(tn).per_track()<<(make_shared<ChainOr>()
		<<(make_shared<ChainCheck>()
			<<[](WTrack&T){return T.Type()==kFDC;}
			<<ForwardHe3Reconstruction("CentralGammas2andHe3",res,he3)
			<<[](WTrack&T){
				He3={.L=lorentz_byEkM(he3.E,Particle::he3().mass(),direction(he3.phi,he3.theta)),.t=T.Time(kFTH1)};
				return true;
			}
		)
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << make_shared<Hist1D>("He3nCentralGammas2","GammaEnergy",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=getParameter(gamma_E_thr);}
		    << [](WTrack&T)->bool{
			gammas.push_back({.L=lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.),.t=T.Time()});
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas2","GammaEnergyCut",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static eta_decay_gg two_gamma{.A=He3,.B=He3};
	static eta_decay_ppp six_gamma{.I={.A=He3,.B=He3},.J={.A=He3,.B=He3},.K={.A=He3,.B=He3}};
	Axis 
	he3mmc([](){return (Ptotal-He3.L).M()-Ptotal.M()+Particle::he3().mass()+Particle::eta().mass();},0.4,0.6,200),
	he3mm([](){return (Ptotal-He3.L).M();},0.4,0.6,200),
	he3me([](){return (Ptotal.E()-He3.L.E());},0.0,1.0,1000),
	he3ggimdiff([](){return (He3.L+two_gamma.L()).M()-Ptotal.M();},-0.5,0.5,500),
	ggim([](){return two_gamma.IM();},0.0,1.0,1000),
	ggimc([](){return two_gamma.IM()-Ptotal.M()+Particle::he3().mass()+Particle::eta().mass();},0.0,1.0,1000),
	ggmm([](){return (Ptotal-two_gamma.L()).M();},0.0,4.0,4000),
        ggt([](){return He3.t-two_gamma.t();},-50,100,150),
	ggdt([](){return two_gamma.dt();},0,50,50),
	measured_eta_angle([](){return direction(two_gamma.L().P()).th()*180/PI();},0,180,180),
        gamma_gamma_cosa([]()->double{
              const auto d1=direction(two_gamma.A.L.P());
              const auto d2=direction(two_gamma.B.L.P());
              return ((d1*1.0)*(d2*1.0));
        },-1.0,1.0,100);

	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas2","Events0",Q_axis_full(res))
	    << make_shared<Hist1D>("He3nCentralGammas2","He3MM0",he3mmc)
	    << make_shared<Hist1D>("He3nCentralGammas2","old_He3MM0",he3mm)
	    <<(make_shared<ChainOr>()
                << ( make_shared<ChainCheck>()                                                                                                             
                    <<[he3mmc](WTrack&T){return (he3mmc(T)>getParameter(he3mm_cut));}
                    <<(make_shared<ChainOr>()
                        << ( make_shared<ChainCheck>()
                                << []()->bool{return gammas.size()>=2;}
				<< (make_shared<ChainOr>()
				    <<(make_shared<ChainCheck>()
	                                << [&res]()->bool{
                                               SortedChain<eta_decay_gg> pairs;
                                 	       for(size_t i=0;i<gammas.size();i++)for(size_t j=i+1;j<gammas.size();j++){
	                                                const auto candidate=eta_decay_gg{
								.A=gammas[i],.B=gammas[j],
								.expected=Ptotal.M()-Particle::he3().mass()
							};
                                                	if((dynamic_cast<const MonteCarlo*>(&res)||(
                                        	           (candidate.dt()<getParameter(time_dt))&&
                                	                   ((He3.t-candidate.t())>getParameter(time_t1))&&((He3.t-candidate.t())<getParameter(time_t2))
                        	                        )))
                	                                        pairs<<candidate;
        	                                }
	                                        if(pairs.size()==0)return false;
                                        	two_gamma=pairs[0];
                                	        return true;
                        	        }
                	                << make_shared<Hist1D>("He3nCentralGammas2","GGcos0",gamma_gamma_cosa)
        	                        << make_shared<Hist1D>("He3nCentralGammas2","t0",ggt)
	                                << make_shared<Hist1D>("He3nCentralGammas2","dt0",ggdt)
				    )
                                    <<(make_shared<ChainCheck>()
                                        << [&res]()->bool{
                                               SortedChain<eta_decay_gg> pairs;
                                               for(size_t i=0;i<gammas.size();i++)for(size_t j=i+1;j<gammas.size();j++){
                                                    const auto candidate=eta_decay_gg{
							.A=gammas[i],.B=gammas[j],
							.expected=Ptotal.M()-Particle::he3().mass()
						    };
                                                    pairs<<candidate;
                                               }
                                               if(pairs.size()==0)return false;
                                               two_gamma=pairs[0];
                                               return true;
                                        }                                                                                                                  
                                        << make_shared<Hist1D>("He3nCentralGammas2","GGcos00",gamma_gamma_cosa)
                                        << make_shared<Hist1D>("He3nCentralGammas2","t00",ggt)
                                        << make_shared<Hist1D>("He3nCentralGammas2","dt00",ggdt)
                                    )                                                                                                                 
			      )
			)
		    )
		)
		<< ( make_shared<ChainCheck>()
		    <<[he3mmc](WTrack&T){return (he3mmc(T)>getParameter(he3mm_cut));}
		    << make_shared<Hist1D>("He3nCentralGammas2","Events1",Q_axis_full(res))
		    << make_shared<Hist1D>("He3nCentralGammas2","He3MM1",he3mmc)
		    << make_shared<Hist1D>("He3nCentralGammas2","GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))
		    <<(make_shared<ChainOr>()
			<< ( make_shared<ChainCheck>()
				<< []()->bool{return gammas.size()>=2;}
				<< [&res]()->bool{
					SortedChain<eta_decay_gg> pairs;
					for(size_t i=0;i<gammas.size();i++)for(size_t j=i+1;j<gammas.size();j++){
						const auto candidate=eta_decay_gg{
							.A=gammas[i],.B=gammas[j],
							.expected=Ptotal.M()-Particle::he3().mass()
						};
						const auto d1=direction(candidate.A.L.P());
						const auto d2=direction(candidate.B.L.P());
						const auto cosa=((d1*1.0)*(d2*1.0));
						if(
							(cosa<getParameter(gg_theta_cut))&&
							(dynamic_cast<const MonteCarlo*>(&res)||(
								(candidate.dt()<getParameter(time_dt))&&
								((He3.t-candidate.t())>getParameter(time_t1))&&((He3.t-candidate.t())<getParameter(time_t2))
							))
						)
							pairs<<candidate;
					}
					if(pairs.size()==0)return false;
					two_gamma=pairs[0];
					return true;
				}

				<< make_shared<Hist1D>("He3nCentralGammas2","GGcos2",gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","t2",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","dt2",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","He3MM2",he3mmc)
				<< make_shared<Hist1D>("He3nCentralGammas2","GMM2",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","GIM2",ggimc)
				<< make_shared<Hist1D>("He3nCentralGammas2","TIM2",he3ggimdiff)
				<< make_shared<Hist1D>("He3nCentralGammas2","ET2",measured_eta_angle)

				<<[gamma_gamma_cosa](WTrack&T)->bool{return gamma_gamma_cosa(T)<0;}//moved to pair selection
				<< make_shared<Hist1D>("He3nCentralGammas2","GGcos3",gamma_gamma_cosa)
                                << make_shared<Hist1D>("He3nCentralGammas2","t3",ggt)                                                                      
                                << make_shared<Hist1D>("He3nCentralGammas2","dt3",ggdt)                                                                    
                                << make_shared<Hist1D>("He3nCentralGammas2","He3MM3",he3mmc)                                                               
                                << make_shared<Hist1D>("He3nCentralGammas2","GMM3",ggmm)                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","GIM3",ggimc)                                                                  
                                << make_shared<Hist1D>("He3nCentralGammas2","TIM3",he3ggimdiff)                                                            
                                << make_shared<Hist1D>("He3nCentralGammas2","ET3",measured_eta_angle)

				<<[measured_eta_angle](WTrack&T)->bool{return measured_eta_angle(T)<getParameter(eta_theta_thr);}
				<< make_shared<Hist1D>("He3nCentralGammas2","GGcos4",gamma_gamma_cosa)
                                << make_shared<Hist1D>("He3nCentralGammas2","t4",ggt)                                                                      
                                << make_shared<Hist1D>("He3nCentralGammas2","dt4",ggdt)                                                                    
                                << make_shared<Hist1D>("He3nCentralGammas2","He3MM4",he3mmc)                                                               
                                << make_shared<Hist1D>("He3nCentralGammas2","GMM4",ggmm)                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","GIM4",ggimc)                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","TIM4",he3ggimdiff)                                                            
                                << make_shared<Hist1D>("He3nCentralGammas2","ET4",measured_eta_angle)

				<<[ggmm](WTrack&T)->bool{return (ggmm(T)>getParameter(gamma_mm_lo))&&(ggmm(T)<getParameter(gamma_mm_hi));}
				<< make_shared<Hist1D>("He3nCentralGammas2","GGcos5",gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","t5",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","dt5",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","He3MM5",he3mmc)
				<< make_shared<Hist1D>("He3nCentralGammas2","GMM5",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","GIM5",ggimc)
				<< make_shared<Hist1D>("He3nCentralGammas2","TIM5",he3ggimdiff)
				<< make_shared<Hist1D>("He3nCentralGammas2","ET5",measured_eta_angle)

				<<[ggimc](WTrack&T)->bool{return (ggimc(T)>getParameter(gamma_im_lo))&&(ggimc(T)<getParameter(gamma_im_hi));}
				<< make_shared<Hist1D>("He3nCentralGammas2","GGcos6",gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","Events6",Q_axis_full(res))                                 
				<< make_shared<Hist1D>("He3nCentralGammas2","t6",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","dt6",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","He3MM6",he3mmc)
				<< make_shared<Hist1D>("He3nCentralGammas2","GMM6",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","GIM6",ggimc)
				<< make_shared<SetOfHists1D>("He3nCentralGammas2","TIM6",Q_axis_full(res),he3ggimdiff)
				<< make_shared<Hist1D>("He3nCentralGammas2","ET6",measured_eta_angle)
				<< make_shared<SetOfHists1D>("He3nCentralGammas2","He3ME6",Q_axis_full(res),he3me)
			)
		    )
		)

	    )
	);
}
