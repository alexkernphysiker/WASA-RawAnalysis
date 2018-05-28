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
	inline eta_decay_gg&operator=(const eta_decay_gg&src){A=src.A;B=src.B;return *this;}
	inline const double dt()const{return abs(A.t-B.t);}
	inline const double t()const{return (A.t<B.t)?A.t:B.t;}
	inline const LorentzVector<> L()const{return A.L+B.L;}
	inline const double IM()const{return L().M();}
	inline const double diff()const{return abs(IM()-Particle::eta().mass());}
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

		    << [](WTrack&T)->bool{return T.Edep()>=0.01;}
		    << [](WTrack&T)->bool{
			gammas.push_back({.L=lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.),.t=T.Time()});
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas6","GammaEnergyCut",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static eta_decay_ppp six_gamma{.I={.A=He3,.B=He3},.J={.A=He3,.B=He3},.K={.A=He3,.B=He3}};
	Axis 
	he3mmc([](){return (Ptotal-He3.L).M()-Ptotal.M()+He3.L.M()+Particle::eta().mass();},0.4,0.6,200),
	he3me([](){return (Ptotal.E()-He3.L.E());},0.0,4.0,4000),
	ggggggdiff([]()->double{return six_gamma.diff();},0.0,0.2,200),
        he3ggggggimdiff([](){return (He3.L+six_gamma.L()).M()-Ptotal.M();},-0.5,0.5,500),
	ggggggimc([](){return six_gamma.IM()-Ptotal.M()+He3.L.M()+Particle::eta().mass();},0.0,1.0,1000),
        ggggggmm([](){return (Ptotal-six_gamma.L()).M();},0.0,4.0,4000),
        ggggggt([](){return He3.t-six_gamma.t();},-50,100,150),
        ggggggdt([](){return six_gamma.dt();},0,50,50),
	measured_eta_angle([](){return direction(six_gamma.L().P()).th()*180/PI();},0,180,180);

	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas6","Events0",Q_axis_full(res))
	    << make_shared<Hist1D>("He3nCentralGammas6","He3MM0",he3mmc)
	    <<[he3mmc](WTrack&T){return (he3mmc(T)>0.51);}
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
						bool passed=true;
						for(size_t r=0;(r<check.size())&&passed;r++)for(size_t q=r+1;(q<check.size())&&passed;q++){
							const auto d1=direction(check[r].L.P());
							const auto d2=direction(check[q].L.P());
							const auto sina=((d1*1.0)^(d2*1.0)).M();
							//const auto cosa=(d1*1.0)*(d2*1.0);
							if((sina<0.15))passed=false;
						}
						if(
							passed &&
							(dynamic_cast<const MonteCarlo*>(&res)||(
                                                                (candidate.dt()<20.)&&
                                                                ((He3.t-candidate.t())>-5.)&&((He3.t-candidate.t())<35.)
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
			<< make_shared<Hist1D>("He3nCentralGammas6","t2",ggggggt)
			<< make_shared<Hist1D>("He3nCentralGammas6","dt2",ggggggdt)
                        << make_shared<Hist1D>("He3nCentralGammas6","He3MM2",he3mmc)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMM2",ggggggmm)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff2",ggggggdiff)
                        << make_shared<Hist1D>("He3nCentralGammas6","GIM2",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET2",measured_eta_angle)
			<< make_shared<Hist1D>("He3nCentralGammas6","TIM2",he3ggggggimdiff)

			<<[measured_eta_angle](WTrack&T)->bool{return measured_eta_angle(T)<40;}
                        << make_shared<Hist1D>("He3nCentralGammas6","t3",ggggggt)                                                                          
                        << make_shared<Hist1D>("He3nCentralGammas6","dt3",ggggggdt)                                                                        
                        << make_shared<Hist1D>("He3nCentralGammas6","He3MM3",he3mmc)                                                                       
                        << make_shared<Hist1D>("He3nCentralGammas6","GMM3",ggggggmm)                                                                       
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff3",ggggggdiff)                                                                
                        << make_shared<Hist1D>("He3nCentralGammas6","GIM3",ggggggimc)                                                                       
                        << make_shared<Hist1D>("He3nCentralGammas6","ET3",measured_eta_angle)                                                              
                        << make_shared<Hist1D>("He3nCentralGammas6","TIM3",he3ggggggimdiff)

                        <<[ggggggdiff](WTrack&T)->bool{return ggggggdiff(T)<0.12;}
                        << make_shared<Hist1D>("He3nCentralGammas6","t4",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt4",ggggggdt)
			<< make_shared<Hist1D>("He3nCentralGammas6","He3MM4",he3mmc)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMM4",ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff4",ggggggdiff)
			<< make_shared<Hist1D>("He3nCentralGammas6","GIM4",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET4",measured_eta_angle)
			<< make_shared<Hist1D>("He3nCentralGammas6","TIM4",he3ggggggimdiff)

                        <<[ggggggmm](WTrack&T)->bool{return (ggggggmm(T)>2.65)&&(ggggggmm(T)<3.00);}
                        << make_shared<Hist1D>("He3nCentralGammas6","t5",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt5",ggggggdt)
			<< make_shared<Hist1D>("He3nCentralGammas6","He3MM5",he3mmc)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMM5",ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff5",ggggggdiff)
			<< make_shared<Hist1D>("He3nCentralGammas6","GIM5",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET5",measured_eta_angle)
			<< make_shared<Hist1D>("He3nCentralGammas6","TIM5",he3ggggggimdiff)

                        <<[ggggggimc](WTrack&T)->bool{return (ggggggimc(T)>0.35)&&(ggggggimc(T)<0.65);}
                        << make_shared<Hist1D>("He3nCentralGammas6","Events6",Q_axis_full(res))
                        << make_shared<Hist1D>("He3nCentralGammas6","t6",ggggggt)
                        << make_shared<Hist1D>("He3nCentralGammas6","dt6",ggggggdt)
			<< make_shared<Hist1D>("He3nCentralGammas6","He3MM6",he3mmc)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMM6",ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff6",ggggggdiff)
			<< make_shared<Hist1D>("He3nCentralGammas6","GIM6",ggggggimc)
			<< make_shared<Hist1D>("He3nCentralGammas6","ET6",measured_eta_angle)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM6",Q_axis_full(res),he3ggggggimdiff)

			<< make_shared<Hist2D>("He3nCentralGammas6","He3MME6",he3mmc,he3me)
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

		    << [](WTrack&T)->bool{return T.Edep()>=0.01;}
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
	he3mmc([](){return (Ptotal-He3.L).M()-Ptotal.M()+He3.L.M()+Particle::eta().mass();},0.4,0.6,200),
	he3mm([](){return (Ptotal-He3.L).M();},0.4,0.6,200),
	he3me([](){return (Ptotal.E()-He3.L.E());},0.0,8.0,8000),
	he3ggimdiff([](){return (He3.L+two_gamma.L()).M()-Ptotal.M();},-0.5,0.5,500),
	ggim([](){return two_gamma.IM();},0.0,1.0,1000),
	ggimc([](){return two_gamma.IM()-Ptotal.M()+He3.L.M()+Particle::eta().mass();},0.0,1.0,1000),
	ggmm([](){return (Ptotal-two_gamma.L()).M();},0.0,4.0,4000),
        ggt([](){return He3.t-two_gamma.t();},-50,100,150),
	ggdt([](){return two_gamma.dt();},0,50,50),
	measured_eta_angle([](){return direction(two_gamma.L().P()).th()*180/PI();},0,180,180),
	gamma_gamma_sina([]()->double{
              const auto d1=direction(two_gamma.A.L.P());
              const auto d2=direction(two_gamma.B.L.P());
              return ((d1*1.0)^(d2*1.0)).M();
        },0.0,1.0,100),
        gamma_gamma_cosa([]()->double{
              const auto d1=direction(two_gamma.A.L.P());
              const auto d2=direction(two_gamma.B.L.P());
              return ((d1*1.0)*(d2*1.0));
        },-1.0,1.0,200);

	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas2","Events0",Q_axis_full(res))
	    << make_shared<Hist1D>("He3nCentralGammas2","He3MM0",he3mmc)
	    << make_shared<Hist1D>("He3nCentralGammas2","old_He3MM0",he3mm)
	    <<(make_shared<ChainOr>()

		<< ( make_shared<ChainCheck>()
		    <<[he3mmc](WTrack&T){return (he3mmc(T)>0.51);}
		    << make_shared<Hist1D>("He3nCentralGammas2","Events1",Q_axis_full(res))
		    << make_shared<Hist1D>("He3nCentralGammas2","He3MM1",he3mmc)
		    << make_shared<Hist1D>("He3nCentralGammas2","GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))
		    <<(make_shared<ChainOr>()
			<< ( make_shared<ChainCheck>()
				<< []()->bool{return gammas.size()>=2;}
				<< [&res]()->bool{
					SortedChain<eta_decay_gg> pairs;
					for(size_t i=0;i<gammas.size();i++)for(size_t j=i+1;j<gammas.size();j++){
						const auto candidate=eta_decay_gg{.A=gammas[i],.B=gammas[j]};
						const auto d1=direction(candidate.A.L.P());
						const auto d2=direction(candidate.B.L.P());
						const auto sina=((d1*1.0)^(d2*1.0)).M();
						//const auto cosa=((d1*1.0)*(d2*1.0));
						if(
							(sina>0.15)&&
							(dynamic_cast<const MonteCarlo*>(&res)||(
								(candidate.dt()<20.)&&
								((He3.t-candidate.t())>-5.)&&((He3.t-candidate.t())<35.)
							))
						)
							pairs<<candidate;
					}
					if(pairs.size()==0)return false;
					two_gamma=pairs[0];
					return true;
				}

				<< make_shared<Hist2D>("He3nCentralGammas2","GGangle2",gamma_gamma_sina,gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","GGsin2",gamma_gamma_sina)
				<< make_shared<Hist1D>("He3nCentralGammas2","Events2",Q_axis_full(res))
				<< make_shared<Hist1D>("He3nCentralGammas2","t2",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","dt2",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","He3MM2",he3mmc)
				<< make_shared<Hist1D>("He3nCentralGammas2","GMM2",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","GIM2",ggimc)
				<< make_shared<Hist1D>("He3nCentralGammas2","TIM2",he3ggimdiff)
				<< make_shared<Hist1D>("He3nCentralGammas2","ET2",measured_eta_angle)

				<<[gamma_gamma_sina](WTrack&T)->bool{return gamma_gamma_sina(T)<0.9;}
                                << make_shared<Hist2D>("He3nCentralGammas2","GGangle3",gamma_gamma_sina,gamma_gamma_cosa)                                  
				<< make_shared<Hist1D>("He3nCentralGammas2","GGsin3",gamma_gamma_sina)
                                << make_shared<Hist1D>("He3nCentralGammas2","Events3",Q_axis_full(res))                                                    
                                << make_shared<Hist1D>("He3nCentralGammas2","t3",ggt)                                                                      
                                << make_shared<Hist1D>("He3nCentralGammas2","dt3",ggdt)                                                                    
                                << make_shared<Hist1D>("He3nCentralGammas2","He3MM3",he3mmc)                                                               
                                << make_shared<Hist1D>("He3nCentralGammas2","GMM3",ggmm)                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","GIM3",ggimc)                                                                  
                                << make_shared<Hist1D>("He3nCentralGammas2","TIM3",he3ggimdiff)                                                            
                                << make_shared<Hist1D>("He3nCentralGammas2","ET3",measured_eta_angle)

				<<[measured_eta_angle](WTrack&T)->bool{return measured_eta_angle(T)<40;}
				<< make_shared<Hist2D>("He3nCentralGammas2","GGangle4",gamma_gamma_sina,gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","GGsin4",gamma_gamma_sina)
                                << make_shared<Hist1D>("He3nCentralGammas2","Events4",Q_axis_full(res))                                                    
                                << make_shared<Hist1D>("He3nCentralGammas2","t4",ggt)                                                                      
                                << make_shared<Hist1D>("He3nCentralGammas2","dt4",ggdt)                                                                    
                                << make_shared<Hist1D>("He3nCentralGammas2","He3MM4",he3mmc)                                                               
                                << make_shared<Hist1D>("He3nCentralGammas2","GMM4",ggmm)                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","GIM4",ggimc)                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","TIM4",he3ggimdiff)                                                            
                                << make_shared<Hist1D>("He3nCentralGammas2","ET4",measured_eta_angle)

				<<[ggmm](WTrack&T)->bool{return (ggmm(T)>2.65)&&(ggmm(T)<3.00);}
				<< make_shared<Hist2D>("He3nCentralGammas2","GGangle5",gamma_gamma_sina,gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","GGsin5",gamma_gamma_sina)
				<< make_shared<Hist1D>("He3nCentralGammas2","Events5",Q_axis_full(res))
				<< make_shared<Hist1D>("He3nCentralGammas2","t5",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","dt5",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","He3MM5",he3mmc)
				<< make_shared<Hist1D>("He3nCentralGammas2","GMM5",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","GIM5",ggimc)
				<< make_shared<Hist1D>("He3nCentralGammas2","TIM5",he3ggimdiff)
				<< make_shared<Hist1D>("He3nCentralGammas2","ET5",measured_eta_angle)

				<<[ggimc](WTrack&T)->bool{return (ggimc(T)>0.45)&&(ggimc(T)<0.65);}
				<< make_shared<Hist2D>("He3nCentralGammas2","GGangle6",gamma_gamma_sina,gamma_gamma_cosa)
				<< make_shared<Hist1D>("He3nCentralGammas2","GGsin6",gamma_gamma_sina)
				<< make_shared<Hist1D>("He3nCentralGammas2","Events6",Q_axis_full(res))                                 
				<< make_shared<Hist1D>("He3nCentralGammas2","t6",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","dt6",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","He3MM6",he3mmc)
				<< make_shared<Hist1D>("He3nCentralGammas2","GMM6",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","GIM6",ggimc)
				<< make_shared<SetOfHists1D>("He3nCentralGammas2","TIM6",Q_axis_full(res),he3ggimdiff)
				<< make_shared<Hist1D>("He3nCentralGammas2","ET6",measured_eta_angle)
				<< make_shared<Hist2D>("He3nCentralGammas2","He3MME6",he3mmc,he3me)
			)
		    )
		)

		<< ( make_shared<ChainCheck>()
		    <<[he3mm](WTrack&T){return (he3mm(T)>0.49)&&(he3mm(T)<0.56);}
		    << make_shared<Hist1D>("He3nCentralGammas2","old_Events1",Q_axis_full(res))
		    << make_shared<Hist1D>("He3nCentralGammas2","old_He3MM1",he3mm)
		    << make_shared<Hist1D>("He3nCentralGammas2","old_GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))
		    <<(make_shared<ChainOr>()
			<< ( make_shared<ChainCheck>()
				<< []()->bool{return gammas.size()>=2;}
				<< [&res]()->bool{
					SortedChain<eta_decay_gg> pairs;
					for(size_t i=0;i<gammas.size();i++)for(size_t j=i+1;j<gammas.size();j++){
						pairs<<eta_decay_gg{.A=gammas[i],.B=gammas[j]};
					}
					two_gamma=pairs[0];
					return true;
				}
				<< make_shared<Hist1D>("He3nCentralGammas2","old_Events2",Q_axis_full(res))
				<< make_shared<Hist1D>("He3nCentralGammas2","old_t2",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_dt2",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_He3MM2",he3mm)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_GMM2",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_GIM2",ggim)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_TIM2",he3ggimdiff)
				<<[ggmm](WTrack&T)->bool{return (ggmm(T)>2.6)&&(ggmm(T)<3.0);}
				<< make_shared<Hist1D>("He3nCentralGammas2","old_Events3",Q_axis_full(res))
				<< make_shared<Hist1D>("He3nCentralGammas2","old_t3",ggt)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_dt3",ggdt)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_He3MM3",he3mm)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_GMM3",ggmm)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_GIM3",ggim)
				<< make_shared<Hist1D>("He3nCentralGammas2","old_TIM3",he3ggimdiff)
                                <<[he3ggimdiff](WTrack&T)->bool{return (he3ggimdiff(T)>0.020);}
                                << make_shared<Hist1D>("He3nCentralGammas2","old_Events4",Q_axis_full(res))                                                                 
                                << make_shared<Hist1D>("He3nCentralGammas2","old_t4",ggt)                                                                                   
                                << make_shared<Hist1D>("He3nCentralGammas2","old_dt4",ggdt)                                                                                 
                                << make_shared<Hist1D>("He3nCentralGammas2","old_He3MM4",he3mm)                                                                             
                                << make_shared<Hist1D>("He3nCentralGammas2","old_GMM4",ggmm)                                                                                
                                << make_shared<Hist1D>("He3nCentralGammas2","old_GIM4",ggim)                                                                                
                                << make_shared<Hist1D>("He3nCentralGammas2","old_TIM4",he3ggimdiff)
                                <<[ggdt](WTrack&T)->bool{return (ggdt(T)<30);}                                                                                                                                           
                                << make_shared<Hist1D>("He3nCentralGammas2","old_Events5",Q_axis_full(res))                                                                
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","old_t5",Q_axis_full(res),ggt)
                                << make_shared<Hist1D>("He3nCentralGammas2","old_dt5",ggdt)                                                                                
                                << make_shared<Hist1D>("He3nCentralGammas2","old_He3MM5",he3mm)                                                                            
                                << make_shared<Hist1D>("He3nCentralGammas2","old_GMM5",ggmm)                                                                               
                                << make_shared<Hist1D>("He3nCentralGammas2","old_GIM5",ggim)                                                                               
                                << make_shared<Hist1D>("He3nCentralGammas2","old_TIM5",he3ggimdiff)
			)
		    )
		)

	    )
	);
}
