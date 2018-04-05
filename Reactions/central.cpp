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
		Ptotal=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass())+lorentz_byPM(Zero(),Particle::d().mass());
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

		    << [](WTrack&T)->bool{return T.Edep()>=0.04;}
		    << [](WTrack&T)->bool{
			gammas.push_back({.L=lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.),.t=T.Time()});
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas6","GammaEnergyCut",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static eta_decay_ppp six_gamma{.I={.A=He3,.B=He3},.J={.A=He3,.B=He3},.K={.A=He3,.B=He3}};
	Axis 
	he3mm([](){return (Ptotal-He3.L).M();},0.4,0.6,200),
	ggggggdiff([]()->double{return six_gamma.diff();},0.0,0.2,200),
        he3ggggggimdiff([](){return (He3.L+six_gamma.L()).M()-Ptotal.M();},-0.5,0.5,500),
	ggggggim([](){return six_gamma.IM();},0.0,1.0,1000),
        ggggggmm([](){return (Ptotal-six_gamma.L()).M();},0.0,4.0,4000),
        ggggggt([](){return He3.t-six_gamma.t();},-50,50,500),
        ggggggdt([](){return six_gamma.dt();},0,50,250);

	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas6","Events0",Q_axis_full(res))
	    << make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM0",Q_axis_full(res),he3mm)
	    <<[he3mm](WTrack&T){
		return (he3mm(T)>0.49)&&(he3mm(T)<0.56);
	    }
	    << make_shared<Hist1D>("He3nCentralGammas6","Events1",Q_axis_full(res))
	    << make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM1",Q_axis_full(res),he3mm)
	    << make_shared<Hist1D>("He3nCentralGammas6","GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))

	    <<(make_shared<ChainOr>()
	        << ( make_shared<ChainCheck>()
			<<[]()->bool{return gammas.size()>=6;}
	                <<[]()->bool{
	                        SortedChain<eta_decay_ppp> combinations;
        	                for(size_t i=0;i<(gammas.size()-5);++i)for(size_t j=i+1;j<gammas.size();++j){
	                                for(size_t k=i+1;k<(gammas.size()-3);(++k)+=((k==j)?1:0))
        	                        for(size_t l=k+1;l<gammas.size();(++l)+=((l==j)?1:0)){
                                        	for(size_t o=k+1;o<(gammas.size()-1);(++o)+=(((o==j)||(o==l))?1:0))
	                                        for(size_t p=o+1;p<gammas.size();(++p)+=(((p==j)||(p==l))?1:0)){
                                	            combinations<<eta_decay_ppp{
							.I={.A=gammas[i],.B=gammas[j]},
							.J={.A=gammas[k],.B=gammas[l]},
							.K={.A=gammas[o],.B=gammas[p]}
						    };
        	                                }
                	                }
	                        }
	                        six_gamma=combinations[0];
	                        return true;
	                }
			<< make_shared<Hist1D>("He3nCentralGammas6","Events2",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","t2",Q_axis_full(res),ggggggt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","dt2",Q_axis_full(res),ggggggdt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM2",Q_axis_full(res),he3mm)                                        
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","GMM2",Q_axis_full(res),ggggggmm)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff2",ggggggdiff)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","GIM2",Q_axis_full(res),ggggggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM2",Q_axis_full(res),he3ggggggimdiff)
			<< [ggggggt,ggggggdt](WTrack&T){
				return (ggggggt(T)>-50.)&&(ggggggt(T)<50.)&&(ggggggdt(T)<50.);
			}
			<< make_shared<Hist1D>("He3nCentralGammas6","Events3",Q_axis_full(res))
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","t3",Q_axis_full(res),ggggggt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","dt3",Q_axis_full(res),ggggggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM3",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM3",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff3",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM3",Q_axis_full(res),ggggggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM3",Q_axis_full(res),he3ggggggimdiff)
                        <<[ggggggdiff](WTrack&T)->bool{
                                return ggggggdiff(T)<0.080;
                        }
			<< make_shared<Hist1D>("He3nCentralGammas6","Events4",Q_axis_full(res))
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","t4",Q_axis_full(res),ggggggt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","dt4",Q_axis_full(res),ggggggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM4",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM4",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff4",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM4",Q_axis_full(res),ggggggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM4",Q_axis_full(res),he3ggggggimdiff)
                        <<[ggggggmm](WTrack&T)->bool{                                                                 
                                return (ggggggmm(T)>2.6)&&(ggggggmm(T)<3.0);                                          
                        }                                                                                             
			<< make_shared<Hist1D>("He3nCentralGammas6","Events5",Q_axis_full(res))
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","t5",Q_axis_full(res),ggggggt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","dt5",Q_axis_full(res),ggggggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM5",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM5",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff5",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM5",Q_axis_full(res),ggggggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM5",Q_axis_full(res),he3ggggggimdiff)
                        <<[he3ggggggimdiff](WTrack&T)->bool{                                                          
                                return (he3ggggggimdiff(T)>-0.4)&&(he3ggggggimdiff(T)<0.4);                           
                        }
                        << make_shared<Hist1D>("He3nCentralGammas6","Events6",Q_axis_full(res))
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","t6",Q_axis_full(res),ggggggt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","dt6",Q_axis_full(res),ggggggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM6",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM6",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff6",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM6",Q_axis_full(res),ggggggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM6",Q_axis_full(res),he3ggggggimdiff)

                        <<[ggggggdt](WTrack&T)->bool{                                                          
                                return (ggggggdt(T)<20);         
                        }
                        << make_shared<Hist1D>("He3nCentralGammas6","Events7",Q_axis_full(res))
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","t7",Q_axis_full(res),ggggggt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","dt7",Q_axis_full(res),ggggggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM7",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM7",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff7",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM7",Q_axis_full(res),ggggggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","TIM7",Q_axis_full(res),he3ggggggimdiff)
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
		Ptotal=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass())+lorentz_byPM(Zero(),Particle::d().mass());
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

		    << [](WTrack&T)->bool{return T.Edep()>=0.05;}
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
	he3mm([](){return (Ptotal-He3.L).M();},0.4,0.6,200),
	he3ggimdiff([](){return (He3.L+two_gamma.L()).M()-Ptotal.M();},-0.5,0.5,500),
	ggim([](){return two_gamma.IM();},0.0,1.0,1000),
	ggmm([](){return (Ptotal-two_gamma.L()).M();},0.0,4.0,4000),
        ggt([](){return He3.t-two_gamma.t();},-50,50,500),
	ggdt([](){return two_gamma.dt();},0,50,250);

	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas2","Events0",Q_axis_full(res))
	    << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM0",Q_axis_full(res),he3mm)
	    <<[he3mm](WTrack&T){
		return (he3mm(T)>0.49)&&(he3mm(T)<0.56);
	    }
	    << make_shared<Hist1D>("He3nCentralGammas2","Events1",Q_axis_full(res))
	    << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM1",Q_axis_full(res),he3mm)
	    << make_shared<Hist1D>("He3nCentralGammas2","GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))

	    <<(make_shared<ChainOr>()
		<< ( make_shared<ChainCheck>()
			<< []()->bool{return gammas.size()>=2;}
			<< []()->bool{
				SortedChain<eta_decay_gg> pairs;
				for(size_t i=0;i<gammas.size();i++)for(size_t j=i+1;j<gammas.size();j++)
					pairs<<eta_decay_gg{.A=gammas[i],.B=gammas[j]};
				two_gamma=pairs[0];
				return true;
			}
                        << make_shared<Hist1D>("He3nCentralGammas2","Events2",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","t2",Q_axis_full(res),ggt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt2",Q_axis_full(res),ggdt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM2",Q_axis_full(res),he3mm)
                        << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM2",Q_axis_full(res),ggim)
                        << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM2",Q_axis_full(res),ggmm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","TIM2",Q_axis_full(res),he3ggimdiff)
			
		    	<<[ggt,ggdt](WTrack&T){return (ggt(T)>-50.)&&(ggt(T)<50.)&&(ggdt(T)<50.);}
			<< make_shared<Hist1D>("He3nCentralGammas2","Events3",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","t3",Q_axis_full(res),ggt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt3",Q_axis_full(res),ggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM3",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GIM3",Q_axis_full(res),ggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GMM3",Q_axis_full(res),ggmm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","TIM3",Q_axis_full(res),he3ggimdiff)
		    
                        <<[ggmm](WTrack&T)->bool{return (ggmm(T)>2.6)&&(ggmm(T)<3.0);}
			<< make_shared<Hist1D>("He3nCentralGammas2","Events4",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","t4",Q_axis_full(res),ggt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt4",Q_axis_full(res),ggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM4",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GIM4",Q_axis_full(res),ggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GMM4",Q_axis_full(res),ggmm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","TIM4",Q_axis_full(res),he3ggimdiff)

		        <<(make_shared<ChainOr>()
                            << ( make_shared<ChainCheck>()                                                                                                                                                  
                                <<[he3ggimdiff](WTrack&T)->bool{return (he3ggimdiff(T)>-0.040);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events50",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t50",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt50",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM50",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM50",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM50",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM50",Q_axis_full(res),he3ggimdiff)
                                <<[ggdt](WTrack&T)->bool{return (ggdt(T)<20);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events60",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t60",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt60",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM60",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM60",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM60",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM60",Q_axis_full(res),he3ggimdiff)
                            )
                            << ( make_shared<ChainCheck>()                                                                                                                                                  
                                <<[he3ggimdiff](WTrack&T)->bool{return (he3ggimdiff(T)>-0.020);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events51",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t51",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt51",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM51",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM51",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM51",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM51",Q_axis_full(res),he3ggimdiff)
                                <<[ggdt](WTrack&T)->bool{return (ggdt(T)<20);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events61",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t61",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt61",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM61",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM61",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM61",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM61",Q_axis_full(res),he3ggimdiff)
                            )
                            << ( make_shared<ChainCheck>()                                                                                                                                                  
                                <<[he3ggimdiff](WTrack&T)->bool{return (he3ggimdiff(T)> 0.000);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events52",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t52",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt52",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM52",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM52",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM52",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM52",Q_axis_full(res),he3ggimdiff)
                                <<[ggdt](WTrack&T)->bool{return (ggdt(T)<20);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events62",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t62",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt62",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM62",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM62",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM62",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM62",Q_axis_full(res),he3ggimdiff)
                            )
                            << ( make_shared<ChainCheck>()                                                                                                                                                  
                                <<[he3ggimdiff](WTrack&T)->bool{return (he3ggimdiff(T)>+0.020);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events53",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t53",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt53",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM53",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM53",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM53",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM53",Q_axis_full(res),he3ggimdiff)
                                <<[ggdt](WTrack&T)->bool{return (ggdt(T)<20);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events63",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t63",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt63",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM63",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM63",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM63",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM63",Q_axis_full(res),he3ggimdiff)
                            )
                            << ( make_shared<ChainCheck>()                                                                                                                                                  
                                <<[he3ggimdiff](WTrack&T)->bool{return (he3ggimdiff(T)>+0.040);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events54",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t54",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt54",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM54",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM54",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM54",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM54",Q_axis_full(res),he3ggimdiff)
                                <<[ggdt](WTrack&T)->bool{return (ggdt(T)<20);}
                                << make_shared<Hist1D>("He3nCentralGammas2","Events64",Q_axis_full(res))                                                                                                    
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","t64",Q_axis_full(res),ggt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","dt64",Q_axis_full(res),ggdt)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM64",Q_axis_full(res),he3mm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM64",Q_axis_full(res),ggim)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM64",Q_axis_full(res),ggmm)
                                << make_shared<SetOfHists1D>("He3nCentralGammas2","TIM64",Q_axis_full(res),he3ggimdiff)
                            )
                       )
		)
	    )
	);
}
