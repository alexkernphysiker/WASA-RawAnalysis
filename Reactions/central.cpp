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
	track_info&A;
	track_info&B;
	inline eta_decay_gg&operator=(const eta_decay_gg&src){A=src.A;B=src.B;return *this;}
	inline const double t()const{return (A.t+B.t)/2.0;}
	inline const LorentzVector<> L()const{return A.L+B.L;}
	inline const double IM()const{return L().M();}
	inline const double diff()const{return abs(IM()-Particle::eta().mass());}
	inline const bool operator<(const eta_decay_gg&other)const{return diff()<other.diff();}
	inline const bool operator>(const eta_decay_gg&other)const{return diff()>other.diff();}
};
struct pi0_decay{
        track_info&A;
        track_info&B;
	inline pi0_decay&operator=(const pi0_decay&src){A=src.A;B=src.B;return *this;}
	inline const double t()const{return (A.t+B.t)/2.0;}
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
	inline const double t()const{return (I.t()+J.t()+K.t())/3.0;}
        inline const LorentzVector<> L()const{return I.L()+J.L()+K.L();}
        inline const double IM()const{return L().M();}
        inline const double diff()const{return sqrt(I.diff()+J.diff()+K.diff());}
        inline const bool operator<(const eta_decay_ppp&other)const{return diff()<other.diff();}
        inline const bool operator>(const eta_decay_ppp&other)const{return diff()>other.diff();}
};
void SearchGamma(Analysis&res){
}
void SearchHe3nGamma(Analysis&res){
	const auto&tn=trigger_he3_forward.number;
	static auto Ptotal=LorentzVector<>::zero();
	static vector<track_info> gammas;
	static particle_kine he3;
	static track_info He3{.L=LorentzVector<>::zero(),.t=INFINITY};
	res.Trigger(tn).pre()<<(make_shared<ChainOr>()
	    <<make_shared<Hist1D>("He3nCentralGammas","0-Reference",Q_axis_full(res))
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
			<<ForwardHe3Reconstruction("CentralGammasandHe3",res,he3)
			<<[](WTrack&T){
				He3={.L=lorentz_byEkM(he3.E,Particle::he3().mass(),direction(he3.phi,he3.theta)),.t=T.Time()};
				return true;
			}
		)
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << make_shared<Hist1D>("He3nCentralGammas","GammaEnergy",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=0.05;}
		    << [](WTrack&T)->bool{
			gammas.push_back({.L=lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.),.t=T.Time()});
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas","GammaEnergyCut",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static eta_decay_gg two_gamma{.A=He3,.B=He3};
	static eta_decay_ppp six_gamma{.I={.A=He3,.B=He3},.J={.A=He3,.B=He3},.K={.A=He3,.B=He3}};
	Axis 
	he3mm([](){return (Ptotal-He3.L).M();},0.4,0.6,200),
	ggim([](){return two_gamma.IM();},0.0,1.0,1000),
	ggmm([](){return (Ptotal-two_gamma.L()).M();},0.0,4.0,4000),
	ggdt([](){return two_gamma.t()-He3.t;},-200,200,400),
	ggggggdiff([]()->double{return six_gamma.diff();},0.0,0.2,200),
        ggggggim([](){return six_gamma.IM();},0.0,1.0,1000),
        ggggggmm([](){return (Ptotal-six_gamma.L()).M();},0.0,4.0,4000),
        ggggggdt([](){return six_gamma.t()-He3.t;},-200,200,400);

	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << [](){return isfinite(He3.t);}
	    << make_shared<Hist1D>("He3nCentralGammas","Events0",Q_axis_full(res))
	    << make_shared<SetOfHists1D>("He3nCentralGammas","He3MM0",Q_axis_full(res),he3mm)
	    <<[&res,he3mm](WTrack&T){
		const double Q=He3eta.P2Q(res.PBeam());
		return (he3mm(T)>0.51+Q)&&(he3mm(T)<0.56+Q);
	    }
	    << make_shared<Hist1D>("He3nCentralGammas","Events1",Q_axis_full(res))
	    << make_shared<SetOfHists1D>("He3nCentralGammas","He3MM1",Q_axis_full(res),he3mm)
	    << make_shared<Hist1D>("He3nCentralGammas","GammaCount",Axis([]()->double{return gammas.size();},-0.5,9.5,10))

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
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt2",Q_axis_full(res),ggdt)
                        << make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM2",Q_axis_full(res),he3mm)
                        << make_shared<SetOfHists1D>("He3nCentralGammas2","GIM2",Q_axis_full(res),ggim)
                        << make_shared<SetOfHists1D>("He3nCentralGammas2","GMM2",Q_axis_full(res),ggmm)
			 
			<< make_shared<Hist1D>("He3nCentralGammas2","Events3",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt3",Q_axis_full(res),ggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM3",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GIM3",Q_axis_full(res),ggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GMM3",Q_axis_full(res),ggmm)
                        <<[ggmm](WTrack&T)->bool{
                                return (ggmm(T)>2.6)&&(ggmm(T)<2.9);
                        }
			<< make_shared<Hist1D>("He3nCentralGammas2","Events4",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt4",Q_axis_full(res),ggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM4",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GIM4",Q_axis_full(res),ggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GMM4",Q_axis_full(res),ggmm)
                        <<[&res,ggim](WTrack&T)->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (ggim(T)>0.45+Q)&&(ggim(T)<0.65+Q);
                        }
			<< make_shared<Hist1D>("He3nCentralGammas2","Events5",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","dt5",Q_axis_full(res),ggdt)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","He3MM5",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GIM5",Q_axis_full(res),ggim)
			<< make_shared<SetOfHists1D>("He3nCentralGammas2","GMM5",Q_axis_full(res),ggmm)
		)
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
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM2",Q_axis_full(res),he3mm)                                        
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","GMM2",Q_axis_full(res),ggggggmm)
                        << make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff2",ggggggdiff)
                        << make_shared<SetOfHists1D>("He3nCentralGammas6","GIM2",Q_axis_full(res),ggggggim)

			<< make_shared<Hist1D>("He3nCentralGammas6","Events3",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM3",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM3",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff3",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM3",Q_axis_full(res),ggggggim)
                        <<[&res,ggggggmm](WTrack&T)->bool{
                                return (ggggggmm(T)>2.6)&&(ggggggmm(T)<2.9);
                        }
			<< make_shared<Hist1D>("He3nCentralGammas6","Events4",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM4",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM4",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff4",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM4",Q_axis_full(res),ggggggim)
	                <<[ggggggdiff](WTrack&T)->bool{
				return ggggggdiff(T)<0.080;
			}
			<< make_shared<Hist1D>("He3nCentralGammas6","Events5",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM5",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM5",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff5",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM5",Q_axis_full(res),ggggggim)
                        <<[&res,ggggggim](WTrack&T)->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (ggggggim(T)>0.40+Q)&&(ggggggim(T)<0.70+Q);
                        }
			<< make_shared<Hist1D>("He3nCentralGammas6","Events6",Q_axis_full(res))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM6",Q_axis_full(res),he3mm)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM6",Q_axis_full(res),ggggggmm)
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff6",ggggggdiff)
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM6",Q_axis_full(res),ggggggim)
	        )
	    )
	);
}
