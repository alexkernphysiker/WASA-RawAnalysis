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
void SearchGamma(Analysis&res){
	const auto&tn=trigger_gammas_central.number;
	static vector<LorentzVector<>> registered;
	static auto Ptotal=LorentzVector<>::zero();
	res.Trigger(tn).pre()<<(make_shared<ChainOr>()
	    <<make_shared<Hist1D>("OnlyCentralGammas","0-Reference",Q_axis_full(res))
	    << [&res](){
		Ptotal=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass())+lorentz_byPM(Zero(),Particle::d().mass());
		registered.clear();
		return true;
	    }
	);
	res.Trigger(tn).per_track()<<(make_shared<ChainOr>()
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << make_shared<Hist1D>("OnlyCentralGammas","GammaEnergy",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=0.08;}
		    << [](WTrack&T)->bool{
			registered.push_back(lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.));
			return true;
		    }
		    << make_shared<Hist1D>("OnlyCentralGammas","GammaEnergy6",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static double MM6G,IM6G,IMdiff6G3Pi;
	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << []()->bool{return registered.size()>0;}
	    << make_shared<Hist1D>("OnlyCentralGammas","GammaCount1For6",Axis([]()->double{return registered.size();},-0.5,9.5,10))
	    <<(make_shared<ChainOr>()
	        << ( make_shared<ChainCheck>()
			<<[]()->bool{return registered.size()==6;}
	                <<[]()->bool{
        	                static const auto M=Particle::pi0().mass();
	                        SortedChain<double> selector;
        	                for(size_t i=0;i<(registered.size()-5);++i)for(size_t j=i+1;j<registered.size();++j){
                	                const auto pair1=registered[i]+registered[j];
                                	const double diff1=pow(pair1.M()-M,2);
	                                for(size_t k=i+1;k<(registered.size()-3);(++k)+=((k==j)?1:0))
        	                        for(size_t l=k+1;l<registered.size();(++l)+=((l==j)?1:0)){
                	                        const auto pair2=registered[k]+registered[l];
                                	        const double diff2=pow(pair2.M()-M,2);
                                        	for(size_t o=k+1;o<(registered.size()-1);(++o)+=(((o==j)||(o==l))?1:0))
	                                        for(size_t p=o+1;p<registered.size();(++p)+=(((p==j)||(p==l))?1:0)){
        	                                    const auto pair3=registered[o]+registered[p];
                        	                    const double diff3=pow(pair3.M()-M,2);
                                	            selector<<sqrt(diff1+diff2+diff3);
        	                                }
                	                }
	                        }
	                        IMdiff6G3Pi=selector[0];
	                        return true;
	                }
			<<[]()->bool{
				auto G=LorentzVector<>::zero();
				for(const auto&g:registered)G+=g;
				IM6G=G.M();
				MM6G=(Ptotal-G).M();				
				return true;
			}
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GMM2",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("OnlyCentralGammas6","GMMPDiff2",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GIM2",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
	                <<[&res]()->bool{
				const double Q=He3eta.P2Q(res.PBeam());
				return (MM6G>2.6+Q)&&(MM6G<2.9+Q);
			}
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GMM3",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("OnlyCentralGammas6","GMMPDiff3",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GIM3",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
	                <<[]()->bool{return IMdiff6G3Pi<0.020;}
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GMM4",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("OnlyCentralGammas6","GMMPDiff4",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GIM4",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
                        <<[&res]()->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (IM6G>0.40+Q)&&(IM6G<0.70+Q);
                        }
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GMM5",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("OnlyCentralGammas6","GMMPDiff5",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("OnlyCentralGammas6","GIM5",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
	        )
	    )
	);
}
void SearchHe3nGamma(Analysis&res){
	const auto&tn=trigger_he3_forward.number;
	static vector<LorentzVector<>> registered_for_two_1,registered_for_two_2,registered;
	static auto Ptotal=LorentzVector<>::zero();
	static particle_kine he3;
	static double he3MM;
	res.Trigger(tn).pre()<<(make_shared<ChainOr>()
	    <<make_shared<Hist1D>("He3nCentralGammas","0-Reference",Q_axis_full(res))
	    << [&res](){
		Ptotal=lorentz_byPM(Z()*res.PBeam(),Particle::p().mass())+lorentz_byPM(Zero(),Particle::d().mass());
		registered_for_two_1.clear();
		registered_for_two_2.clear();
		registered.clear();
		he3MM=INFINITY;
		return true;
	    }
	);
	res.Trigger(tn).per_track()<<(make_shared<ChainOr>()
		<<(make_shared<ChainCheck>()
			<<[](WTrack&T){return T.Type()==kFDC;}
			<<ForwardHe3Reconstruction("CentralGammasandHe3",res,he3)
			<<[&res](){he3MM=He3eta.MissingMass({{.index=0,.E=he3.E,.theta=he3.theta,.phi=he3.phi}},res.PBeam());return true;}
		)
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << make_shared<Hist1D>("He3nCentralGammas","GammaEnergy",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=0.08;}
		    << [](WTrack&T)->bool{
			registered.push_back(lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.));
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas","GammaEnergy6",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=0.25;}
		    << [](WTrack&T)->bool{
			registered_for_two_1.push_back(lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.));
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas","GammaEnergy21",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))

		    << [](WTrack&T)->bool{return T.Edep()>=0.30;}
		    << [](WTrack&T)->bool{
			registered_for_two_2.push_back(lorentz_byPM(direction(T.Phi(),T.Theta())*T.Edep(),0.));
			return true;
		    }
		    << make_shared<Hist1D>("He3nCentralGammas","GammaEnergy22",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static double MM2G_1,IM2G_1,MM2G_2,IM2G_2,MM6G,IM6G,IMdiff6G3Pi;
	res.Trigger(tn).post()<<(make_shared<ChainCheck>()

	    << [](){return isfinite(he3MM);}
	    << make_shared<SetOfHists1D>("He3nCentralGammas","He3MM0",Q_axis_full(res),Axis([]()->double{return he3MM;},0.4,0.6,200))
	    <<[&res](){
		const double Q=He3eta.P2Q(res.PBeam());
		return (he3MM>0.51+Q)&&(he3MM<0.56+Q);
	    }
	    << make_shared<SetOfHists1D>("He3nCentralGammas","He3MM1",Q_axis_full(res),Axis([]()->double{return he3MM;},0.4,0.6,200))
	    << make_shared<Hist1D>("He3nCentralGammas","GammaCount1For21",Axis([]()->double{return registered_for_two_1.size();},-0.5,9.5,10))
	    << make_shared<Hist1D>("He3nCentralGammas","GammaCount1For22",Axis([]()->double{return registered_for_two_2.size();},-0.5,9.5,10))
	    << make_shared<Hist1D>("He3nCentralGammas","GammaCount1For6",Axis([]()->double{return registered.size();},-0.5,9.5,10))

	    <<(make_shared<ChainOr>()
		<< ( make_shared<ChainCheck>()
			<< []()->bool{return registered_for_two_1.size()==2;}
			<< [&res]()->bool{
				IM2G_1=(registered_for_two_1[0]+registered_for_two_1[1]).M();
				MM2G_1=(Ptotal-registered_for_two_1[0]-registered_for_two_1[1]).M();
				return true;
			}
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","He3MM2",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","GIM2",Q_axis_full(res),Axis([](){return IM2G_1;},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","GMM2",Q_axis_full(res),Axis([](){return MM2G_1;},0.0,4.0,4000))
                        <<[&res]()->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (MM2G_1>(2.6+Q))&&(MM2G_1<(2.9+Q));
                        }
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","He3MM3",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","GIM3",Q_axis_full(res),Axis([](){return IM2G_1;},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","GMM3",Q_axis_full(res),Axis([](){return MM2G_1;},0.0,4.0,4000))
                        <<[&res]()->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (IM2G_1>0.45+Q)&&(IM2G_1<0.65+Q);
                        }
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","He3MM4",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","GIM4",Q_axis_full(res),Axis([](){return IM2G_1;},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_1","GMM4",Q_axis_full(res),Axis([](){return MM2G_1;},0.0,4.0,4000))
		)
		<< ( make_shared<ChainCheck>()
			<< []()->bool{return registered_for_two_2.size()==2;}
			<< [&res]()->bool{
				IM2G_2=(registered_for_two_2[0]+registered_for_two_2[1]).M();
				MM2G_2=(Ptotal-registered_for_two_2[0]-registered_for_two_2[1]).M();
				return true;
			}
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","He3MM2",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","GIM2",Q_axis_full(res),Axis([](){return IM2G_2;},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","GMM2",Q_axis_full(res),Axis([](){return MM2G_2;},0.0,4.0,4000))
                        <<[&res]()->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (MM2G_2>2.6+Q)&&(MM2G_2<2.9+Q);
                        }
			<< []()->bool{return (MM2G_2>2.6)&&(MM2G_2<2.9);}
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","He3MM3",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","GIM3",Q_axis_full(res),Axis([](){return IM2G_2;},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","GMM3",Q_axis_full(res),Axis([](){return MM2G_2;},0.0,4.0,4000))
		        <<[&res]()->bool{
                		const double Q=He3eta.P2Q(res.PBeam());
                		return (IM2G_2>0.45+Q)&&(IM2G_2<0.65+Q);
            		}
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","He3MM4",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","GIM4",Q_axis_full(res),Axis([](){return IM2G_2;},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("He3nCentralGammas2_2","GMM4",Q_axis_full(res),Axis([](){return MM2G_2;},0.0,4.0,4000))
		)
	        << ( make_shared<ChainCheck>()
			<<[]()->bool{return registered.size()==6;}
	                <<[]()->bool{
        	                static const auto M=Particle::pi0().mass();
	                        SortedChain<double> selector;
        	                for(size_t i=0;i<(registered.size()-5);++i)for(size_t j=i+1;j<registered.size();++j){
                	                const auto pair1=registered[i]+registered[j];
                                	const double diff1=pow(pair1.M()-M,2);
	                                for(size_t k=i+1;k<(registered.size()-3);(++k)+=((k==j)?1:0))
        	                        for(size_t l=k+1;l<registered.size();(++l)+=((l==j)?1:0)){
                	                        const auto pair2=registered[k]+registered[l];
                                	        const double diff2=pow(pair2.M()-M,2);
                                        	for(size_t o=k+1;o<(registered.size()-1);(++o)+=(((o==j)||(o==l))?1:0))
	                                        for(size_t p=o+1;p<registered.size();(++p)+=(((p==j)||(p==l))?1:0)){
        	                                    const auto pair3=registered[o]+registered[p];
                        	                    const double diff3=pow(pair3.M()-M,2);
                                	            selector<<sqrt(diff1+diff2+diff3);
        	                                }
                	                }
	                        }
	                        IMdiff6G3Pi=selector[0];
	                        return true;
	                }
			<<[]()->bool{
				auto G=LorentzVector<>::zero();
				for(const auto&g:registered)G+=g;
				IM6G=G.M();
				MM6G=(Ptotal-G).M();				
				return true;
			}
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM2",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM2",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff2",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM2",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
                        <<[&res]()->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (MM6G>2.6+Q)&&(MM6G<2.9+Q);
                        }
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM3",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM3",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff3",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM3",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
	                <<[]()->bool{return IMdiff6G3Pi<0.025;}
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM4",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM4",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff4",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM4",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
                        <<[&res]()->bool{
                                const double Q=He3eta.P2Q(res.PBeam());
                                return (IM6G>0.40+Q)&&(IM6G<0.70+Q);
                        }
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","He3MM5",Q_axis_full(res),Axis([](){return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GMM5",Q_axis_full(res),Axis([]()->double{return MM6G;},0.0,4.0,4000))
			<< make_shared<Hist1D>("He3nCentralGammas6","GMMPDiff5",Axis([]()->double{return IMdiff6G3Pi;},0.0,0.2,200))
			<< make_shared<SetOfHists1D>("He3nCentralGammas6","GIM5",Q_axis_full(res),Axis([]()->double{return IM6G;},0.0,1.0,1000))
	        )
	    )
	);
}
