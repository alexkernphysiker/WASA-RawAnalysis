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
	const auto&tn=trigger_he3_forward.number;
	static vector<Vector4<double>> registered;
	static auto Ptotal=Vector4<double>::zero();
	static double TotalE;
	static double AcceptedE;
	static particle_kine he3;
	static double he3MM;
	res.Trigger(tn).pre()<<(make_shared<ChainOr>()
	    <<make_shared<Hist1D>("CentralGammas","0-Reference",Q_axis_full(res))
	    << [&res](){
		Ptotal=Get4Vector({.particle=Particle::p(),.E=Particle::he3().P2E(res.PBeam()),.theta=0,.phi=0})
		    +Get4Vector({.particle=Particle::d(),.E=0,.theta=0,.phi=0});
		registered.clear();
		TotalE=0;
		AcceptedE=0;
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
		    << [](WTrack&T)->bool{return T.Edep()>=0.050;}
		    << [](WTrack&T)->bool{
			registered.push_back(Get4Vector({.particle=Particle::gamma(),.E=T.Edep(),.theta=T.Theta(),.phi=T.Phi()}));
			TotalE+=T.Edep();
			return true;
		    }
		    << make_shared<Hist1D>("CentralGammas","GammaEnergy2",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.6,800))
		)
	);
	static point<double> gamma_pair(INFINITY,INFINITY);
	static point<double> gamma_pair_mm(INFINITY,INFINITY);
	static point<double> pi0_pair(INFINITY,INFINITY);
	static point<double> pi0_pair_mm(INFINITY,INFINITY);
	static point<double> pi0_triple(INFINITY,INFINITY);
	static point<double> pi0_triple_mm(INFINITY,INFINITY);
	res.Trigger(tn).post()<<(make_shared<ChainCheck>()
	    << []()->bool{return registered.size()>1;}
	    <<[&res](){
		if(!isfinite(he3MM))return false;
		const double Q=He3eta.P2Q(res.PBeam());
		return (he3MM>0.50+Q*10.)&&(he3MM<0.55+Q*10.);
	    }
	    << make_shared<SetOfHists1D>("CentralGammas","He3MM",Q_axis_full(res),Axis([]()->double{return he3MM;},0.4,0.6,200))
	    << make_shared<Hist1D>("CentralGammas","GammaCount",Axis([]()->double{return registered.size();},-0.5,9.5,10))
	    << make_shared<Hist1D>("CentralGammas","GammaTotalEnergy",Axis([](){return TotalE;},0.0,1.6,800))
	    <<(make_shared<ChainOr>()
		<< ( make_shared<ChainCheck>()
			<< []()->bool{return registered.size()==2;}
			<< []()->bool{AcceptedE=0;return true;}
			<< []()->bool{
				static const auto M=Particle::eta().mass();
				SortedPoints<double> table,table2,e_table;
				for(size_t i=0;i<(registered.size()-1);i++)for(size_t j=i+1;j<registered.size();j++){
				    const double im=(registered[i]+registered[j]).length4();
				    const double mm=(Ptotal-(registered[i]+registered[j])).length4();
				    const double e=registered[i].time_component()+registered[j].time_component();
				    const double diff=pow(im-M,2);
				    table<<point<double>(sqrt(diff),im);
				    table2<<point<double>(sqrt(diff),mm);
				    e_table<<point<double>(diff,e);
				}
				gamma_pair=table[0];
				gamma_pair_mm=table2[0];
				AcceptedE=e_table[0].Y();
				return true;
			}
			<< make_shared<Hist1D>("CentralGammas2","TotalEnergy",Axis([](){return AcceptedE;},0.0,1.6,800))
			<< make_shared<SetOfHists1D>("CentralGammas2","He3MM",Q_axis_full(res),Axis([]()->double{return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("CentralGammas2","GInvMass",Q_axis_full(res),Axis([]()->double{return gamma_pair.Y();},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("CentralGammas2","GMMass",Q_axis_full(res),Axis([]()->double{return gamma_pair_mm.Y();},0.0,4.0,4000))
	        	<< make_shared<SetOfHists1D>("CentralGammas2","GMMDiff",Q_axis_full(res),Axis([]()->double{return gamma_pair.X();},0.0,0.2,200))
		)
	        << ( make_shared<ChainCheck>()
			<<[]()->bool{return registered.size()==6;}
	                <<[]()->bool{AcceptedE=0;return true;}
	                <<[]()->bool{
        	                static const auto M=Particle::pi0().mass();
	                        SortedPoints<double> selector,selector_mm,e_table;
        	                for(size_t i=0;i<(registered.size()-5);++i)for(size_t j=i+1;j<registered.size();++j){
                	                const auto pair1=registered[i]+registered[j];
                        	        const auto e1=registered[i].time_component()+registered[j].time_component();
                                	const double diff1=pow(pair1.length4()-M,2);
	                                for(size_t k=i+1;k<(registered.size()-3);(++k)+=((k==j)?1:0))
        	                        for(size_t l=k+1;l<registered.size();(++l)+=((l==j)?1:0)){
                	                        const auto pair2=registered[k]+registered[l];
                        	                const auto e2=registered[k].time_component()+registered[l].time_component();
                                	        const double diff2=pow(pair2.length4()-M,2);
                                        	for(size_t o=k+1;o<(registered.size()-1);(++o)+=(((o==j)||(o==l))?1:0))
	                                        for(size_t p=o+1;p<registered.size();(++p)+=(((p==j)||(p==l))?1:0)){
        	                                    const auto pair3=registered[o]+registered[p];
                	                            const auto e3=registered[o].time_component()+registered[p].time_component();
                        	                    const double diff3=pow(pair3.length4()-M,2);
                                	            selector<<point<double>(sqrt(diff1+diff2+diff3),(pair1+pair2+pair3).length4());
                                        	    selector_mm<<point<double>(sqrt(diff1+diff2+diff3),(Ptotal-(pair1+pair2+pair3)).length4());
	                                            e_table<<point<double>(diff1+diff2+diff3,e1+e2+e3);
        	                                }
                	                }
	                        }
	                        pi0_triple=selector[0];
	                        pi0_triple_mm=selector_mm[0];
	                        AcceptedE=e_table[0].Y();
	                        return true;
	                }
	                <<[]()->bool{return pi0_triple.X()<0.030;}
			<< make_shared<Hist1D>("CentralGammas6","TotalEnergy",Axis([](){return AcceptedE;},0.0,1.6,800))
			<< make_shared<SetOfHists1D>("CentralGammas6","He3MM",Q_axis_full(res),Axis([]()->double{return he3MM;},0.4,0.6,200))
			<< make_shared<SetOfHists1D>("CentralGammas6","GInvMass",Q_axis_full(res),Axis([]()->double{return pi0_triple.Y();},0.0,1.0,1000))
			<< make_shared<SetOfHists1D>("CentralGammas6","GMMass",Q_axis_full(res),Axis([]()->double{return pi0_triple_mm.Y();},0.0,4.0,4000))
			<< make_shared<SetOfHists1D>("CentralGammas6","GMMDiff",Q_axis_full(res),Axis([]()->double{return pi0_triple.X();},0.0,0.2,200))
	        )
	    )
	);
}
