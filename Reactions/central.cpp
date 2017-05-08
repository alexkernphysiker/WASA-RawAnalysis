// this file is distributed under 
// GPL license
#include <Wasa.hh>
#include <TCutG.h>
#include <math_h/functions.h>
#include <math_h/tabledata.h>
#include <math_h/error.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include <trackprocessing.h>
#include <detectors.h>
#include <reconstruction.h>
#include <data.h>
#include "forward.h"
#include "central.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;
const Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
Axis Q_axis_full(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta.P2Q(res.PBeam());},-70.0,30.0,40);}

void SearchGamma(Analysis&res){
	const auto&tn=trigger_gammas_central.number;
	static vector<Vector4<double>> registered;
	res.Trigger(tn).pre()
	    << [](){registered.clear(); return true;}
	;
	res.Trigger(tn).per_track()
		<<(make_shared<ChainCheck>()
		    << [](WTrack&T)->bool{return T.Type()==kCDN;}
		    << [](WTrack&T)->bool{return T.Edep()>=0.010;}
		    << [](WTrack&T)->bool{
			registered.push_back(Get4Vector({.particle=Particle::gamma(),.E=T.Edep(),.theta=T.Theta(),.phi=T.Phi()}));
			return true;
		    }
		    << make_shared<Hist1D>("CentralGammas","GammaEnergy",Axis([](WTrack&T)->double{return T.Edep();},0.0,1.0,1000))
		)
	;
	
	static point<double> gamma_pair(INFINITY,INFINITY);
	static point<double> pi0_pair(INFINITY,INFINITY);
	res.Trigger(tn).post() 
		<< make_shared<Hist1D>("CentralGammas","GammaCount",Axis([]()->double{return registered.size();},-0.5,9.5,10))
		<< make_shared<Hist1D>("CentralGammas","GammaTotalEnergy",Axis([]()->double{
			double res=0;
			for(const auto&g:registered)res+=g.space_component().mag();
			return res;
		},0.,1.5,1500))
		<< (make_shared<ChainCheck>()
		    << []()->bool{return registered.size()>=2;}
		    << []()->bool{
			SortedPoints<double> table;
			for(size_t i=0;i<(registered.size()-1);i++)for(size_t j=i+1;j<registered.size();j++){
			    double im=(registered[i]+registered[j]).length4();
			    double rest=0;
			    for(size_t k=0;k<registered.size();k++)if((k!=i)&&(k!=j))
				rest+=registered[k].space_component().mag();
			    table<<point<double>(rest,im);
			}
			gamma_pair=table[0];
			return true;
		    }
		    <<[]()->bool{return gamma_pair.X()<0.200;}
                    << make_shared<SetOfHists1D>("CentralGammas","InvMass2Gamma",
			Q_axis_full(res),
			Axis([]()->double{return gamma_pair.Y();},0.0,1.0,1000)
		    )
                    << make_shared<SetOfHists1D>("CentralGammas","RestAfter2Gamma",
			Q_axis_full(res),
			Axis([]()->double{return gamma_pair.X();},0.0,0.2,200)
		    )
		)
		<< (make_shared<ChainCheck>()
			<<[]()->bool{return registered.size()>=4;}
			<<[]()->bool{
				static auto M=Particle::pi0().mass();
				SortedPoints<double> selector;
				for(size_t i=0;i<(registered.size()-3);++i)
				for(size_t j=i+1;j<registered.size();++j){
					const auto pair1=registered[i]+registered[j];
					const double diff1=pow(pair1.length4()-M,2);
					for(size_t k=i+1;k<(registered.size()-1);(++k)+=((k==j)?1:0))
					for(size_t l=k+1;l<registered.size();(++l)+=((l==j)?1:0)){
						const auto pair2=registered[k]+registered[l];
						const double diff2=pow(pair2.length4()-M,2);
						selector<<point<double>(
							sqrt(diff1+diff2),
							(pair1+pair2).length4()
						);
					}
				}
				pi0_pair=selector[0];
				return true;
			}
			<<[]()->bool{return pi0_pair.X()<0.200;}
			<< make_shared<SetOfHists1D>("CentralGammas","InvMass2Pairs",
				Q_axis_full(res),
				Axis([]()->double{return pi0_pair.Y();},0.0,1.0,1000)
			)
			<< make_shared<SetOfHists1D>("CentralGammas","Diff2Pairs",
				Q_axis_full(res),
				Axis([]()->double{return pi0_pair.X();},0.0,0.2,200)
			)
		)
	;
}
