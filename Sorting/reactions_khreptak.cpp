// this file is distributed under 
// GPL license
#include <Wasa.hh>
#include <TCutG.h>
#include <vector>
#include <memory>
#include "math_h/functions.h"
#include "math_h/vectors.h"
#include "Experiment/experiment_conv.h"
#include "Kinematics/reactions.h"
#include "trackprocessing.h"
#include "data.h"

#include "reactions_khreptak.h"
namespace Khreptak{
    using namespace std;
    using namespace MathTemplates;
    using namespace TrackAnalyse;
    const Reaction& He3eta(){
	static Reaction res(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return res;
    }
    Axis Q_axis(const Analysis&res){return Axis([&res]()->double{return 1000.0*He3eta().P2Q(res.PBeam());},-70.0,30.0,40);}
    void SetupForPPNsp(Analysis&res){
	auto count_ch_CD=make_shared<size_t>(0);
        auto count_ch_FD=make_shared<size_t>(0);
        auto Theta=make_shared<vector<double>>();
        auto Phi=make_shared<vector<double>>();

	//also pay attention at trigges number
	// for Monte Carlo it will work with any,
	// but for data you must choose the one
        res.Trigger(0).pre()
	    <<[count_ch_CD,count_ch_FD,Theta,Phi]()->bool{
		(*count_ch_CD)=0;
		(*count_ch_FD)=0;
		Theta->clear();
		Phi->clear();
		return true;
	    };

        res.Trigger(0).per_track()
	    <<[count_ch_CD, count_ch_FD, Theta, Phi](WTrack&T)->bool{
		if (T.Type()==kCDC)
		    (*count_ch_CD)++;
		if (T.Type()==kFDC)
		    (*count_ch_FD)++;
		if ((T.Type()==kCDC) || (T.Type()==kFDC)){
		    Theta->push_back(T.Theta());//!!!! methods names were incorrect
		    Phi->push_back(T.Phi());//!!!! and also were used without brackets
		}
		return true;
	    };

        auto count_compl=make_shared<size_t>(0);
	res.Trigger(0).post()
	    <<[Theta,Phi,count_compl]()->bool{
		(*count_compl)=0;
		vector<Vector3<double>> P;
		for(size_t i=0; i<Theta->size(); i++){
		    P.push_back(Vector3<double>::Polar(1, Theta->operator[](i), Phi->operator[](i)));
		}
		for(size_t i=0; i<P.size(); i++){
		    for(size_t j=i+1; j<P.size(); j++){
			//name compl is used somewhere else
			if(pow(Vector3<double>::basis_z().VecP(P[i])*P[j],2)<0.01)
			    (*count_compl)++;
		    }
		}
		return true;//!!!! all these function must return bool
	    }

	    <<make_shared<Hist1D>("ppn","ChargedCentral",Axis([count_ch_CD]()->double{return *count_ch_CD;},-0.5,9.5,10))
	    <<make_shared<Hist1D>("ppn","ChargedForward",Axis([count_ch_FD]()->double{return *count_ch_FD;},-0.5,9.5,10))
	    <<make_shared<Hist1D>("ppn","Complanarity"  ,Axis([count_compl]()->double{return *count_compl;},-0.5,9.5,10));
    }
}
