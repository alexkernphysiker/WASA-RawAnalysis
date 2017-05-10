// this file is distributed under 
// GPL license
#include <Wasa.hh>
#include <SorterConfig.hh>
#include <math_h/error.h>
#include <analysis-setup.h>
#include <data.h>
#include <montecarlo.h>
#include <Reactions/analyses.h>
#include <Reconstruction/forward.h>
using namespace std;
using namespace MathTemplates;
int main(int argc, char** argv) {
    if(argc<2)
	return -1;
    int new_c=argc-1;
    char *args[new_c+1];
    args[0]=argv[0];
    string type(argv[1]);
    SetAnalysisType([type](){
	static particle_kinematics p,d,he3;
	Analysis* res=nullptr;
	if("Data"==type)
	    res=new RealData();
	else{
	    if(
		("RE_ppn_qf"==type)
	    )
		res=new MonteCarlo(false);
	    else 
		res=new MonteCarlo(true);
	}
	if(
	    ("RE_He3eta"==type)||
	    ("RE_He3pi0"==type)||
	    ("RE_He3pi0pi0"==type)||
	    ("RE_He3pi0pi0pi0"==type)
	){
	    res->Trigger(0).per_track()<<ForwardHe3Reconstruction(*res,he3);
	}
	if("RE_pd"==type){
	    res->Trigger(0).per_track()<<ForwardPReconstruction(*res,p);
	    res->Trigger(0).per_track()<<ForwardDReconstruction(*res,d);
	}
	if("RE_ppn_qf"==type){
	    res->Trigger(0).per_track()<<ForwardPReconstruction(*res,d);
	}
	if(
	    ("Data"==type)||
	    ("MC_He3eta"==type)||
	    ("MC_He3pi0"==type)||
	    ("MC_He3pi0pi0"==type)||
	    ("MC_He3pi0pi0pi0"==type)||
	    ("MC_He3eta6g"==type)||
	    ("MC_He3pi06g"==type)
	){
	    He3_X_analyse(*res);
	    SearchGamma(*res);
	}
	return res;
    });
    for(int i=1;i<=new_c;i++)
	args[i]=argv[i+1];
    gSorterConfig->ReadCmdLine(new_c,args);
    Wasa::Initialize("AnalysisModule","","RootSorter.log");
    gWasa->AddAnalysis("AnalysisModule","Raw");
    gWasa->Run();
    delete gWasa;
}
