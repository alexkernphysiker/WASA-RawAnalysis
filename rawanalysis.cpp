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
    const string type(argv[1]);
    SetAnalysisType([type](){
	Analysis* res=nullptr;
	if(type.substr(0,4)=="Data")
		res=new RealData();
	else{
		if(
			(type.substr(2,6)=="ppn_qf")||
			(type.substr(2,6)=="ppn_qf_")||
			(type.substr(2,5)=="bound")
		)
			res=new MonteCarlo(false);
		else 
			res=new MonteCarlo(true);
	}
	if(type.substr(0,5)=="REHe3")
		res->Trigger(0).per_track()<<ForwardHe3Reconstruction("He3Forward",*res);
	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
		("DataF"==type)
	)
		He3_X_analyse(*res);
	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
		("DataC"==type)
	){
		SearchHe3nGamma(*res);
	}
	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
		("DataCC"==type)
	){
		SearchGamma(*res);
	}
	if(
		("MCpd"==type)||
		("MCppn_qf"==type)||
		("MCppn_qf_"==type)||
		("DataE"==type)
	)
		p_or_d_analyse(*res);
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
