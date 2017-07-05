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
	if(
		("DataL"==type)||
		("DataR"==type)
	)
		res=new RealData();
	else{
		if(
			(type.substr(2,6)=="ppn_qf")||
			(type.substr(2,5)=="bound")
		)
			res=new MonteCarlo(false);
		else 
			res=new MonteCarlo(true);
	}
	if(
		("REHe3eta"==type)||
		("REHe3pi0"==type)||
		("REHe3pi0pi0"==type)||
		("REHe3pi0pi0pi0"==type)
	)
		res->Trigger(0).per_track()<<ForwardHe3Reconstruction(*res);

	if("REpd"==type){
		res->Trigger(0).per_track()<<ForwardPReconstruction(*res);
		res->Trigger(0).per_track()<<ForwardDReconstruction(*res);
	}
	if("REppn_qf"==type)
		res->Trigger(0).per_track()<<ForwardPReconstruction(*res);

	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
		("DataL"==type)
	)
		He3_X_analyse(*res);

	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
		("DataR"==type)
	)
		SearchGamma(*res);

	if(
		("MCpd"==type)||
		("MCppn_qf"==type)||
		("DataL"==type)
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
