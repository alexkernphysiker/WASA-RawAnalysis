// this file is distributed under 
// GPL license
#include <string>
#include <sstream>
#include <Parameters/parameters.h>
#include <Wasa.hh>
#include <SorterConfig.hh>
#include <analysis-setup.h>
#include <data.h>
#include <montecarlo.h>
#include <Reactions/analyses.h>
#include <Reconstruction/forward.h>
using namespace std;
using namespace MathTemplates;
int main(int argc, char** argv) {
    if(argc<3)
	return -1;
    int new_c=argc-2;
    char *args[new_c+1];
    args[0]=argv[0];
    const string type(argv[1]),param_setup(argv[2]);
    if(param_setup!="_"){
	const string param_num=param_setup.substr(0,2);
	const string param_mod=param_setup.substr(2,1);
	istringstream ss(param_num);
	size_t index;
	ss>>index;
	ParameterMode mode=ParameterMode::param_normal;
	if(param_mod=="+")mode=ParameterMode::param_up;
	if(param_mod=="-")mode=ParameterMode::param_down;
	ChangedParameter(index,mode);
    }
    SetAnalysisType([type](){
	Analysis* res=nullptr;
	if(
		(type.substr(0,4)=="Data")||
		(type.substr(0,4)=="PRES")
	){
		 res=new RealData();
	}else{
		if(
			(type.substr(2,3)=="ppn")||
			(type.substr(2,3)=="pd_")||
			(type.substr(2,5)=="bound")||
			(type.substr(2,7)=="He3eta-")
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
                ("DataAll"==type)||
		("DataF"==type)
	)
		He3_X_analyse(*res);
	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
                ("DataAll"==type)||
		("DataC"==type)
	)
		Search3He2Gamma(*res);
	if(
		(type.substr(0,5)=="MCHe3")||
		(type.substr(0,7)=="MCbound")||
                ("DataAll"==type)||
		("DataCC"==type)
	)
		Search3He6Gamma(*res);
	if(
		(type.substr(0,3)=="MCp")||
                ("DataAll"==type)||
		("DataQ"==type)
	)
		qe_central_analysis(*res);
	if(
                ("PRESELECTION"==type)
	)
		Preselection(*res);
	return res;
    });
    for(int i=1;i<=new_c;i++)
	args[i]=argv[i+2];
    gSorterConfig->ReadCmdLine(new_c,args);
    Wasa::Initialize("AnalysisModule","","RootSorter.log");
    gWasa->AddAnalysis("AnalysisModule","Raw");
    gWasa->Run();
    delete gWasa;
}
