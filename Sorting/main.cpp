// this file is distributed under 
// GPL license
#include <Wasa.hh>
#include <SorterConfig.hh>
#include <math_h/error.h>
#include "analysis-setup.h"
#include "data.h"
#include "montecarlo.h"
#include "reactions_rundel.h"
#include "reactions_khreptak.h"
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
	Analysis* res=nullptr;
	if("Data"==type)
	    res=new RealData();
	else
	    res=new MonteCarlo(true);
	if(
	    ("RE_He3eta"==type)||
	    ("RE_He3pi0"==type)||
	    ("RE_He3pi0pi0"==type)||
	    ("RE_He3pi0pi0pi0"==type)
	){
	    Rundel::He3_X_reconstruction(*res);;
	}
	if(
	    ("Data"==type)||
	    ("MC_He3eta"==type)||
	    ("MC_He3pi0"==type)||
	    ("MC_He3pi0pi0"==type)||
	    ("MC_He3pi0pi0pi0"==type)
	){
	    Rundel::He3_X_analyse(*res);
	    Rundel::SearchGammaTracks(*res);
	}
	if(
	    ("MC_ppn"==type)||
	    ("Data_ppn"==type)
	){
	    Khreptak::SetupForPPNsp(*res);
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
