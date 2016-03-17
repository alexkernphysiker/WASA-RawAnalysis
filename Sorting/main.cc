// this file is distributed under 
// MIT license
#include <Wasa.hh>
#include <SorterConfig.hh>
#include <math_h/error.h>
#include "analysis-setup.h"
#include "he3_X.h"
using namespace std;
using namespace MathTemplates;
using namespace ReactionSetup;
int main(int argc, char** argv) {
	if(argc<2)
		return -1;
	int new_c=argc-1;
	char *args[new_c+1];
	args[0]=argv[0];
	string type(argv[1]);
	SetAnalysisType([type](){
		Analysis* res=nullptr;
		if("RE_He3eta"==type)
			res= He3_X_reconstruction(forEta);
		if(
			("RE_He3pi0"==type)||
			("RE_He3pi0pi0"==type)||
			("RE_He3pi0pi0pi0"==type)
		)
			res=He3_X_reconstruction(forPi0);;
		if("MC_He3eta"==type)
			res=He3_X_analyse(forEta);
		if(
			("MC_He3pi0"==type)||
			("MC_He3pi0pi0"==type)||
			("MC_He3pi0pi0pi0"==type)
		)
			res=He3_X_analyse(forPi0);
		if("Data_He3"==type)
			res=He3_X_analyse(forData);
		if(nullptr==res)
			throw Exception<Analysis*>("Cannot create analysis module");
		else return res;
	});
	for(int i=1;i<=new_c;i++)
		args[i]=argv[i+1];
	gSorterConfig->ReadCmdLine(new_c,args);
	Wasa::Initialize("AnalysisModule","","RootSorter.log");
	gWasa->AddAnalysis("AnalysisModule","Raw");
	gWasa->Run();
	delete gWasa;
}