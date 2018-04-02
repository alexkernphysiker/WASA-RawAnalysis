// this file is distributed under 
// GPL license
#include <iostream>
#include <Wasa.hh>
#include <Experiment/experiment_conv.h>
#include <data.h>
#include "analyses.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;
void Preselection(Analysis&res){
	Int_t RunNumber = SorterOption::GetIntValue("RunNumber");
	//ACHTUNG: HARDCODE!!!!!
	TString Tstream1 = Form("file:/home/PRESEL/run_%i.bz2", RunNumber);
	gWasa->AddOutput("run_presel",Tstream1);

	res.Trigger(trigger_he3_forward.number).pre()<< [](){gWasa->SaveEvent("run_presel");return true;};
	res.Trigger(trigger_elastic1.number).pre()<< [](){gWasa->SaveEvent("run_presel");return true;};
	res.Trigger(trigger_elastic2.number).pre()<< [](){gWasa->SaveEvent("run_presel");return true;};
}
