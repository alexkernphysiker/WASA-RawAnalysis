// this file is distributed under 
// GPL license
#include <functional>
#include "analysismodule.hh"
#include "analysis-setup.h"
using namespace std;
function<Analysis*()> static_analysis_create;
void SetAnalysisType(std::function<Analysis*()>create_analysis){
	static_analysis_create=create_analysis;
}
ClassImp(AnalysisModule);
AnalysisModule::AnalysisModule(){}
AnalysisModule::AnalysisModule(const char* name):CAnalysisModule(name){
	m_data=static_analysis_create();
}
AnalysisModule::~AnalysisModule(){
	if(m_data)
		delete m_data;
}
void AnalysisModule::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	m_data->ProcessEvent();
}
void AnalysisModule::Clear(Option_t *option){fProcessed=kFALSE;}
void AnalysisModule::Print(Option_t *option){}
void AnalysisModule::UserCommand(CCommand * command){}
