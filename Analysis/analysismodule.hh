// this file is distributed under 
// GPL license
#ifndef MZQRVXYK
# define MZQRVXYK
#include <Wasa.hh>
#include <CAnalysisModule.hh>
class Analysis;
class AnalysisModule:public CAnalysisModule{
public:
	AnalysisModule();
	explicit AnalysisModule(const char* name);
	virtual ~AnalysisModule();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	Analysis*m_data;
protected:
	ClassDef(AnalysisModule,0);
};

#endif
