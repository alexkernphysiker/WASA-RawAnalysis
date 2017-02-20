// this file is distributed under 
// GPL license
#ifndef WAFXQSAKXMIIWRNB
#define WAFXQSAKXMIIWRNB
#include "analysis.h"
class MonteCarlo:public Analysis{
public:
	MonteCarlo(const bool kostyl=false);
	virtual ~MonteCarlo();
protected:
	virtual bool DataTypeSpecificEventAnalysis()const override;
	virtual bool DataSpecificTriggerCheck(int n)const override;
private:
	bool f_kostyl;
	REventHeader *fHeader;
	REventWmcHeader *fEventHeader;
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
};
#endif
