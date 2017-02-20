// this file is distributed under 
// GPL license
#ifndef GXMEYPWAEFAIYNHM
#define GXMEYPWAEFAIYNHM
#include "analysis.h"
#include "reconstruction.h"
class RealData:public Analysis{
public:
	RealData();
	virtual ~RealData();
protected:
	virtual bool DataTypeSpecificEventAnalysis()const override;
	virtual bool DataSpecificTriggerCheck(int n)const override;
private:
	REventHeader *fHeader;
	InterpolationBasedReconstruction<> BeamMomenta;
};
#endif
