// this file is distributed under 
// MIT license
#include <fstream>
#include <exception>
#include "data.h"
#include "reconstruction.h"
using namespace std;
const double pbeam_measurement_offset=0.0073;
RealData::RealData()
:BeamMomenta("Time.2.PBeam",[this](){return 1000.0*fHeader->GetTimeInCycle();},[](){return INFINITY;})
{
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
}
RealData::~RealData(){}
bool RealData::DataTypeSpecificEventAnalysis()const{
	CachePBeam((BeamMomenta()/1000.0)+pbeam_measurement_offset);
	return true;
}
bool RealData::DataSpecificTriggerCheck(int n)const{
	return fHeader->TriggerNumSet(n);
}