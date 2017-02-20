// this file is distributed under 
// GPL license
#include <fstream>
#include <exception>
#include "data.h"
#include "reconstruction.h"
using namespace std;
const double pbeam_measurement_offset=0.005;
RealData::RealData()
:BeamMomenta("Time.2.PBeam",[this](){return 1000.0*fHeader->GetTimeInCycle();},[](){return INFINITY;}){
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
}
RealData::~RealData(){}
bool RealData::DataTypeSpecificEventAnalysis()const{
	CachePBeam((BeamMomenta()/1000.0)+pbeam_measurement_offset);
	return true;
}
bool RealData::DataSpecificTriggerCheck(int n)const{
	if(0==n)return true;
	return fHeader->TriggerNumSet(n);
}
