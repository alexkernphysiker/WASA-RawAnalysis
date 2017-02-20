// this file is distributed under 
// GPL license
#ifndef QTBEVPHGEADUWGMY
#define QTBEVPHGEADUWGMY
#include <memory>
#include <utility>
#include <vector>
#include <CAnalysisModule.hh>
#include <CDataManager.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include <Wasa.hh>
#include <CAnalysisModule.hh>
#include <REventWmcHeader.hh>
#include <REventHeader.hh>
#include <WTrackBank.hh>
#include <WVertexBank.hh>
#include <FDFTHTracks.hh>
#include <CDTracksSimple.hh>
#include "Kinematics/particles.h"
#include "trackprocessing.h"
#include "analysismodule.hh"
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};
class Analysis{
protected:
	Analysis();
public:
	virtual ~Analysis();
	void ProcessEvent();

	class trig_rec{
	public:
		trig_rec();
		~trig_rec();
		TrackAnalyse::EventProcess&pre();
		TrackAnalyse::TrackProcess&per_track();
		TrackAnalyse::EventProcess&post();
	private:
		TrackAnalyse::TrackProcess m_per_track;
		TrackAnalyse::EventProcess m_pre,m_post;
	};
	trig_rec&Trigger(size_t n);
	
	struct Kinematic{Particle particle;double E,Th,Phi;};
	const std::vector<Kinematic>&Vertex(const size_t index)const;
	double PBeam()const;
protected:
	virtual bool DataTypeSpecificEventAnalysis()const=0;
	virtual bool DataSpecificTriggerCheck(int n)const=0;
	void CachePBeam(const double value)const;
	void ClearVerticesCache()const;
	void CacheVertex(const std::vector<Kinematic>&item)const;
private:
	FDFTHTracks* TrackFinderFD;
	CDTracksSimple* CDTrackFinder;
	WTrackBank *fTrackBankFD,*fTrackBankCD;
	CCardWDET *fDetectorTable;
	
	std::vector<trig_rec> m_triggers;

	struct Cache;
	std::shared_ptr<Cache> m_cache;
};
#endif
