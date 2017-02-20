// this file is distributed under 
// GPL license
#include <iostream>
#include <functional>
#include "math_h/error.h"
#include "analysis.h"
#define static_right(A) (static_cast<decltype(A)&&>(A))
using namespace std;
using namespace MathTemplates;
Analysis::trig_rec::trig_rec(){}
Analysis::trig_rec::~trig_rec(){}
TrackAnalyse::EventProcess& Analysis::trig_rec::pre(){return m_pre;}
TrackAnalyse::TrackProcess& Analysis::trig_rec::per_track(){return m_per_track;}
TrackAnalyse::EventProcess& Analysis::trig_rec::post(){return m_post;}
struct Analysis::Cache{
	double p_beam_cache;
	vector<vector<Analysis::Kinematic>> vertex_cache;
};

Analysis::Analysis(){
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
	m_cache=make_shared<Cache>();
	for(size_t i=0;i<33;i++)
		m_triggers.push_back(trig_rec());
}
Analysis::~Analysis(){}
Analysis::trig_rec& Analysis::Trigger(size_t n){
	if(n>=m_triggers.size())
		throw Exception<Analysis>("Attempt to access bad trigger record");
	return m_triggers[n];
}

void Analysis::ProcessEvent(){
	m_cache->p_beam_cache=INFINITY;
	if(DataTypeSpecificEventAnalysis()){
		vector<bool> trig;
		for(size_t i=0;i<m_triggers.size();i++)
			trig.push_back(DataSpecificTriggerCheck(i));
		for(size_t i=0;i<m_triggers.size();i++)
			if(trig[i])
				trig[i]=m_triggers[i].pre().Process();
		vector<WTrackBank*> BANK;
		BANK.push_back(fTrackBankCD);
		BANK.push_back(fTrackBankFD);
		for(WTrackBank*bank:BANK){
			WTrackIter iterator(bank);
			while(WTrack* track = dynamic_cast<WTrack*> (iterator.Next()))
				for(size_t i=0;i<m_triggers.size();i++)
					if(trig[i])
						m_triggers[i].per_track().Process(*track);
		}
		for(size_t i=0;i<m_triggers.size();i++)
			if(trig[i])
				m_triggers[i].post().Process();
	}
}
double Analysis::PBeam()const{return m_cache->p_beam_cache;}
void Analysis::CachePBeam(const double value)const{
	if(value>0)m_cache->p_beam_cache=value;
	else throw Exception<Analysis>("Wrong beam momentum value");
}
const vector<Analysis::Kinematic>&Analysis::Vertex(const size_t index) const{
	static vector<Analysis::Kinematic> dummy;
	if(index>=m_cache->vertex_cache.size())
		return dummy;
	return m_cache->vertex_cache[index];
}
void Analysis::ClearVerticesCache() const{
	m_cache->vertex_cache.clear();
}
void Analysis::CacheVertex(const vector<Kinematic>& item) const{
	m_cache->vertex_cache.push_back(vector<Kinematic>());
	for(const auto&i:item)
		m_cache->vertex_cache[m_cache->vertex_cache.size()-1].push_back(i);
}



