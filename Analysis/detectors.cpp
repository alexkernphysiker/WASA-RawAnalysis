// this file is distributed under 
// GPL license
#include <vector>
#include <utility>
#include "math_h/error.h"
#include "detectors.h"
namespace TrackAnalyse {
	using namespace std;
	Forward::plane_data::plane_data(ForwardDetectorPlane p,const string&& n, double thr,double upper)
		:m_plane(p),m_name(n),m_thr(thr),m_axis([p](WTrack&T){return T.Edep(p);},0,upper,500){}
	Forward::plane_data::~plane_data(){}
	ForwardDetectorPlane Forward::plane_data::plane() const{return m_plane;}
	const Axis& Forward::plane_data::axis() const{return m_axis;}
	const string& Forward::plane_data::name() const{return m_name;}
	double Forward::plane_data::threshold() const{return m_thr;}
	double Forward::plane_data::Edep(WTrack&T) const{
		vector<double> P;
		return m_axis(T,P);
	}
	
	Forward::Forward(){
		PlaneData.push_back(plane_data(kFWC1,"FWC1",0.002 ,0.05));
		PlaneData.push_back(plane_data(kFWC2,"FWC2",0.002 ,0.05));
		//PlaneData.push_back(plane_data(kFPC,"FPC",0.002 ,0.05));
		PlaneData.push_back(plane_data(kFTH1,"FTH1",0.0015,0.05));
		PlaneData.push_back(plane_data(kFRH1,"FRH1",0.001, 0.5));
		PlaneData.push_back(plane_data(kFRH2,"FRH2",0.001, 0.5));
	}
	Forward::~Forward(){}
	const Forward& Forward::Get(){
		static Forward result;
		return result;
	}
	Forward::const_iterator Forward::begin() const{
		return PlaneData.begin();
	}
	Forward::const_iterator Forward::cbegin() const{
		return PlaneData.cbegin();
	}
	Forward::const_iterator Forward::end() const{
		return PlaneData.end();
	}
	Forward::const_iterator Forward::cend() const{
		return PlaneData.cend();
	}
	size_t Forward::count() const{return PlaneData.size();}
	const Forward::plane_data& Forward::operator[](ForwardDetectorPlane index) const{
		for(const plane_data&plane:PlaneData)
			if(plane.plane()==index)
				return plane;
		throw MathTemplates::Exception<Forward>("Forward layer not found");
	}
	ForwardDetectorPlane Forward::StoppingLayer(WTrack&T)const{
		int res=-1;
		for(size_t i=0;i<PlaneData.size();i++)
			if(PlaneData[i].Edep(T)>PlaneData[i].threshold()){
				if(res==(int(i)-1))res++;
				else res=-1;
			}
		if(res>=0)
			return PlaneData[res].plane();
		return kForwardError;
	}
	shared_ptr<Chain> Forward::CreateMarker(const string&& dir,const string&& name)const{
		auto res=make_shared<Chain>();
		for(size_t i=1;i<PlaneData.size();i++)
			res << make_shared<Hist2D>(
				const_cast<string&&>(dir),name+"-"+PlaneData[i-1].name()+"-vs-"+PlaneData[i].name(),
				PlaneData[i].axis(),PlaneData[i-1].axis()
			);
		return res;
	}
}
