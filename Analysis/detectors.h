// this file is distributed under 
// GPL license
#ifndef HGDOYENRKILTMEYD
#define HGDOYENRKILTMEYD
#include <string>
#include <functional>
#include <TH2F.h>
#include "trackprocessing.h"
enum ForwardDetectorPlane{
	kForwardError=0,
	kFWC1 = 10, kFWC2 = 11, kFTH1 = 1, kFTH2 =2, kFTH3 = 3, 
	kFRH1 = 4, kFRH2 = 5, kFRH3 = 6, kFRH4 = 7,kFRH5 = 8, kFVH =9
};
namespace TrackAnalyse{
	class Forward{
	protected:
		Forward();
	public:
		static const Forward&Get();
		virtual ~Forward();
		class plane_data{
			friend class Forward;
		public:
			plane_data(ForwardDetectorPlane p,const std::string&&n,double thr,double upper);
			~plane_data();
			double Edep(WTrack&)const;
		protected:
			ForwardDetectorPlane plane()const;
			const std::string&name()const;
			const Axis&axis()const;
			double threshold()const;
		private:
			ForwardDetectorPlane m_plane;
			std::string m_name;
			double m_thr;
			Axis m_axis;
		};
		typedef std::vector<plane_data>::const_iterator const_iterator;
		const_iterator begin()const;
		const_iterator cbegin()const;
		const_iterator end()const;
		const_iterator cend()const;
		size_t count()const;
		const plane_data&operator[](ForwardDetectorPlane)const;
		ForwardDetectorPlane StoppingLayer(WTrack&)const;
		std::shared_ptr<Chain> CreateMarker(const std::string&&dir,const std::string&&name)const;
	private:
		std::vector<plane_data> PlaneData;
	};
}
#endif
