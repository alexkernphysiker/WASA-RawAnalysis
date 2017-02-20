// this file is distributed under 
// GPL license
#include "TH1F.h"
#include "TH2F.h"
#include <TLorentzVector.h>
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
#include <math_h/error.h>
#include <Experiment/experiment_conv.h>
#include "trackprocessing.h"
namespace TrackAnalyse{
	using namespace std;
	using namespace MathTemplates;
	Axis::Axis(ValueTrackParamDependent v, double f, double t, unsigned int b){
		value=v;from=f;to=t;bins=b;
		CheckCorrectness();
	}
	Axis::Axis(ValueTrackDependent v, double f, double t, unsigned int b)
		:Axis([v](WTrack&track,const vector<double>&){return v(track);},f,t,b){}
	Axis::Axis(ValueParamDependent v, double f, double t, unsigned int b)
		:Axis([v](WTrack&,const vector<double>&param){return v(param);},f,t,b){}
	Axis::Axis(ValueIndependent v, double f, double t, unsigned int b)
		:Axis([v](WTrack&,const vector<double>&){return v();},f,t,b){}
	Axis::Axis(const Axis& source):Axis(source.value,source.from,source.to,source.bins){}
	Axis::Axis(ValueTrackParamDependent v, const Axis& source):Axis(v,source.from,source.to,source.bins){}
	Axis::Axis(ValueTrackDependent v, const Axis& source):Axis(v,source.from,source.to,source.bins){}
	Axis::Axis(ValueParamDependent v, const Axis& source):Axis(v,source.from,source.to,source.bins){}
	Axis::Axis(ValueIndependent v, const Axis& source):Axis(v,source.from,source.to,source.bins){}
	Axis::~Axis(){}
	void Axis::CheckCorrectness() const{
		if(to<=from)
			throw Exception<Axis>("wrong binning ranges");
		if(0==bins)
			throw Exception<Axis>("there cannot be zero bins");
	}
	unsigned int Axis::count() const{return bins;}
	double Axis::left() const{return from;}
	double Axis::right() const{return to;}
	double Axis::bin_width() const{
		if(0==bins)throw Exception<Axis>("count==0");
		return (to-from)/double(bins);
	}
	double Axis::bin_center(size_t i) const{
		if(i>=bins)throw Exception<Axis>("bin range check error");
		return from+bin_width()*(double(i)+0.5);
	}
	double Axis::operator()(WTrack& T, const vector< double >& P) const{
		return value(T,P);
	}
	ValueTrackParamDependent Axis::valuegetter() const{
		return value;
	}
	bool Axis::FindBinIndex(unsigned int& output, WTrack& T, const vector< double >& P) const{
		double x=operator()(T,P),delta=bin_width()/2.0;
		if(delta<=0)throw Exception<Axis>("delta<=0");
		for(size_t i=0,n=count();i<n;i++){
			double pos=bin_center(i);
			if((x>=(pos-delta))&&(x<(pos+delta))){
				output=i;
				return true;
			}
		}
		return false;
		
	}
	
	bool ITrackParamAnalyse::Process(WTrack&T, vector<double>&P)const{
		Analyse(T,P);
		return true;
	}
	
	Hist1D::Hist1D(const string&&dir,const string&& name, const Axis& x):X(x){
		hist=new TH1F(name.c_str(),"",X.count(),X.left(),X.right());
		gHistoManager->Add(hist,dir.c_str());
	}
	Hist1D::~Hist1D(){}
	void Hist1D::Analyse(WTrack&T, const vector<double>&P)const{hist->Fill(X(T,P));}
	SetOfHists1D::SetOfHists1D(const string&&dir,const string&& name, const Axis& binning, const Axis& x):Z(binning),X(x){
		All=new TH1F((name+"-AllBins").c_str(),"",X.count(),X.left(),X.right());
		gHistoManager->Add(All,dir.c_str());
		OutOfBorder=new TH1F((name+"-OutOfBins").c_str(),"",X.count(),X.left(),X.right());
		gHistoManager->Add(OutOfBorder,dir.c_str());
		for(unsigned int i=0;i<Z.count();i++){
			TH1F*hist=new TH1F((name+"-Bin-"+to_string(i)).c_str(),"",X.count(),X.left(),X.right());
			gHistoManager->Add(hist,dir.c_str());
			Bins.push_back(hist);
		}
	}
	SetOfHists1D::~SetOfHists1D(){}
	void SetOfHists1D::Analyse(WTrack&T, const vector<double>&P)const{
		double x=X(T,P);
		unsigned int i=0;
		if(Z.FindBinIndex(i,T,P)){
			All->Fill(x);
			Bins[i]->Fill(x);
		}else{
			OutOfBorder->Fill(x);
		}
	}
	Hist2D::Hist2D(const string&&dir,const string&& name, const Axis& x, const Axis& y):X(x),Y(y){
		hist=new TH2F(name.c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(hist,dir.c_str());
	}
	Hist2D::~Hist2D(){}
	void Hist2D::Analyse(WTrack&T, const vector<double>&P)const{hist->Fill(X(T,P),Y(T,P));}
	SetOfHists2D::SetOfHists2D(const string&&dir,const string&& name, const Axis& binning, const Axis& x, const Axis& y):Z(binning),X(x),Y(y){
		All=new TH2F((name+"-AllBins").c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(All,dir.c_str());
		OutOfBorder=new TH2F((name+"-OutOfBins").c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(OutOfBorder,dir.c_str());
		for(unsigned int i=0;i<Z.count();i++){
			TH2F*hist=new TH2F((name+"-Bin-"+to_string(i)).c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
			gHistoManager->Add(hist,dir.c_str());
			Bins.push_back(hist);
		}
	}
	SetOfHists2D::~SetOfHists2D(){}
	void SetOfHists2D::Analyse(WTrack&T, const vector<double>&P)const{
		double x=X(T,P);
		double y=Y(T,P);
		unsigned int i=0;
		if(Z.FindBinIndex(i,T,P)){
			All->Fill(x,y);
			Bins[i]->Fill(x,y);
		}else{
			OutOfBorder->Fill(x,y);
		}
	}
	Condition::Condition(ConditionTrackParamDependent func){
		condition=func;
	}
	Condition::Condition(ConditionIndependent func)
		:Condition([func](WTrack&,const vector<double>&){return func();}){}
	Condition::Condition(ConditionParamDependent func)
		:Condition([func](WTrack&,const vector<double>&P){return func(P);}){}
	Condition::Condition(ConditionTrackDependent func)
		:Condition([func](WTrack&T,const vector<double>&){return func(T);}){}
	Condition::~Condition(){}
	bool Condition::Process(WTrack&T, vector<double>&P)const{return condition(T,P);}
	
	Parameter::Parameter(ValueTrackParamDependent f){func=f;}
	Parameter::Parameter(ValueTrackDependent f):Parameter([f](WTrack&T,const vector<double>&){return f(T);}){}
	Parameter::Parameter(ValueParamDependent f):Parameter([f](WTrack&,const vector<double>&P){return f(P);}){}
	Parameter::Parameter(ValueIndependent f):Parameter([f](WTrack&,const vector<double>&){return f();}){}
	Parameter::~Parameter(){}
	bool Parameter::Process(WTrack&T, vector< double >&P)const{
		double x=func(T,P);
		P.push_back(x);
		return true;
	}
	
	AbstractChain::AbstractChain(){}
	AbstractChain::~AbstractChain(){}
	AbstractChain& AbstractChain::operator<<(shared_ptr<ITrackParamProcess> element){
		m_chain.push_back(element);
		return *this;
	}
	const vector<shared_ptr<ITrackParamProcess>>&AbstractChain::chain() const{return m_chain;}
	shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,shared_ptr<ITrackParamProcess>v){
		if(v)
			ch->operator<<(v);
		else
			throw Exception<AbstractChain>("Cannot add process to container");
		return ch;
	}
	bool Chain::Process(WTrack&T, vector< double >&P) const{
		for(auto element:chain())
			element->Process(T,P);
		return true;
	}
	bool ChainCheck::Process(WTrack&T,vector<double>&P) const{
		for(auto element:chain())
			if(!element->Process(T,P))
				return false;
		return true;
	}
	bool ChainAnd::Process(WTrack&T, vector< double >&P) const{
		bool res=true;
		for(auto element:chain())
			res&=element->Process(T,P);
		return res;
	}
	bool ChainOr::Process(WTrack&T, vector< double >&P) const{
		bool res=false;
		for(auto element:chain())
			res|=element->Process(T,P);
		return res;
	}
	ChainBinner::ChainBinner(const Axis& source):m_axis(source){}
	ChainBinner::ChainBinner(Axis&& source):m_axis(source){}
	ChainBinner::~ChainBinner(){}
	bool ChainBinner::Process(WTrack&T, vector<double>&P) const{
		unsigned int i=0;
		if(m_axis.FindBinIndex(i,T,P))
			if(i<chain().size())
				return chain()[i]->Process(T,P);
		return false;
	}
	TrackProcess& TrackProcess::operator<<(shared_ptr< ITrackParamProcess > element){
		m_proc.push_back(element);
		return *this;
	}
	shared_ptr<TrackProcess>operator<<(shared_ptr<TrackProcess>ch,shared_ptr<ITrackParamProcess>v){
		if(v)
			ch->operator<<(v);
		else
			throw Exception<TrackProcess>("Cannot add process to container");
		return ch;
	}
	void TrackProcess::Process(WTrack& T) const{
		vector<double> P;
		for(auto proc:m_proc)
			proc->Process(T,P);
	}
	EventProcess& EventProcess::operator<<(shared_ptr< ITrackParamProcess > element){
		m_proc.push_back(element);
		return *this;
	}
	shared_ptr<EventProcess>operator<<(shared_ptr<EventProcess>ch,shared_ptr<ITrackParamProcess>v){
		if(v)
			ch->operator<<(v);
		else
			throw Exception<EventProcess>("Cannot add process to container");
		return ch;
	}
	bool EventProcess::Process() const{
		WTrack T;
		vector<double> P;
		for(auto proc:m_proc)
			if(!proc->Process(T,P))
				return false;
		return true;
	}
}
