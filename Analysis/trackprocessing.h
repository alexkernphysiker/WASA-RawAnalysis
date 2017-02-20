// this file is distributed under 
// GPL license
#ifndef xoyoyptv
#define xoyoyptv
#include <functional>
#include <vector>
#include <memory>
#include <string>
#include <WTrack.hh>
#include <TH1F.h>
#include <TH2F.h>
namespace TrackAnalyse{
	//WTrack cannot be transfered as const because
	//it does not contain const methods 
	//(even those ones that really should be)
	typedef std::function<double()> ValueIndependent;
	typedef std::function<double(WTrack&)> ValueTrackDependent;
	typedef std::function<double(const std::vector<double>&)> ValueParamDependent;
	typedef std::function<double(WTrack&,const std::vector<double>&)> ValueTrackParamDependent;
	typedef std::function<bool()> ConditionIndependent;
	typedef std::function<bool(WTrack&)> ConditionTrackDependent;
	typedef std::function<bool(const std::vector<double>&)> ConditionParamDependent;
	typedef std::function<bool(WTrack&,const std::vector<double>&)> ConditionTrackParamDependent;
	class Axis{
	public:
		Axis(ValueTrackParamDependent v,double f, double t,unsigned int b);
		Axis(ValueTrackDependent v,double f, double t,unsigned int b);
		Axis(ValueParamDependent v,double f, double t,unsigned int b);
		Axis(ValueIndependent v,double f, double t,unsigned int b);
		Axis(const Axis&source);
		Axis(ValueTrackParamDependent v,const Axis&source);
		Axis(ValueTrackDependent v,const Axis&source);
		Axis(ValueParamDependent v,const Axis&source);
		Axis(ValueIndependent v,const Axis&source);
		~Axis();
		double left()const;
		double right()const;
		unsigned int count()const;
		double operator()(WTrack&T,const std::vector<double>&P)const;
		ValueTrackParamDependent valuegetter()const;
		double bin_width()const;
		double bin_center(size_t i)const;
		bool FindBinIndex(unsigned int&output,WTrack&T,const std::vector<double>&P)const;
	private:
		void CheckCorrectness()const;
		ValueTrackParamDependent value;
		double from;
		double to;
		unsigned int bins;
	};
	class ITrackParamProcess{
	public:
		virtual bool Process(WTrack&,std::vector<double>&)const=0;
		virtual ~ITrackParamProcess(){}
	};
	class ITrackParamAnalyse:public ITrackParamProcess{
	public:
		virtual bool Process(WTrack&,std::vector<double>&)const final;
		virtual ~ITrackParamAnalyse(){}
	protected:
		virtual void Analyse(WTrack&,const std::vector<double>&)const =0;
	};
	
	class Hist1D:public ITrackParamAnalyse{
	public:
		Hist1D(const std::string&&dir,const std::string&&name,const Axis&x);
		virtual ~Hist1D();
	protected:
		virtual void Analyse(WTrack&,const std::vector<double>&)const override;
	private:
		TH1F *hist;
		Axis X;
	};
	class SetOfHists1D:public ITrackParamAnalyse{
	public:
		SetOfHists1D(const std::string&&dir,const std::string&&name,const Axis&binning,const Axis&x);
		virtual ~SetOfHists1D();
	protected:
		virtual void Analyse(WTrack&,const std::vector<double>&)const override;
	private:
		Axis Z,X;
		TH1F *All,*OutOfBorder;
		std::vector<TH1F*> Bins;
	};

	class Hist2D:public ITrackParamAnalyse{
	public:
		Hist2D(const std::string&&dir,const std::string&&name,const Axis&x,const Axis&y);
		virtual ~Hist2D();
	protected:
		virtual void Analyse(WTrack&,const std::vector<double>&)const override;
	private:
		TH2F *hist;
		Axis X,Y;
	};
	class SetOfHists2D:public ITrackParamAnalyse{
	public:
		SetOfHists2D(const std::string&&dir,const std::string&&name,const Axis&binning,const Axis&x,const Axis&y);
		virtual ~SetOfHists2D();
	protected:
		virtual void Analyse(WTrack&,const std::vector<double>&)const override;
	private:
		Axis Z,X,Y;
		TH2F *All,*OutOfBorder;
		std::vector<TH2F*> Bins;
	};
	class Condition:public ITrackParamProcess{
	public:
		Condition(ConditionTrackParamDependent func);
		Condition(ConditionTrackDependent func);
		Condition(ConditionParamDependent func);
		Condition(ConditionIndependent func);
		virtual ~Condition();
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	private:
		ConditionTrackParamDependent condition;
	};
	class Parameter:public ITrackParamProcess{
	public:
		Parameter(ValueTrackParamDependent f);
		Parameter(ValueParamDependent f);
		Parameter(ValueTrackDependent f);
		Parameter(ValueIndependent f);
		virtual ~Parameter();
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	private:
		ValueTrackParamDependent func;
	};
	class AbstractChain:public ITrackParamProcess{
	protected:
		AbstractChain();
	public:
		virtual ~AbstractChain();
		AbstractChain&operator<<(std::shared_ptr<ITrackParamProcess>element);
	protected:
		const std::vector<std::shared_ptr<ITrackParamProcess>>&chain()const;
	private:
		std::vector<std::shared_ptr<ITrackParamProcess>> m_chain;
	};
	std::shared_ptr<AbstractChain>operator<<(std::shared_ptr<AbstractChain>ch,std::shared_ptr<ITrackParamProcess>v);
	inline std::shared_ptr<AbstractChain>operator<<(std::shared_ptr<AbstractChain>ch,ConditionTrackParamDependent f){
		return ch<<std::dynamic_pointer_cast<ITrackParamProcess>(std::make_shared<Condition>(f));
	}
	inline std::shared_ptr<AbstractChain>operator<<(std::shared_ptr<AbstractChain>ch,ConditionTrackDependent f){
		return ch<<std::dynamic_pointer_cast<ITrackParamProcess>(std::make_shared<Condition>(f));
	}
	inline std::shared_ptr<AbstractChain>operator<<(std::shared_ptr<AbstractChain>ch,ConditionParamDependent f){
		return ch<<std::dynamic_pointer_cast<ITrackParamProcess>(std::make_shared<Condition>(f));
	}
	inline std::shared_ptr<AbstractChain>operator<<(std::shared_ptr<AbstractChain>ch,ConditionIndependent f){
		return ch<<std::dynamic_pointer_cast<ITrackParamProcess>(std::make_shared<Condition>(f));
	}
	class Chain:public AbstractChain{
	public:
		Chain(){}
		virtual ~Chain(){}
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	};
	class ChainCheck:public AbstractChain{
	public:
		ChainCheck(){}
		virtual ~ChainCheck(){}
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	};
	class ChainAnd:public AbstractChain{
	public:
		ChainAnd(){}
		virtual ~ChainAnd(){}
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	};
	class ChainOr:public AbstractChain{
	public:
		ChainOr(){}
		virtual ~ChainOr(){}
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	};
	class ChainBinner:public AbstractChain{
	public:
		ChainBinner(const Axis&source);
		ChainBinner(Axis&&source);
		virtual ~ChainBinner();
		virtual bool Process(WTrack&,std::vector<double>&)const override;
	private:
		Axis m_axis;
	};
	class TrackProcess{
	public:
		TrackProcess(){}
		~TrackProcess(){}
		TrackProcess&operator<<(std::shared_ptr<ITrackParamProcess>element);
		void Process(WTrack&T)const;
	private:
		std::vector<std::shared_ptr<ITrackParamProcess>> m_proc;
	};
	std::shared_ptr<TrackProcess>operator<<(std::shared_ptr<TrackProcess>ch,std::shared_ptr<ITrackParamProcess>v);
	inline TrackProcess&operator<<(TrackProcess&ch,ConditionTrackParamDependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	inline TrackProcess&operator<<(TrackProcess&ch,ConditionTrackDependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	inline TrackProcess&operator<<(TrackProcess&ch,ConditionParamDependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	inline TrackProcess&operator<<(TrackProcess&ch,ConditionIndependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	class EventProcess{
	public:
		EventProcess(){}
		~EventProcess(){}
		EventProcess&operator<<(std::shared_ptr<ITrackParamProcess>element);
		bool Process()const;
	private:
		std::vector<std::shared_ptr<ITrackParamProcess>> m_proc;
	};
	std::shared_ptr<EventProcess>operator<<(std::shared_ptr<EventProcess>ch,std::shared_ptr<ITrackParamProcess>v);
	inline EventProcess&operator<<(EventProcess&ch,ConditionTrackParamDependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	inline EventProcess&operator<<(EventProcess&ch,ConditionTrackDependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	inline EventProcess&operator<<(EventProcess&ch,ConditionParamDependent f){
		return ch<<std::make_shared<Condition>(f);
	}
	inline EventProcess&operator<<(EventProcess&ch,ConditionIndependent f){
		return ch<<std::make_shared<Condition>(f);
	}
}
#endif 
