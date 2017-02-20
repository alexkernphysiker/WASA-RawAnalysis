// this file is distributed under 
// GPL license
#ifndef mpOBjkfX
#define mpOBjkfX
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <fstream>
#include <math_h/error.h>
#include <math_h/interpolate.h>
#include <Genetic/paramfunc.h>
#include <config.h>
#include "analysis.h"
const std::string DataFiles=std::string("./")+DataDirectoryName+"/";
template<typename... Args>
class InterpolationBasedReconstruction{
public:
	typedef std::function<double(Args...)> FUNC;
private:
	std::string m_name;
	MathTemplates::LinearInterpolation<double> data;
	std::shared_ptr<std::vector<MathTemplates::point<double>>> out;
	FUNC Experiment,Theory;
public:
	InterpolationBasedReconstruction(const std::string&&name,const FUNC measured,const FUNC theory){
		out=std::make_shared<std::vector<MathTemplates::point<double,double>>>();
		m_name=name;
		Experiment=measured;
		Theory=theory;
		std::ifstream file;
		file.open((DataFiles+name+".calibration.txt").c_str());
		if(file){
			double measured,calculated;
			while(file>>measured>>calculated)
				data<<MathTemplates::point<double>(measured,calculated);
			file.close();
		}
	}
	InterpolationBasedReconstruction(
		const InterpolationBasedReconstruction&source
	){
		m_name=source.m_name;
		data=source.data;
		out=source.out;
		Experiment=source.Experiment;
		Theory=source.Theory;
	}
	virtual ~InterpolationBasedReconstruction(){
		if(out->size()>0){
			std::ofstream file;
			file.open((DataFiles+m_name+".simulation.txt").c_str(),ios_base::app);
			if(file){
				for(auto&p:*out)
					file<<p.X()<<" "<<p.Y()<<std::endl;
				file.close();
			}
		}
	}
	double operator()(Args... args)const{
		if(data.size()>0){
			try{
				return data(Experiment(args...));
			}catch(exception){
				return INFINITY;
			}
		}else{
			out->push_back(MathTemplates::point<double>(Experiment(args...),Theory(args...)));
			return INFINITY;
		}
	}
	FUNC func()const{return [this](Args... args)->double{return operator()(args...);};}
};

template<class FitFunc,typename... Args>
class FitBasedReconstruction{
public:
	typedef std::function<double(Args...)> FUNC;
private:
	std::string m_name;
	FitFunc m_func;
	FUNC Theory;
	std::vector<FUNC> Experiment;
	Genetic::ParamSet P;
	std::shared_ptr<std::vector<Genetic::ParamSet>> data;
public:
	FitBasedReconstruction(const std::string&&name,const std::vector<FUNC>&&measured,const FUNC theory){
		using namespace Genetic;
		m_name=name;
		Theory=theory;
		for(auto f:measured)Experiment.push_back(f);
		data=std::make_shared<std::vector<Genetic::ParamSet>>();
		std::ifstream file;
		file.open((DataFiles+name+".fit.txt").c_str());
		if(file){
			file>>P;
			file.close();
		}
	}
	FitBasedReconstruction(const FitBasedReconstruction&source){
		m_name=source.m_name;
		for(const auto f:source.Experiment)Experiment.push_back(f);
		Theory=source.Theory;
		P=source.P;
		data=source.data;
	}
	virtual ~FitBasedReconstruction(){
		using namespace Genetic;
		if(data->size()>0){
			std::ofstream file;
			file.open((DataFiles+m_name+".simulation.txt").c_str(),ios_base::app);
			if(file){
				for(const ParamSet&p:*data)file<<p<<std::endl;
				file.close();
			}
		}
	}
	double operator()(Args... args)const{
		using namespace Genetic;
		ParamSet X;
		for(FUNC f:Experiment)X<<f(args...);
		if(P.size()>0){
			return m_func(X,P);
		}else{
			X<<Theory(args...);
			data->push_back(X);			
			return INFINITY;
		}
	}
	FUNC func(){return [this](Args... args)->double{return operator()(args...);};}
};
#endif 
