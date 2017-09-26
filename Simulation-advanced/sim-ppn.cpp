// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const BiSortedPoints<> ReadCrossSection(){
    BiSortedPoints<> result(ChainWithCount(181, 0., PI<>()), ChainWithCount(13, 1.000, 2.200));
    for (size_t degree = 0; degree <= 180; degree++) {
        ifstream file("pp/Theta_" + to_string(degree) + ".txt");
        for (double E = 0, C = 0; (file >> E >> C); result.Bin(degree, (size_t(E) - 1000) / 100) = C);
        file.close();
    }
    return result;
}
const vector<RandomValueTableDistr<>> PrepareCrossSection(const BiLinearInterpolation<>&source){
	vector<RandomValueTableDistr<>> res;
	for(double p=1.0;p<2.1;p+=0.001){
		SortedPoints<> chain;
		for(const auto&angle:source.X())chain<<make_point(angle,source(angle,p));
		res.push_back(RandomValueTableDistr<>(chain*[](const double&th){return sin(th);} ));
	}
	return res;
}
int main(){
	RANDOM RG;
	Plotter<>::Instance().SetOutput(".","sim-ppn");
	const RandomUniform<>Pb_distr(p_beam_low,p_beam_hi);
	const auto Pf_dens=Plotter<>::Instance().GetPoints<double>("pp/pfermi");
	Plot<>().Line(Pf_dens);
	const RandomValueTableDistr<>Pf_distr=Pf_dens;
	const auto CS=PrepareCrossSection(ReadCrossSection());
	const auto THETA=[&CS](RANDOM&RG,double p)->double{
		size_t index=size_t(((p-1.0)*100.0)+0.5);
		static PlotDistr1D<> p_dep("Pp(pt) index","",BinsByCount(101,-0.5,100.5));
                p_dep.Fill(index);
		if(index<CS.size())
			return CS[index](RG);
		else
			return 0;
	};
	Simulate("ppn_qf_",[&RG,&Pb_distr,&Pf_distr,&THETA]()->list<particle_sim>{
                const LorentzVector<> d_lab=Particle::d().mass();
                const auto nt_lab=lorentz_byPM(RandomIsotropicDirection3<>(RG)*Pf_distr(RG),Particle::n().mass());
                const auto pt_lab=d_lab-nt_lab;
                const auto pr_lab=lorentz_byPM(Z<>()*Pb_distr(RG),Particle::p().mass());
                const auto PP=pr_lab+pt_lab;
                const double ppr_dep=pr_lab.Lorentz(pt_lab.Beta()).space_component().mag();
                static PlotDistr1D<> p_dep("Pp(pt)","",BinsByCount(100,1.0,2.0));
                p_dep.Fill(ppr_dep);
                const static RandomUniform<> PHI(0.0,2.0*PI<>());
		if(PP.length4()<=(2.0*Particle::p().mass()))return {};
                const auto final_cm=binaryDecay(PP.length4(),Particle::p().mass(),Particle::p().mass(),THETA(RG,ppr_dep),PHI(RG));
                const auto p1_lab=final_cm.first.Lorentz(-PP.Beta());
                const auto p2_lab=final_cm.second.Lorentz(-PP.Beta());
                return {
                        {.type=Particle::n(),.P=nt_lab.space_component()},
                        {.type=Particle::p(),.P=p1_lab.space_component()},
                        {.type=Particle::p(),.P=p2_lab.space_component()}
                };
        });
	return 0;
}
