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
		const auto ppr_t =pr_lab.Transform(pt_lab.Beta());
		const auto ppr_cm=pr_lab.Transform(PP.Beta());
		static PlotDistr1D<> p_dep_x("Pp(pt)_x","",BinsByCount(100,-0.5,0.5));
		static PlotDistr1D<> p_dep_y("Pp(pt)_y","",BinsByCount(100,-0.5,0.5));
		static PlotDistr1D<> p_dep_z("Pp(pt)_z","",BinsByCount(100,1.0,2.0));
                const static RandomUniform<> PHI(-PI(),PI());
		if(PP.M()<=(2.0*Particle::p().mass()))return {};
		p_dep_x.Fill(ppr_t.S().x());
		p_dep_y.Fill(ppr_t.S().y());
		p_dep_z.Fill(ppr_t.S().z());
                const auto final_chm=binaryDecay(PP.M(),Particle::p().mass(),Particle::p().mass(),Angles(ppr_cm.S()));
		const auto transform=[&ppr_cm,&ppr_t,&THETA,&RG](const LorentzVector<>&P)->const LorentzVector<>{
			return P.Rotate(ppr_cm.S()^X<>(),THETA(RG,ppr_t.S().M())).Rotate(ppr_cm.S(),PHI(RG));
		};
		const auto final_cm=make_pair(transform(final_chm.first),transform(final_chm.second));
                const auto p1_lab=final_cm.first.Transform(-PP.Beta());
                const auto p2_lab=final_cm.second.Transform(-PP.Beta());
                return {
                        {.type=Particle::n(),.P=nt_lab.S()},
                        {.type=Particle::p(),.P=p1_lab.S()},
                        {.type=Particle::p(),.P=p2_lab.S()}
                };
        });
	return 0;
}
