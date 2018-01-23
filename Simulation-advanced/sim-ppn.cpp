// this file is distributed under
// GPL license
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
const BiSortedPoints<> ReadCrossSection(){
    BiSortedPoints<> result(ChainWithCount(91, 0., PI<>()/2.0), ChainWithCount(13, 1.000, 2.200));
    for (size_t degree = 0; degree <= 90; degree++) {
        ifstream file("pp/Theta_" + to_string(degree) + ".txt");
        for (double E = 0, C = 0; (file >> E >> C); result.Bin(degree, (size_t(E) - 1000) / 100) = C);
        file.close();
    }
    return result;
}
const vector<RandomValueTableDistr<>> AngularDistribution(const BiLinearInterpolation<>&source){
	vector<RandomValueTableDistr<>> res;
	for(double p=1.0;p<2.1;p+=0.001){
		SortedPoints<> chain;
		for(const auto&angle:source.X())chain<<make_point(angle,source(angle,p));
		if(size_t(p*1000.)%100==0)Plot().Line(chain)<<"set title'P="+to_string(p)+", cross section"<<"set log y";
		chain*=[](const double&th){return sin(th);};
		if(size_t(p*1000.)%100==0)Plot().Line(chain)<<"set title'P="+to_string(p)+", theta distribution density'"<<"unset log y";
		res.push_back(RandomValueTableDistr<>(chain));
	}
	return res;
}
int main(){
	Plotter::Instance().SetOutput(".","sim-ppn");
	const RandomUniform<>Pb_distr(p_beam_low,p_beam_hi);
	const auto Pf_dens=Plotter::Instance().GetPoints<double>("pp/pfermi");
	Plot().Line(Pf_dens);
	const RandomValueTableDistr<>Pf_distr=Pf_dens;
	const auto CS=AngularDistribution(ReadCrossSection());
	const auto THETA=[&CS](double p)->double{
		int index=size_t(((p-1.)*100.)+0.5);
		static PlotDistr1D<> p_dep("Pp(pt) index","",BinsByCount(101,-0.5,100.5));
                p_dep.Fill(index);
		if(index<0){
			return CS[0]();
		}else{
			if(size_t(index)<CS.size())
				return CS[index]();
			else
				return CS[CS.size()-1]();
		}
	};
	Simulate("ppn_qf_",[&Pb_distr,&Pf_distr,&THETA]()->list<particle_sim>{
                const auto d_lab=lorentz_byPM(Zero(),Particle::d().mass());
                const auto nt_lab=lorentz_byPM(randomIsotropic<3>()*Pf_distr(),Particle::n().mass());
                const auto pt_lab=d_lab-nt_lab;
		static PlotDistr1D<> p_mass("p_t mass","",BinsByCount(100,0.8,1.0));
		p_mass.Fill(pt_lab.M());
                const auto pr_lab=lorentz_byPM(Z()*Pb_distr(),Particle::p().mass());
                const auto PP=pr_lab+pt_lab;
		const auto ppr_t =pr_lab.Transform(pt_lab.Beta());
		const auto ppr_cm=pr_lab.Transform(PP.Beta());
		static PlotDistr1D<> p_dep_x("Pp(pt)_x","",BinsByCount(100,-0.5,0.5));
		static PlotDistr1D<> p_dep_y("Pp(pt)_y","",BinsByCount(100,-0.5,0.5));
		static PlotDistr1D<> p_dep_z("Pp(pt)_z","",BinsByCount(100,1.0,2.0));
		static PlotDistr1D<> IM_plot_before("IM before","",BinsByCount(500,1.8,2.8));
		static PlotDistr1D<> IM_plot_after("IM","",BinsByCount(500,1.8,2.8));
		IM_plot_before.Fill(PP.M());
		if(PP.M()<=(2.0*Particle::p().mass()))return {};
		IM_plot_after.Fill(PP.M());
		p_dep_x.Fill(ppr_t.P().x());
		p_dep_y.Fill(ppr_t.P().y());
		p_dep_z.Fill(ppr_t.P().z());
                const auto final_chm=binaryDecay(PP.M(),Particle::p().mass(),Particle::p().mass(),direction(ppr_cm.P()));
                const static RandomUniform<> PHI(-PI(),PI());
		const auto theta_cm=THETA(ppr_t.P().M());
		const auto phi_cm=PHI();
		const auto transform=[&theta_cm,&phi_cm,&ppr_cm](const LorentzVector<>&P)->const LorentzVector<>{
			return P.Rotate(direction(ppr_cm.P()^X()),theta_cm).Rotate(direction(ppr_cm.P()),phi_cm);
		};
		const auto final_cm=make_pair(transform(final_chm.first),transform(final_chm.second));
                const auto p1_lab=final_cm.first.Transform(-PP.Beta());
                const auto p2_lab=final_cm.second.Transform(-PP.Beta());
                return {
                        {.type=Particle::n(),.P=nt_lab.P()},
                        {.type=Particle::p(),.P=p1_lab.P()},
                        {.type=Particle::p(),.P=p2_lab.P()}
                };
        },40);
	return 0;
}
