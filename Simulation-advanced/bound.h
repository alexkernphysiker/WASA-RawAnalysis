// this file is distributed under
// GPL license
#include <math_h/vectors.h>
#include <math_h/randomfunc.h>
#include <math_h/functions.h>
#include <math_h/interpolate.h>
#include "runsim.h"
const MathTemplates::LinearInterpolation<> ReadPfFromFile(const std::string&name);
struct etamesic{MathTemplates::LorentzVector<>he3;MathTemplates::LorentzVector<>eta_;};
const etamesic Compound(
	MathTemplates::RANDOM&RG,
	const MathTemplates::RandomValueGenerator<>&Pb_distr,
	const MathTemplates::RandomValueGenerator<>&Pf_distr,
	const double&s_thr=0
);
