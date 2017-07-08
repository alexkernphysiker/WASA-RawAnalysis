// this file is distributed under
// GPL license
#include <math_h/vectors.h>
#include <math_h/randomfunc.h>
#include <math_h/functions.h>
#include "runsim.h"
const MathTemplates::SortedPoints<> ReadPfFromFile(const std::string&name);
const std::pair<MathTemplates::Vector4<>,MathTemplates::Vector4<>> Compound(
	MathTemplates::RANDOM&RG,
	const MathTemplates::IFunction<double,MathTemplates::RANDOM&>&Pb_distr,
	const MathTemplates::IFunction<double,MathTemplates::RANDOM&>&Pf_distr,
	const double&s_thr=0
);
