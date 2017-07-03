// this file is distributed under
// GPL license
#include <math_h/vectors.h>
#include <math_h/randomfunc.h>
#include <math_h/functions.h>
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
using namespace MathTemplates;
const EventGenerator Simulation4Test(mt19937&RG){
	return [&RG]()->list<particle_sim>{
		const auto g1=Vector4<double>::SpaceLength4(Vector3<double>::RandomIsotropicDirection(RG)*Particle::eta().mass()/2.0,0.0);
		const auto g2=Vector4<double>::SpaceLength4(-g1.space_component(),0.0);
		return{
			{.type=Particle::gamma(),.P=g1.space_component()},
			{.type=Particle::gamma(),.P=g2.space_component()}
		};
	};
}
int main(){
	mt19937 RG;
	Simulate("test",Simulation4Test(RG));
	return 0;
}
