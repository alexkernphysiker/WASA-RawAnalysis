// this file is distributed under
// GPL license
#include <Experiment/experiment_conv.h>
#include "runsim.h"
using namespace std;
int main(){
    Simulate("bound-2g",[](mt19937&RG)->list<particle_sim>{
	return {};
    });
    return 0;
}
