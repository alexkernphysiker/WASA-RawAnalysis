// this file is distributed under
// GPL license
#include <Experiment/experiment_conv.h>
#include "runpluto.h"
int main(){
    SimulatePluto("He3 eta [pi0 [g g] pi0 [g g] pi0 [g g]]","He3eta6g",1.573,p_beam_hi);
    SimulatePluto("He3 pi0 [g g]  pi0 [g g] pi0 [g g]","He3pi06g",p_beam_low,p_beam_hi);
    return 0;
}
