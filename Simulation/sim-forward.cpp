// this file is distributed under
// GPL license
#include <Experiment/experiment_conv.h>
#include "runpluto.h"
int main(){
    SimulatePluto("He3 pi0","He3pi0",p_beam_low,p_beam_hi);
    SimulatePluto("He3 pi0 pi0","He3pi0pi0",p_beam_low,p_beam_hi);
    SimulatePluto("He3 pi0 pi0 pi0","He3pi0pi0pi0",p_beam_low,p_beam_hi);
    return 0;
}
