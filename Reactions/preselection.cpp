// this file is distributed under 
// GPL license
#include <iostream>
#include <Wasa.hh>
#include <Experiment/experiment_conv.h>
#include <data.h>
#include "analyses.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;
void Preselection(Analysis&res){
	Int_t RunNumber = SorterOption::GetIntValue("RunNumber");
	//ACHTUNG: HARDCODE!!!!!
	TString Tstream1 = Form("file:/home/PRESEL/run_%i.bz2", RunNumber);
	gWasa->AddOutput("run_presel",Tstream1);

        //Preselection for 3He X reactions
        static long he3count;
        res.Trigger(trigger_he3_forward.number).pre()<<[](){he3count=0;return true;};
        static particle_kine he3;
        res.Trigger(trigger_he3_forward.number).per_track()<<(make_shared<ChainCheck>()
	    <<ForwardHe3Reconstruction("He3Forward",res,he3)
	    <<[](){he3count++;return true;}
        );
        res.Trigger(trigger_he3_forward.number).post()<<[](){
            if(he3count>0)gWasa->SaveEvent("run_presel");
            return true;
        };

        //Preselection for quasielastic scattering
        static long c_tracks;
	res.Trigger(trigger_elastic1.number).pre()<<[](){c_tracks=0;return true;};
	res.Trigger(trigger_elastic1.number).per_track()<<(make_shared<ChainCheck>()
                <<[](WTrack&T){return T.Type()==kCDC;}
                <<[](WTrack&T){return T.Theta()*180./PI()>23.;}
                <<[](WTrack&T){return T.Edep()>0.03;}
                <<[](WTrack&T){c_tracks++;return true;}
	);
        res.Trigger(trigger_elastic1.number).post()<<[](){
                if(c_tracks>=2)gWasa->SaveEvent("run_presel");
                return true;
        };
}
