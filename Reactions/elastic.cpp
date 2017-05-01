// this file is distributed under
// GPL license
#include <Wasa.hh>
#include <TCutG.h>
#include <math_h/functions.h>
#include <math_h/tabledata.h>
#include <math_h/error.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include <trackprocessing.h>
#include <detectors.h>
#include <reconstruction.h>
#include <data.h>
#include "elastic.h"
using namespace std;
using namespace MathTemplates;
using namespace TrackAnalyse;

void ReconstructPD(Analysis&res){
	res.Trigger(trigger_he3_forward.number).pre()
	;
	static const string recdmsg="ForwardD_Reconstruction";
	
        res.Trigger(trigger_he3_forward.number).per_track()
		<<(make_shared<ChainCheck>()
                	<<[](WTrack&T)->bool{return T.Type()==kFDC;}
                	<<Forward::Get().CreateMarker(recdmsg,"1-ChargedTracks")
        	)
	;
}
