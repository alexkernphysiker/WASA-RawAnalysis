// this file is distributed under 
// GPL license
#include <map>
#include "math_h/vectors.h"
#include "Kinematics/particles.h"
#include "montecarlo.h"
using namespace std;
using namespace MathTemplates;
enum ParticleType{
	kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,
	kNeutron=13,kProton=14,kEta=17,kDeuteron=45,kTriton=46,kHe3=49,
	kDDummy=50
};
map<ParticleType,Particle> dictionary;
MonteCarlo::MonteCarlo(const bool kostyl):Analysis(),f_kostyl(kostyl){
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
	dictionary[kGamma]=Particle::gamma();
	dictionary[kPi0]=Particle::pi0();
	dictionary[kPiMinus]=Particle::pi_minus();
	dictionary[kPiPlus]=Particle::pi_plus();
	dictionary[kNeutron]=Particle::n();
	dictionary[kProton]=Particle::p();
	dictionary[kEta]=Particle::eta();
	dictionary[kDeuteron]=Particle::d();
	dictionary[kHe3]=Particle::he3();
}
MonteCarlo::~MonteCarlo(){}
bool MonteCarlo::DataTypeSpecificEventAnalysis()const{
	if(gWasa->IsAnalysisMode(Wasa::kMCRaw)||
		gWasa->IsAnalysisMode(Wasa::kMCReco)||
		gWasa->IsAnalysisMode(Wasa::kMC)
	){
		ClearVerticesCache();
		WVertexIter iterator(fMCVertexBank);
		if((!f_kostyl)||(dynamic_cast<WVertex*>(iterator.Next())))
			while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
				vector<Kinematic> vertex_rec;
				for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
					WParticle *p=vertex->GetParticle(particleindex);
					if((p->GetType()>kDummy)&&(p->GetType()<kDDummy))
						vertex_rec.push_back({
							.particle=dictionary[ParticleType(p->GetType())],
							.E=p->GetEkin(),
							.Th=p->GetTheta(),
							.Phi=p->GetPhi()
						});
				}
				CacheVertex(vertex_rec);
			}
		auto P=Vector3<double>::zero();
		for(const Kinematic&item:Vertex(0))
			P+=Vector3<double>::Polar(item.particle.E2P(item.E),item.Th,item.Phi);
		CachePBeam(P.mag());
		return true;
	}
	return false;
}
bool MonteCarlo::DataSpecificTriggerCheck(int n)const{
	return true;
}
