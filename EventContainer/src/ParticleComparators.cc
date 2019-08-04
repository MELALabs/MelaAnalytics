#include "ParticleComparators.h"


namespace ParticleComparators{
  double jetDeltaR = 0.4;
  double jetEtaCut = 4.7;
  double jetPtCut = 30.;
  double electronEtaCut = 2.5;
  double electronPtCut = 7.;
  double muonEtaCut = 2.4;
  double muonPtCut = 5.;
  double ghostDeltaRCut = 0.02;
  double mV1LowCut = 40.;
  double mV2LowCut = 12.;
  double mllLowCut = 4.;
  double mV1HighCut = 120.;
  double mV2HighCut = 120.;
  double mllHighCut = 120.;
}

void ParticleComparators::setJetDeltaR(double cut){ jetDeltaR = cut; }
void ParticleComparators::setJetEtaCut(double cut){ jetEtaCut = cut; }
void ParticleComparators::setJetPtCut(double cut){ jetPtCut = cut; }
void ParticleComparators::setElectronEtaCut(double cut){ electronEtaCut = cut; }
void ParticleComparators::setElectronPtCut(double cut){ electronPtCut = cut; }
void ParticleComparators::setMuonEtaCut(double cut){ muonEtaCut = cut; }
void ParticleComparators::setMuonPtCut(double cut){ muonPtCut = cut; }
void ParticleComparators::setGhostDeltaRCut(double cut){ ghostDeltaRCut = cut; }

void ParticleComparators::setMV1LowCut(double cut){ mV1LowCut = cut; }
void ParticleComparators::setMV2LowCut(double cut){ mV2LowCut = cut; }
void ParticleComparators::setMllLowCut(double cut){ mllLowCut = cut; }
void ParticleComparators::setMV1HighCut(double cut){ mV1HighCut = cut; }
void ParticleComparators::setMV2HighCut(double cut){ mV2HighCut = cut; }
void ParticleComparators::setMllHighCut(double cut){ mllHighCut = cut; }

