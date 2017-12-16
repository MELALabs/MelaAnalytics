#include "ParticleComparators.h"


namespace ParticleComparators{
  double jetDeltaR = 0.5;
  double jetEtaAcceptanceCut = 4.7;
  double jetPTCut = 30.;
  double electronEtaAcceptanceCut = 2.5;
  double electronPTCut = 7.;
  double muonEtaAcceptanceCut = 2.4;
  double muonPTCut = 5.;
  double ghostDeltaRCut = 0.02;
  double mV1LowCut = 40.;
  double mV2LowCut = 12.;
  double mllLowCut = 4.;
  double mV12HighCut = 120.;
}

void ParticleComparators::setJetDeltaR(double cut){ jetDeltaR = cut; }
void ParticleComparators::setJetEtaAcceptanceCut(double cut){ jetEtaAcceptanceCut = cut; }
void ParticleComparators::setJetPTCut(double cut){ jetPTCut = cut; }
void ParticleComparators::setElectronEtaAcceptanceCut(double cut){ electronEtaAcceptanceCut = cut; }
void ParticleComparators::setElectronPTCut(double cut){ electronPTCut = cut; }
void ParticleComparators::setMuonEtaAcceptanceCut(double cut){ muonEtaAcceptanceCut = cut; }
void ParticleComparators::setMuonPTCut(double cut){ muonPTCut = cut; }
void ParticleComparators::setGhostDeltaRCut(double cut){ ghostDeltaRCut = cut; }

void ParticleComparators::setMV1LowCut(double cut){ mV1LowCut = cut; }
void ParticleComparators::setMV2LowCut(double cut){ mV2LowCut = cut; }
void ParticleComparators::setMllLowCut(double cut){ mllLowCut = cut; }
void ParticleComparators::setMV12HighCut(double cut){ mV12HighCut = cut; }

