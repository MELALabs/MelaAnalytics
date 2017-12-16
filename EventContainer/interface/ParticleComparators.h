#ifndef PARTICLECOMPARATORS_H
#define PARTICLECOMPARATORS_H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

namespace ParticleComparators{
  extern double jetDeltaR;
  extern double jetEtaAcceptanceCut;
  extern double jetPTCut;
  extern double electronEtaAcceptanceCut;
  extern double electronPTCut;
  extern double muonEtaAcceptanceCut;
  extern double muonPTCut;
  extern double ghostDeltaRCut;

  extern double mV1LowCut;
  extern double mV2LowCut;
  extern double mllLowCut;
  extern double mV12HighCut;

  void setJetDeltaR(double cut);
  void setJetEtaAcceptanceCut(double cut);
  void setJetPTCut(double cut);
  void setElectronEtaAcceptanceCut(double cut);
  void setElectronPTCut(double cut);
  void setMuonEtaAcceptanceCut(double cut);
  void setMuonPTCut(double cut);
  void setGhostDeltaRCut(double cut);

  void setMV1LowCut(double cut);
  void setMV2LowCut(double cut);
  void setMllLowCut(double cut);
  void setMV12HighCut(double cut);
}

#endif
