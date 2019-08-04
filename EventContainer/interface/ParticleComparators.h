#ifndef PARTICLECOMPARATORS_H
#define PARTICLECOMPARATORS_H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

namespace ParticleComparators{
  extern double jetDeltaR;
  extern double jetEtaCut;
  extern double jetPtCut;
  extern double electronEtaCut;
  extern double electronPtCut;
  extern double muonEtaCut;
  extern double muonPtCut;
  extern double ghostDeltaRCut;

  extern double mV1LowCut;
  extern double mV2LowCut;
  extern double mllLowCut;
  extern double mV1HighCut;
  extern double mV2HighCut;
  extern double mllHighCut;

  void setJetDeltaR(double cut);
  void setJetEtaCut(double cut);
  void setJetPtCut(double cut);
  void setElectronEtaCut(double cut);
  void setElectronPtCut(double cut);
  void setMuonEtaCut(double cut);
  void setMuonPtCut(double cut);
  void setGhostDeltaRCut(double cut);

  void setMV1LowCut(double cut);
  void setMV2LowCut(double cut);
  void setMllLowCut(double cut);
  void setMV1HighCut(double cut);
  void setMV2HighCut(double cut);
  void setMllHighCut(double cut);
}

#endif
