#ifndef MELACANDIDATERECASTER_H
#define MELACANDIDATERECASTER_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "Mela.h"


class MELACandidateRecaster{
public:
  static double getVffEquivalentCoupling(int iferm, int jferm);
  static void readCandidate(
    MELACandidate* cand,
    SimpleParticleCollection_t& mothers,
    SimpleParticleCollection_t& daughters,
    SimpleParticleCollection_t& associated
    );
  static MELAParticle* getBestAssociatedV(MELACandidate* cand, TVar::Production production, double* bestScore=nullptr);
  static void adjustForIncomingMomenta(
    SimpleParticleCollection_t& mothers,
    SimpleParticleCollection_t& daughters,
    SimpleParticleCollection_t& associated
    );

protected:
  constexpr static const double GeVunit=1e-2;
  constexpr static const double GeVsqunit=1e-4;

  TVar::Production candScheme;
  bool doDotincomingParticles;
  bool protectVStrict;
  unsigned int nAssociated;

  std::vector<MELAParticle*> extraParticles;

  double getQGMergeScore(MELAParticle* p1, MELAParticle* p2) const;
  double getGGMergeScore(MELAParticle* g1, MELAParticle* g2) const;
  double getMergeScore(MELAParticle* p1, MELAParticle* p2) const;

  MELAParticle* mergeTwoGluons(MELAParticle* glu1, MELAParticle* glu2);

  bool merge2Qto1G(
    MELACandidate*& cand,
    std::vector<MELAParticle*>& gluons,
    std::vector<MELAParticle*>& quarks,
    MELAParticle*& protectV
    );

  std::vector<std::pair<MELAParticle*, std::vector<MELAParticle*>>> mapGluonsToQuarks(
    const std::vector<MELAParticle*>& gluons,
    const std::vector<MELAParticle*>& quarks,
    const std::vector<int>& permutation
    );

  double getBestVHConfig(
    MELAParticle* protectV,
    const std::vector<std::pair<MELAParticle*, std::vector<MELAParticle*>>>& quarks,
    std::vector<int>* qordered=nullptr,
    int* swapconfig=nullptr
    );

  double getBestVBFConfig(
    const std::vector<std::pair<MELAParticle*, std::vector<MELAParticle*>>>& quarks,
    std::vector<int>* qordered=nullptr,
    int* swapconfig=nullptr
    );

  double getBestHJJConfig(
    const std::vector<std::pair<MELAParticle*, std::vector<MELAParticle*>>>& partons,
    std::vector<int>* qordered=nullptr,
    int* swapconfig=nullptr
  );

  double getVBFLikelihood(
    const std::vector<MELAParticle*>& gluons,
    const std::vector<MELAParticle*>& quarks,
    const std::vector<int>& permutation
    );

  double getMergeOrder_GluonsIntoQuarks(
    std::vector<MELAParticle*>& gluons,
    std::vector<MELAParticle*>& quarks,
    std::vector<int>& bestConfig,
    bool doMergeGluons,
    std::vector<MELAParticle*>* VmassMonitor=nullptr
    );

  MELAParticle* getProtectedV(MELACandidate* cand, double* bestScore=nullptr){ return getBestAssociatedV(cand, candScheme, bestScore); }

  void clearExtraParticles(){ for (auto& part:extraParticles) delete part; }

public:

  MELACandidateRecaster(const TVar::Production candScheme_);
  MELACandidateRecaster() : candScheme(TVar::nProductions) {}
  // Delete the copy and move constructors
  MELACandidateRecaster(const MELACandidateRecaster& other) = delete;
  MELACandidateRecaster(const MELACandidateRecaster&& other) = delete;

  ~MELACandidateRecaster();

  void copyCandidate(MELACandidate* cand, MELACandidate*& candModified, bool adjustIncoming=false);

  void reduceJJtoQuarks(MELACandidate*& cand);

  void deduceLOVHTopology(MELACandidate*& cand);
  void deduceLOVBFTopology(MELACandidate*& cand); // Alias for reduceJJtoQuarks
  void deduceLOHJJTopology_JetMerge(MELACandidate*& cand);
  void deduceLOHJJTopology_LeadingPt(MELACandidate*& cand);

};


#endif
