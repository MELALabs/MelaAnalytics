#include "HiggsComparators.h"
#include "MELAStreamHelpers.hh"


MELACandidate* HiggsComparators::matchAHiggsToParticle(MELAEvent& ev, MELAParticle* genH){
  MELACandidate* cand=nullptr;
  if (!genH) return cand;
  for (MELACandidate* tmpCand:ev.getCandidates()){
    if (!tmpCand) continue;
    double genmassquant = genH->m()+genH->pt()+fabs(genH->z());
    double massdiff = fabs(genmassquant-tmpCand->m()-tmpCand->pt()-fabs(tmpCand->z()));
    double massratio = 0;
    if (genmassquant>0.) massratio = massdiff / genmassquant;
    if (massratio<0.001){
      if (!cand) cand = tmpCand;
      else{
        TLorentzVector const& vGen = genH->p4;
        TLorentzVector const& vTmp = tmpCand->p4;
        TLorentzVector const& vCur = cand->p4;

        double dot_tmp = vTmp.Dot(vGen);
        double dot_curr = vCur.Dot(vGen);
        if (fabs(dot_tmp-vGen.M2())<fabs(dot_curr - vGen.M2())) cand = tmpCand;
      }
    }
  }
  return cand;
}

MELACandidate* HiggsComparators::candidateSelector(MELAEvent& ev, HiggsComparators::CandidateSelection scheme, MELAEvent::CandidateVVMode VVMode){
  MELACandidate* cand=nullptr;
  for (MELACandidate* tmpCand:ev.getCandidates()){
    if (!tmpCand) continue;
    if (!tmpCand->passSelection) continue;
    if (!cand) cand=tmpCand;
    else cand = HiggsComparators::candComparator(cand, tmpCand, scheme, VVMode);
  }
  return cand;
}

MELACandidate* HiggsComparators::candComparator(MELACandidate* cand1, MELACandidate* cand2, HiggsComparators::CandidateSelection scheme, MELAEvent::CandidateVVMode VVMode){
  MELACandidate* theChosenOne=nullptr;

  if (!cand1 && !cand2) return theChosenOne;
  else if (cand1 && !cand2) return cand1;
  else if (!cand1 && cand2) return cand2;

  TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
  double HVVmass = PDGHelpers::Zeromass;
  if (VVMode==MELAEvent::WWMode){
    PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_WW);
    HVVmass = PDGHelpers::Wmass;
  }
  else if (VVMode==MELAEvent::ZZMode){
    PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    HVVmass = PDGHelpers::Zmass;
  }
  else if (VVMode==MELAEvent::ZGammaMode){
    PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZG);
    HVVmass = PDGHelpers::Zmass;
  }
  else if (VVMode==MELAEvent::GammaGammaMode) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_GG);
  else PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ff);

  if (VVMode==MELAEvent::UndecayedMode){
    if (
      (
        (cand2->m()>cand1->m())
        ||
        (cand2->m()==cand1->m() && cand2->pt()>cand1->pt())
      ) && PDGHelpers::isAHiggs(cand2->id)
      ) theChosenOne = cand2;
    else if (PDGHelpers::isAHiggs(cand1->id)) theChosenOne = cand1;
  }
  else if (scheme==HiggsComparators::BestZ1ThenZ2ScSumPt){
    double diffmass1 = fabs(cand1->getSortedV(0)->m()-HVVmass);
    double diffmass2 = fabs(cand2->getSortedV(0)->m()-HVVmass);

    double Z1scsumpt_cand1=0, Z1scsumpt_cand2=0;
    MELAParticle* c1_11 = cand1->getSortedV(0)->getDaughter(0);
    MELAParticle* c1_12 = cand1->getSortedV(0)->getDaughter(1);
    MELAParticle* c2_11 = cand2->getSortedV(0)->getDaughter(0);
    MELAParticle* c2_12 = cand2->getSortedV(0)->getDaughter(1);
    if (c1_11) Z1scsumpt_cand1 += c1_11->pt();
    if (c1_12) Z1scsumpt_cand1 += c1_12->pt();
    if (c2_11) Z1scsumpt_cand2 += c2_11->pt();
    if (c2_12) Z1scsumpt_cand2 += c2_12->pt();

    double Z2scsumpt_cand1=0, Z2scsumpt_cand2=0;
    MELAParticle* c1_21 = cand1->getSortedV(1)->getDaughter(0);
    MELAParticle* c1_22 = cand1->getSortedV(1)->getDaughter(1);
    MELAParticle* c2_21 = cand2->getSortedV(1)->getDaughter(0);
    MELAParticle* c2_22 = cand2->getSortedV(1)->getDaughter(1);
    if (c1_21) Z2scsumpt_cand1 += c1_21->pt();
    if (c1_22) Z2scsumpt_cand1 += c1_22->pt();
    if (c2_21) Z2scsumpt_cand2 += c2_21->pt();
    if (c2_22) Z2scsumpt_cand2 += c2_22->pt();

    if (
      diffmass1>diffmass2
      ||
      (
        diffmass1==diffmass2 && (
          ((VVMode==MELAEvent::WWMode || VVMode==MELAEvent::ZZMode) && Z2scsumpt_cand2>Z2scsumpt_cand1)
          ||
          (Z1scsumpt_cand2>Z1scsumpt_cand1 || (Z1scsumpt_cand2==Z1scsumpt_cand1 && Z2scsumpt_cand2>Z2scsumpt_cand1))
        )
      )
      ) theChosenOne = cand2;
    else theChosenOne = cand1;
  }
  else MELAStreamHelpers::MELAerr << "HiggsComparators::candComparator: Scheme " << scheme << " is not defined!" << std::endl;

  PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
  return theChosenOne;
}


