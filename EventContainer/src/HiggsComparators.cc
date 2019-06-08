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

MELACandidate* HiggsComparators::candidateSelector(MELAEvent& ev, HiggsComparators::CandidateSelection scheme, int isZZ){
  MELACandidate* cand=nullptr;
  for (MELACandidate* tmpCand:ev.getCandidates()){
    if (!tmpCand) continue;
    if (!tmpCand->passSelection) continue;
    if (!cand) cand=tmpCand;
    else cand = HiggsComparators::candComparator(cand, tmpCand, scheme, isZZ);
  }
  return cand;
}

MELACandidate* HiggsComparators::candComparator(MELACandidate* cand1, MELACandidate* cand2, HiggsComparators::CandidateSelection scheme, int isZZ){
  MELACandidate* theChosenOne=nullptr;

  if (!cand1 && !cand2) return theChosenOne;
  else if (cand1 && !cand2) return cand1;
  else if (!cand1 && cand2) return cand2;

  TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
  if (isZZ==0) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_WW);
  else if (isZZ==1) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  else if (isZZ==3) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZG);
  else if (isZZ==4) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_GG);
  else PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ff);

  double HVVmass = PDGHelpers::Zeromass;
  if (isZZ==0){
    HVVmass = PDGHelpers::Wmass;
  }
  else if (isZZ==1){
    HVVmass = PDGHelpers::Zmass;
  }

  if (isZZ==-1){
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
    double Z2scsumpt_cand1=0, Z2scsumpt_cand2=0;
    MELAParticle* c11 = cand1->getSortedV(1)->getDaughter(0);
    MELAParticle* c12 = cand1->getSortedV(1)->getDaughter(1);
    MELAParticle* c21 = cand2->getSortedV(1)->getDaughter(0);
    MELAParticle* c22 = cand2->getSortedV(1)->getDaughter(1);
    if (c11!=0) Z2scsumpt_cand1 += c11->pt();
    if (c12!=0) Z2scsumpt_cand1 += c12->pt();
    if (c21!=0) Z2scsumpt_cand2 += c21->pt();
    if (c22!=0) Z2scsumpt_cand2 += c22->pt();
    if (
      (diffmass1>diffmass2)
      ||
      (diffmass1==diffmass2 && Z2scsumpt_cand2>Z2scsumpt_cand1)
      ) theChosenOne = cand2;
    else theChosenOne = cand1;
  }
  else MELAStreamHelpers::MELAerr << "HiggsComparators::candComparator: Scheme " << scheme << " is not defined!" << std::endl;

  PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
  return theChosenOne;
}


