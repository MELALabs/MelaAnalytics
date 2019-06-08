#include "TopComparators.h"
#include "MELAStreamHelpers.hh"


MELATopCandidate_t* TopComparators::matchATopToParticle(MELAEvent& ev, MELAParticle* genT){
  MELATopCandidate_t* cand=nullptr;
  if (!genT) return cand;
  for (MELATopCandidate_t* tmpCand:ev.getTopCandidates()){
    if (!tmpCand) continue;
    double genmassquant = genT->m()+genT->pt()+fabs(genT->z());
    double massdiff = fabs(genmassquant-tmpCand->m()-tmpCand->pt()-fabs(tmpCand->z()));
    double massratio = 0;
    if (genmassquant>0.) massratio = massdiff / genmassquant;
    if (massratio<0.001){
      if (!cand) cand = tmpCand;
      else{
        TLorentzVector const& vGen = genT->p4;
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

MELATopCandidate_t* TopComparators::candidateSelector(MELAEvent& ev, TopComparators::CandidateSelection scheme){
  MELATopCandidate_t* cand=nullptr;
  for (MELATopCandidate_t* tmpCand:ev.getTopCandidates()){
    if (!tmpCand) continue;
    if (!tmpCand->passSelection) continue;
    if (!cand) cand=tmpCand;
    else cand = TopComparators::candComparator(cand, tmpCand, scheme);
  }
  return cand;
}

MELATopCandidate_t* TopComparators::candComparator(MELATopCandidate_t* cand1, MELATopCandidate_t* cand2, TopComparators::CandidateSelection scheme){
  MELATopCandidate_t* theChosenOne=nullptr;

  if (!cand1 && !cand2) return theChosenOne;
  else if (cand1 && !cand2) return cand1;
  else if (!cand1 && cand2) return cand2;

  constexpr double GeVunit = 0.01;
  double const mT = PDGHelpers::Topmass*GeVunit;
  double const mW = PDGHelpers::Wmass*GeVunit;
  double const GaT = PDGHelpers::Topwidth*GeVunit;
  double const GaW = PDGHelpers::Wwidth*GeVunit;
  double const mTsq = mT*mT;
  double const mWsq = mW*mW;
  double const mT_GaTsq = pow(mT*GaT, 2);
  double const mW_GaWsq = pow(mW*GaW, 2);

  if (scheme==TopComparators::BestTopWMass){
    double BWinv1 = pow(pow(cand1->m()*GeVunit, 2)-mTsq, 2)+mT_GaTsq;
    double BWinv2 = pow(pow(cand2->m()*GeVunit, 2)-mTsq, 2)+mT_GaTsq;
    double const mW1 = cand1->getWmass();
    double const mW2 = cand2->getWmass();
    if (mW1>=0. && mW2>=0.){
      BWinv1 *= pow(pow(mW1*GeVunit, 2)-mWsq, 2)+mW_GaWsq;
      BWinv2 *= pow(pow(mW2*GeVunit, 2)-mWsq, 2)+mW_GaWsq;
    }
    theChosenOne = (BWinv1<=BWinv2 ? cand1 : cand2);
  }
  else MELAStreamHelpers::MELAerr << "TopComparators::candComparator: Scheme " << scheme << " is not defined!" << std::endl;

  return theChosenOne;
}


