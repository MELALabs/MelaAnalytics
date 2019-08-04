#ifndef HIGGSCOMPARATORS_H
#define HIGGSCOMPARATORS_H

#include <iostream>
#include "MELAEvent.h"

namespace HiggsComparators{
  enum CandidateSelection{
    BestZ1ThenZ2ScSumPt,
    nCandidateSelections
  };

  MELACandidate* matchAHiggsToParticle(MELAEvent& ev, MELAParticle* genH);
  MELACandidate* candidateSelector(MELAEvent& ev, HiggsComparators::CandidateSelection scheme, MELAEvent::CandidateVVMode VVMode);
  MELACandidate* candComparator(MELACandidate* cand1, MELACandidate* cand2, HiggsComparators::CandidateSelection scheme, MELAEvent::CandidateVVMode VVMode);
}

#endif
