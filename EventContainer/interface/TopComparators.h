#ifndef TOPCOMPARATORS_H
#define TOPCOMPARATORS_H

#include <iostream>
#include "MELAEvent.h"

namespace TopComparators{
  enum CandidateSelection{
    BestTopWMass,
    nCandidateSelections
  };

  MELATopCandidate_t* matchATopToParticle(MELAEvent& ev, MELAParticle* genT);
  MELATopCandidate_t* candidateSelector(MELAEvent& ev, TopComparators::CandidateSelection scheme);
  MELATopCandidate_t* candComparator(MELATopCandidate_t* cand1, MELATopCandidate_t* cand2, TopComparators::CandidateSelection scheme);
}

#endif
