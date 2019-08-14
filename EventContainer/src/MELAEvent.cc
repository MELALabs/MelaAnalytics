#include <cassert>
#include <iostream>
#include <algorithm>
#include "MELAEvent.h"
#include "TopComparators.h"
#include "MELAStreamHelpers.hh"

using namespace PDGHelpers;
using namespace ParticleComparators;
using namespace MELAStreamHelpers;


void MELAEvent::applyParticleSelection(){
  applyLeptonSelection();
  applyNeutrinoSelection();
  applyPhotonSelection();
  applyJetSelection();
  applyTopCandidateSelection();
  applyCandidateSelection(); // Order matters here
}
void MELAEvent::applyLeptonSelection(){
  for (std::vector<MELAParticle*>::iterator it = leptons.begin(); it!=leptons.end(); it++){
    // Trigger and acceptance
    bool passAcceptance = true;
    if (std::abs((*it)->id)==11 && ((*it)->pt()<=electronPtCut || std::abs((*it)->eta())>=electronEtaCut)) passAcceptance = false;
    else if (std::abs((*it)->id)==13 && ((*it)->pt()<=muonPtCut || std::abs((*it)->eta())>=muonEtaCut)) passAcceptance = false;
    else if (std::abs((*it)->id)==15) passAcceptance = false;
    for (std::vector<MELAParticle*>::iterator it2 = leptons.begin(); it2!=leptons.end(); it2++){
      if ((*it2)==(*it)) continue; // Every particle is their own ghost.
      else if ((*it)->deltaR((*it2)->p4)<=ghostDeltaRCut){ passAcceptance = false; break; } // Ghost removal
    }
    (*it)->setSelected(passAcceptance);
  }
}
void MELAEvent::applyNeutrinoSelection(){
  for (std::vector<MELAParticle*>::iterator it = neutrinos.begin(); it!=neutrinos.end(); it++) (*it)->setSelected(false);
}
void MELAEvent::applyPhotonSelection(){
  for (std::vector<MELAParticle*>::iterator it = photons.begin(); it!=photons.end(); it++) (*it)->setSelected(true); // For now...
}
void MELAEvent::applyJetSelection(){
  for (std::vector<MELAParticle*>::iterator it = jets.begin(); it!=jets.end(); it++){
    bool passAcceptance = true;
    if ((*it)->pt()<=jetPtCut || std::abs((*it)->eta())>=jetEtaCut) passAcceptance = false;
    for (std::vector<MELAParticle*>::iterator it2 = leptons.begin(); it2!=leptons.end(); it2++){ // Clean from selected leptons
      if ((*it2)->passSelection){ // If it is not selected at all, why would I care?
        if ((*it)->deltaR((*it2)->p4)<=jetDeltaR){ passAcceptance = false; break; }
      }
    }
    (*it)->setSelected(passAcceptance);
  }
}
void MELAEvent::applyTopCandidateSelection(){
  for (std::vector<MELATopCandidate_t*>::iterator it = topcandidates.begin(); it!=topcandidates.end(); it++){
    (*it)->testPreSelectedDaughters();
    //if (!(*it)->passSelection) continue;
  }
}
void MELAEvent::applyCandidateSelection(){
  for (std::vector<MELACandidate*>::iterator it = candidates.begin(); it!=candidates.end(); it++){
    (*it)->testPreSelectedDaughters();
    if (!(*it)->passSelection) continue;

    bool passAcceptance = true;
    if ((*it)->getSortedV(0) && ((*it)->getSortedV(0)->m()<=mV1LowCut || (*it)->getSortedV(0)->m()>=mV1HighCut)){
      passAcceptance = false; (*it)->getSortedV(0)->setSelected(passAcceptance);
    } // Z1 selection
    if ((*it)->getSortedV(1) && ((*it)->getSortedV(1)->m()<=mV2LowCut || (*it)->getSortedV(1)->m()>=mV2HighCut)){
      passAcceptance = false; (*it)->getSortedV(1)->setSelected(passAcceptance);
    } // Z2 selection

    // Extra V selection, no effect on the main candidate
    for (MELAParticle* extraV:(*it)->getAssociatedSortedVs()){
      for (MELAParticle* dauV:extraV->getDaughters()){ if (!dauV->passSelection) extraV->setSelected(false); break; }
      if (!extraV->passSelection) continue;

      if (isAZBoson(extraV->id) && extraV->getDaughter(0)){
        if (
          (isALepton(extraV->getDaughter(0)->id) && (extraV->m()<=mllLowCut || extraV->m()>=mllHighCut))
          ) extraV->setSelected(false);
      }
    }
    TLorentzVector pLOC[2];
    pLOC[0]=(*it)->getAlternativeVMomentum(0);
    pLOC[1]=(*it)->getAlternativeVMomentum(1);
    if (pLOC[0].M()<=mllLowCut || pLOC[1].M()<=mllLowCut) passAcceptance=false;

    (*it)->setSelected(passAcceptance);
  }
}

void MELAEvent::addCandidate(MELACandidate*& myParticle){
  bool isIdentical = false;
  bool isSame = false;
  for (MELACandidate* testCand:getCandidates()){
    if (testCand==myParticle){ isIdentical=isSame=true; break; }

    bool hasSameDaughters = (myParticle->getNDaughters()==testCand->getNDaughters() && myParticle->getDecayMode()==testCand->getDecayMode() && myParticle->id==testCand->id);
    for (MELAParticle* partD:myParticle->getDaughters()){
      if (!hasSameDaughters) break;
      hasSameDaughters &= testCand->hasDaughter(partD);
    }
    if (hasSameDaughters){ isIdentical=true; break; }
  }
  if (!isIdentical) candidates.push_back(myParticle);
  else if (!isSame){ delete myParticle; myParticle=nullptr; }
}

void MELAEvent::addTopCandidate(MELATopCandidate_t*& myParticle){
  bool isIdentical = false;
  bool isSame = false;
  for (MELATopCandidate_t* testCand:getTopCandidates()){
    if (testCand==myParticle){ isIdentical=isSame=true; break; }

    MELAParticle* testPartnerPart = testCand->getPartnerParticle();
    MELAParticle* testWf = testCand->getWFermion();
    MELAParticle* testWfb = testCand->getWAntifermion();

    MELAParticle* candPartnerPart = myParticle->getPartnerParticle();
    MELAParticle* candWf = myParticle->getWFermion();
    MELAParticle* candWfb = myParticle->getWAntifermion();

    bool hasSameId = (myParticle->id == testCand->id);
    bool hasSamePartnerPart = (candPartnerPart == testPartnerPart);
    bool hasSameWf = (candWf == testWf);
    bool hasSameWfb = (candWfb == testWfb);

    if (hasSameId && hasSamePartnerPart && hasSameWf && hasSameWfb){
      if (!candPartnerPart && !candWf && !candWfb){ // Case of undecayed tops
        double eucdot = myParticle->euclidean_dot(testCand);
        double eucdot_self = myParticle->euclidean_dot(myParticle);
        if (fabs(eucdot - eucdot_self)<1e-5*eucdot_self) isIdentical=true;
      }
      else isIdentical=true;
      if (isIdentical) break;
    }
  }
  if (!isIdentical) topcandidates.push_back(myParticle);
  else if (!isSame){ delete myParticle; myParticle=nullptr; }
}

MELAEvent::CandidateVVMode MELAEvent::getCandidateVVModeFromString(std::string const& s){
  std::string str = s;
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  if (str == "undecayed") return MELAEvent::UndecayedMode;
  else if (str == "ww") return MELAEvent::WWMode;
  else if (str == "zz") return MELAEvent::ZZMode;
  else if (str == "hff" || str == "hffb" || str == "ff" || str == "ffb" || str == "yukawa") return MELAEvent::YukawaMode;
  else if (str == "zgam" || str == "zgamma") return MELAEvent::ZGammaMode;
  else if (str == "gamgam" || str == "gammagamma") return MELAEvent::GammaGammaMode;
  else if (str == "zjets" || str == "zjet" || str == "z" || str == "zff" || str == "zffb") return MELAEvent::ZJetsMode;
  else{
    MELAerr << "MELAEvent::getCandidateVVModeFromString: Cannot convert string " << s << " to a MELAEvent::CandidateVVMode enumerator." << std::endl;
    assert(0);
    return MELAEvent::nCandidateVVModes;
  }
}
void MELAEvent::printCandidateDecayModeDescriptions(){
  MELAout << "MELAEvent supports the following decay mode / final state type combinations:\n"
    << "fstype    / Undecayed / ZZ    / WW     / Yukawa    / Zgam   / gamgam   / ZJets \n"
    << "fstype=0:   -         / 4l    / lnulnu / 2l        / 2l     / gam      / 2l    \n"
    << "fstype=1:   -         / 4q    / qq'QQ' / 2q        / 2q     / -        / 2q    \n"
    << "fstype=2:   -         / 2l2q  / lnuqq' / -         / -      / -        / -     \n"
    << "fstype=3:   -         / 2l2nu / -      / -         / -      / -        / -     \n"
    << "fstype=4:   -         / 2q2nu / -      / -         / -      / -        / -     \n"
    << "fstype=5:   -         / 4nu   / -      / -         / 2nu    / -        / 2nu   \n"
    << "fstype=-1:            / Any                                                    \n"
    << "fstype=-2:  -         / 2l2X  / lnuXX'                                         \n"
    << "fstype=-3:  -         / 2nu2X / -                                              \n"
    << "fstype=-4:  -         / 2q2X  / qq'XX'                                         "
    << std::endl;
}

void MELAEvent::constructVVCandidates(CandidateVVMode const& VVMode, int const& fstype){
  /*
  fstype    / Undecayed / ZZ    / WW     / Yukawa    / Zgam   / gamgam   / ZJets
  fstype=0:   -         / 4l    / lnulnu / 2l        / 2l     / gam      / 2l
  fstype=1:   -         / 4q    / qq'QQ' / 2q        / 2q     / -        / 2q
  fstype=2:   -         / 2l2q  / lnuqq' / -         / -      / -        / -
  fstype=3:   -         / 2l2nu / -      / -         / -      / -        / -
  fstype=4:   -         / 2q2nu / -      / -         / -      / -        / -
  fstype=5:   -         / 4nu   / -      / -         / 2nu    / -        / 2nu
  fstype=-1:            / Any
  fstype=-2:  -         / 2l2X  / lnuXX'
  fstype=-3:  -         / 2nu2X / -
  fstype=-4:  -         / 2q2X  / qq'XX'
  */

  if (
    (VVMode==UndecayedMode && fstype!=-1)
    ||
    (VVMode==WWMode && fstype>2)
    ||
    (VVMode==ZZMode && fstype>5)
    ||
    (VVMode==YukawaMode && fstype>1)
    ||
    (VVMode==ZGammaMode && (fstype>1 && fstype!=5))
    ||
    (VVMode==GammaGammaMode && fstype>0)
    ||
    (VVMode==ZJetsMode && (fstype>1 && fstype!=5))
    ||
    VVMode==nCandidateVVModes
    ||
    (fstype<-1 && !(VVMode==WWMode || VVMode==ZZMode))
    ||
    (fstype==-3 && VVMode==WWMode)
    ||
    fstype<-4
    ){
    if (VVMode==WWMode) std::cerr << "No " << "WW" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (VVMode==ZZMode) std::cerr << "No " << "ZZ" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (VVMode==YukawaMode) std::cerr << "No " << "f-fbar" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (VVMode==ZGammaMode) std::cerr << "No " << "Zgamma" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (VVMode==GammaGammaMode) std::cerr << "No " << "gammagamma" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (VVMode==ZJetsMode) std::cerr << "No " << "Z+(n)jets" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (VVMode==UndecayedMode) std::cerr << "No " << "undecayed" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else std::cerr << "Unknown candidate decay mode " << VVMode << " and final state " << fstype << "!" << std::endl;
    MELAEvent::printCandidateDecayModeDescriptions();
    return;
  }

  std::vector<MELAParticle*> lepMinusPlus[3][2]; // l-, l+
  std::vector<MELAParticle*> lepNuNubar[3][2]; // nu, nub
  std::vector<MELAParticle*> quarkAntiquark[7][2]; // q, qb

  for (std::vector<MELAParticle*>::iterator it = leptons.begin(); it!=leptons.end(); it++){ // Leptons
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==11) iFirst = 0;
    else if (abs((*it)->id)==13) iFirst = 1;
    else if (abs((*it)->id)==15) iFirst = 2;
    if ((*it)->id<0) iSecond=1;
    lepMinusPlus[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<MELAParticle*>::iterator it = neutrinos.begin(); it!=neutrinos.end(); it++){ // Neutrinos
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==12) iFirst = 0;
    else if (abs((*it)->id)==14) iFirst = 1;
    else if (abs((*it)->id)==16) iFirst = 2;
    if ((*it)->id<0) iSecond=1;
    lepNuNubar[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<MELAParticle*>::iterator it = jets.begin(); it!=jets.end(); it++){ // Jets
    int iFirst=abs((*it)->id); // Yes, 0-6, 0 being unknown
    if (PDGHelpers::isAGluon(iFirst)) continue;
    int iSecond=0;
    if ((*it)->id<0) iSecond=1;
    quarkAntiquark[iFirst][iSecond].push_back(*it);
  }
  for (auto* top:topcandidates){ // Tops
    int iFirst=abs(top->id); // 0 or 6
    if (!PDGHelpers::isAnUnknownJet(iFirst) && !PDGHelpers::isATopQuark(iFirst)) continue;
    int iSecond=0;
    if (top->id<0) iSecond=1;
    quarkAntiquark[iFirst][iSecond].push_back(top);
  }

  std::vector<MELAParticle*> tmpVhandle;

  if (VVMode==ZZMode || VVMode==ZGammaMode){ // ZZ or Zgam

    if (fstype<0 || (VVMode==ZZMode && (fstype==0 || fstype==2 || fstype==3)) || (VVMode==ZGammaMode && fstype==0)){ // Z->2l
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepMinusPlus[c][0]){
          for (MELAParticle* F2:lepMinusPlus[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pV = F1->p4 + F2->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(F1);
            V->addDaughter(F2);
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || (VVMode==ZZMode && (fstype==3 || fstype==4 || fstype==5)) || (VVMode==ZGammaMode && fstype==5)){ // Z->2nu
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepNuNubar[c][0]){
          for (MELAParticle* F2:lepNuNubar[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pV = F1->p4 + F2->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(F1);
            V->addDaughter(F2);
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || (VVMode==ZZMode && (fstype==1 || fstype==2 || fstype==4)) || (VVMode==ZGammaMode && fstype==1)){ // Z->2q
      for (int c=1; c<7; c++){
        for (MELAParticle* F1:quarkAntiquark[c][0]){
          for (MELAParticle* F2:quarkAntiquark[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pV = F1->p4 + F2->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(F1);
            V->addDaughter(F2);
            tmpVhandle.push_back(V);
          }
        }
      }
    }

  }
  else if (VVMode==WWMode){ // WW

    if (fstype<0 || fstype==0 || fstype==2){ // W->lnu
      for (int c=0; c<3; c++){
        for (int signW=0; signW<2; signW++){ // ==0: W+, ==1: W-
          for (MELAParticle* F1:lepMinusPlus[c][1-signW]){
            for (MELAParticle* F2:lepNuNubar[c][signW]){
              if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
              TLorentzVector pV = F1->p4 + F2->p4;
              MELAParticle* V = new MELAParticle(24*(1-2*signW), pV);
              V->addDaughter(F1);
              V->addDaughter(F2);
              tmpVhandle.push_back(V);
            }
          }
        }
      }
    }
    if (fstype<0 || fstype==1 || fstype==2){ // W->2q
      for (int c=1; c<7; c++){
        for (int d=1; d<7; d++){
          if (d==c) continue;
          for (MELAParticle* F1:quarkAntiquark[c][0]){
            for (MELAParticle* F2:quarkAntiquark[d][1]){
              if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
              int totalcharge = F1->charge() + F2->charge();
              if (abs(totalcharge)!=1) continue;

              TLorentzVector pV = F1->p4 + F2->p4;
              MELAParticle* V = new MELAParticle(24*totalcharge, pV);
              V->addDaughter(F1);
              V->addDaughter(F2);
              tmpVhandle.push_back(V);
            }
          }
        }
      }
    }
  }
  else if (VVMode==YukawaMode){ // H->f fbar

    if (fstype<0 || fstype==0){ // H->2l
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepMinusPlus[c][0]){
          for (MELAParticle* F2:lepMinusPlus[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pH = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(25, pH, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
            setCandidateDecayMode(TVar::CandidateDecay_ff);
            cand->sortDaughters();
            setCandidateDecayMode(defaultHDecayMode);
            addCandidate(cand);
          }
        }
      }
    }
    if (fstype<0 || fstype==1){ // H->2q
      for (int c=1; c<7; c++){
        for (MELAParticle* F1:quarkAntiquark[c][0]){
          for (MELAParticle* F2:quarkAntiquark[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pH = F1->p4+F2->p4;
            MELACandidate* cand = new MELACandidate(25, pH, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
            setCandidateDecayMode(TVar::CandidateDecay_ff);
            cand->sortDaughters();
            setCandidateDecayMode(defaultHDecayMode);
            addCandidate(cand);
          }
        }
      }
    }

  }
  else if (VVMode==ZJetsMode){ // Z->f fbar

    if (fstype<0 || fstype==0){ // Z->2l
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepMinusPlus[c][0]){
          for (MELAParticle* F2:lepMinusPlus[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pCand = F1->p4 + F2->p4;
            MELACandidate* cand = new MELACandidate(23, pCand, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
            setCandidateDecayMode(TVar::CandidateDecay_ff);
            cand->sortDaughters();
            setCandidateDecayMode(defaultHDecayMode);
            addCandidate(cand);
          }
        }
      }
    }
    if (fstype<0 || fstype==1){ // Z->2q
      for (int c=1; c<7; c++){
        for (MELAParticle* F1:quarkAntiquark[c][0]){
          for (MELAParticle* F2:quarkAntiquark[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pCand = F1->p4 + F2->p4;
            MELACandidate* cand = new MELACandidate(23, pCand, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
            setCandidateDecayMode(TVar::CandidateDecay_ff);
            cand->sortDaughters();
            setCandidateDecayMode(defaultHDecayMode);
            addCandidate(cand);
          }
        }
      }
    }
    if (fstype<0 || fstype==5){ // Z->2nu
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepNuNubar[c][0]){
          for (MELAParticle* F2:lepNuNubar[c][1]){
            if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
            TLorentzVector pCand = F1->p4 + F2->p4;
            MELACandidate* cand = new MELACandidate(23, pCand, true);
            cand->addDaughter(F1);
            cand->addDaughter(F2);

            TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
            setCandidateDecayMode(TVar::CandidateDecay_ff);
            cand->sortDaughters();
            setCandidateDecayMode(defaultHDecayMode);
            addCandidate(cand);
          }
        }
      }
    }

  }
  else if (VVMode==UndecayedMode){ // Undecayed
    for (std::vector<MELAParticle*>::iterator it = intermediates.begin(); it!=intermediates.end(); it++){ // Add directly
      if (isAHiggs((*it)->id)){
        MELACandidate* cand = new MELACandidate(25, (*it)->p4, true);

        TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
        setCandidateDecayMode(TVar::CandidateDecay_Stable);
        cand->sortDaughters();
        setCandidateDecayMode(defaultHDecayMode);
        addCandidate(cand);
      }
    }
  }

  if (debugVars::debugFlag) std::cout << "Number of V/ZZ before sorting photons: " << tmpVhandle.size() << " " << getNCandidates() << std::endl;

  if (VVMode==ZGammaMode || VVMode==GammaGammaMode){
    for (MELAParticle* part:photons){ // Copy the photons
      MELAParticle* V = new MELAParticle(part->id, part->p4);
      V->addDaughter(part); // Photon is its own daughter!
      tmpVhandle.push_back(V);
    }
  }

  if (debugVars::debugFlag) std::cout << "Number of V/ZZ after sorting photons: " << tmpVhandle.size() << " " << getNCandidates() << std::endl;

  if (
    ((fstype<0 || fstype==1 || fstype==2 || fstype==4) && (VVMode==WWMode || VVMode==ZZMode)) // W/Z->2j reco.-level
    ||
    ((fstype<0 || fstype==1) && VVMode==YukawaMode) // H->2j reco.-level
    ||
    ((fstype<0 || fstype==1) && VVMode==ZGammaMode) // H->Zgam with Z->2j
    ||
    ((fstype<0 || fstype==1) && VVMode==ZJetsMode) // Z+2jets with Z->2j
    ){
    for (std::vector<MELAParticle*>::iterator it1 = quarkAntiquark[0][0].begin(); it1!=quarkAntiquark[0][0].end(); it1++){
      MELAParticle* F1 = *it1;
      if (F1->id!=0) continue;
      for (std::vector<MELAParticle*>::iterator it2=it1; it2!=quarkAntiquark[0][0].end(); it2++){
        if (it1==it2) continue;
        MELAParticle* F2 = *it2;
        if (F2->id!=0) continue;

        if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;

        if (VVMode==WWMode || VVMode==ZZMode || VVMode==ZGammaMode){
          TLorentzVector pV = F1->p4 + F2->p4;
          MELAParticle* V = new MELAParticle(0, pV);
          V->addDaughter(F1);
          V->addDaughter(F2);
          tmpVhandle.push_back(V);
        }
        else if (VVMode==YukawaMode || VVMode==ZJetsMode){
          TLorentzVector pCand = F1->p4 + F2->p4;
          int candId = (VVMode==ZJetsMode ? 23 : 25);
          MELACandidate* cand = new MELACandidate(candId, pCand, true);
          cand->addDaughter(F1);
          cand->addDaughter(F2);

          TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
          setCandidateDecayMode(TVar::CandidateDecay_ff);
          cand->sortDaughters();
          setCandidateDecayMode(defaultHDecayMode);
          addCandidate(cand);
        }
      }
    }
    if (debugVars::debugFlag) std::cout << "Number of V/ZZ after sorting reco. jets: " << tmpVhandle.size() << " " << getNCandidates() << std::endl;
  }



  for (std::vector<MELAParticle*>::iterator it1 = tmpVhandle.begin(); it1!=tmpVhandle.end(); it1++){
    MELAParticle* Vi = *it1;
    for (std::vector<MELAParticle*>::iterator it2 = it1; it2!=tmpVhandle.end(); it2++){
      if (it2==it1) continue;
      MELAParticle* Vj = *it2;
      if (Vj==Vi) continue;
      if ((Vi->charge() + Vj->charge()) != 0) continue;

      MELAParticle* Vi1 = Vi->getDaughter(0);
      MELAParticle* Vi2 = Vi->getDaughter(1);
      MELAParticle* Vj1 = Vj->getDaughter(0);
      MELAParticle* Vj2 = Vj->getDaughter(1);
      /*
      std::cout << "11: " << Vi1->id << '\t' << Vi1->x() << '\t' << Vi1->y() << '\t' << Vi1->z() << '\t' << Vi1->t() << '\t' << std::endl;
      std::cout << "12: " << Vi2->id << '\t' << Vi2->x() << '\t' << Vi2->y() << '\t' << Vi2->z() << '\t' << Vi2->t() << '\t' << std::endl;
      std::cout << "21: " << Vj1->id << '\t' << Vj1->x() << '\t' << Vj1->y() << '\t' << Vj1->z() << '\t' << Vj1->t() << '\t' << std::endl;
      std::cout << "22: " << Vj2->id << '\t' << Vj2->x() << '\t' << Vj2->y() << '\t' << Vj2->z() << '\t' << Vj2->t() << '\t' << std::endl;
      */
      if (Vi1==Vj1 || (Vi2 && Vi2==Vj2)) continue;
      bool createCandidate=true;
      if ((VVMode==WWMode || VVMode==ZZMode) && fstype<-1){
        unsigned int partcounter=0;
        if (fstype==-2){ // Count leptons
          if (Vi1 && PDGHelpers::isALepton(Vi1->id)) partcounter++;
          if (Vi2 && PDGHelpers::isALepton(Vi2->id)) partcounter++;
          if (Vj1 && PDGHelpers::isALepton(Vj1->id)) partcounter++;
          if (Vj2 && PDGHelpers::isALepton(Vj2->id)) partcounter++;
        }
        if ((fstype==-3 && VVMode==ZZMode) || (fstype==-2 && VVMode==WWMode)){ // Count neutrinos
          if (Vi1 && PDGHelpers::isANeutrino(Vi1->id)) partcounter++;
          if (Vi2 && PDGHelpers::isANeutrino(Vi2->id)) partcounter++;
          if (Vj1 && PDGHelpers::isANeutrino(Vj1->id)) partcounter++;
          if (Vj2 && PDGHelpers::isANeutrino(Vj2->id)) partcounter++;
        }
        if (fstype==-4){ // Count jets
          if (Vi1 && PDGHelpers::isAJet(Vi1->id)) partcounter++;
          if (Vi2 && PDGHelpers::isAJet(Vi2->id)) partcounter++;
          if (Vj1 && PDGHelpers::isAJet(Vj1->id)) partcounter++;
          if (Vj2 && PDGHelpers::isAJet(Vj2->id)) partcounter++;
        }
        if (partcounter<2) createCandidate=false;
      }
      if (!createCandidate) continue;

      if (debugVars::debugFlag){
        if (Vi1) std::cout << "Vi1 not zero. Id: " << Vi1->id << std::endl;
        if (Vi2) std::cout << "Vi2 not zero. Id: " << Vi2->id << std::endl;
        if (Vj1) std::cout << "Vj1 not zero. Id: " << Vj1->id << std::endl;
        if (Vj2) std::cout << "Vj2 not zero. Id: " << Vj2->id << std::endl;
      }

      TLorentzVector pH(0, 0, 0, 0);
      if (Vi1) pH = pH + Vi1->p4;
      if (Vi2) pH = pH + Vi2->p4;
      if (Vj1) pH = pH + Vj1->p4;
      if (Vj2) pH = pH + Vj2->p4;
      MELACandidate* cand = new MELACandidate(25, pH, true);

      if (Vi1) cand->addDaughter(Vi1);
      if (Vi2) cand->addDaughter(Vi2);
      if (Vj1) cand->addDaughter(Vj1);
      if (Vj2) cand->addDaughter(Vj2);

      TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
      if (VVMode==WWMode) setCandidateDecayMode(TVar::CandidateDecay_WW);
      else if (VVMode==ZZMode) setCandidateDecayMode(TVar::CandidateDecay_ZZ);
      else if (VVMode==ZGammaMode) setCandidateDecayMode(TVar::CandidateDecay_ZG);
      else if (VVMode==GammaGammaMode) setCandidateDecayMode(TVar::CandidateDecay_GG);
      else setCandidateDecayMode(TVar::CandidateDecay_ff);
      if (debugVars::debugFlag) std::cout << "Sorting daughters..." << std::endl;
      cand->sortDaughters();
      if (debugVars::debugFlag) std::cout << "Sorted daughters successfully!" << std::endl;
      setCandidateDecayMode(defaultHDecayMode);
      addCandidate(cand);
    }
  }

  for (MELAParticle*& tmpV:tmpVhandle) delete tmpV;
  if (debugVars::debugFlag) std::cout << "tmpVhandle deletion step is done." << std::endl;
  tmpVhandle.clear();
}

void MELAEvent::constructTopCandidates(){
  std::vector<MELAParticle*> lepMinusPlus[3][2]; // l-, l+
  std::vector<MELAParticle*> lepNuNubar[3][2]; // nu, nub
  std::vector<MELAParticle*> quarkAntiquark[6][2]; // q, qb

  for (std::vector<MELAParticle*>::iterator it = leptons.begin(); it!=leptons.end(); it++){ // Leptons
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==11) iFirst = 0;
    else if (abs((*it)->id)==13) iFirst = 1;
    else if (abs((*it)->id)==15) iFirst = 2;
    if ((*it)->id<0) iSecond=1;
    lepMinusPlus[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<MELAParticle*>::iterator it = neutrinos.begin(); it!=neutrinos.end(); it++){ // Neutrinos
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==12) iFirst = 0;
    else if (abs((*it)->id)==14) iFirst = 1;
    else if (abs((*it)->id)==16) iFirst = 2;
    if ((*it)->id<0) iSecond=1;
    lepNuNubar[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<MELAParticle*>::iterator it = jets.begin(); it!=jets.end(); it++){ // Jets
    int iFirst=abs((*it)->id); // Yes, 0-6, 0 being unknown
    if (PDGHelpers::isAGluon(iFirst)) continue;
    if (PDGHelpers::isATopQuark(iFirst)) continue;
    int iSecond=0;
    if ((*it)->id<0) iSecond=1;
    quarkAntiquark[iFirst][iSecond].push_back(*it);
  }

  std::vector<MELAParticle*> tmpVhandle;

  for (int c=0; c<3; c++){
    for (int signW=0; signW<2; signW++){ // ==0: W+, ==1: W-
      for (MELAParticle* F1:lepMinusPlus[c][1-signW]){
        for (MELAParticle* F2:lepNuNubar[c][signW]){
          if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
          TLorentzVector pV = F1->p4 + F2->p4;
          MELAParticle* V = new MELAParticle(24*(1-2*signW), pV);
          V->addDaughter(F1);
          V->addDaughter(F2);
          tmpVhandle.push_back(V);
        }
      }
    }
  }
  for (int c=1; c<6; c++){
    for (int d=1; d<6; d++){
      if (d==c) continue;
      for (MELAParticle* F1:quarkAntiquark[c][0]){
        for (MELAParticle* F2:quarkAntiquark[d][1]){
          if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;
          int totalcharge = F1->charge() + F2->charge();
          if (abs(totalcharge)!=1) continue;

          TLorentzVector pV = F1->p4 + F2->p4;
          MELAParticle* V = new MELAParticle(24*totalcharge, pV);
          V->addDaughter(F1);
          V->addDaughter(F2);
          tmpVhandle.push_back(V);
        }
      }
    }
  }

  // Reco. jet Vs
  for (std::vector<MELAParticle*>::iterator it1 = quarkAntiquark[0][0].begin(); it1!=quarkAntiquark[0][0].end(); it1++){
    MELAParticle* F1 = *it1;
    if (F1->id!=0) continue;
    for (std::vector<MELAParticle*>::iterator it2=it1; it2!=quarkAntiquark[0][0].end(); it2++){
      if (it1==it2) continue;
      MELAParticle* F2 = *it2;
      if (F2->id!=0) continue;

      if (MELAParticle::checkDeepDaughtership(F1, F2)) continue;

      TLorentzVector pV = F1->p4 + F2->p4;
      MELAParticle* V = new MELAParticle(0, pV);
      V->addDaughter(F1);
      V->addDaughter(F2);
      tmpVhandle.push_back(V);
    }
  }


  std::vector<MELATopCandidate_t*> alltops;
  std::vector<MELATopCandidate_t*> topsToRemove;
  // Match quarks to Vs
  for (int c=1; c<6; c++){
    for (int qaq=0; qaq<2; qaq++){
      for (MELAParticle* quark:quarkAntiquark[c][qaq]){
        for (MELAParticle* partV:tmpVhandle){
          //if (partV->id==0) continue;
          if (partV->charge()*quark->charge()>0.) continue;
          if (MELAParticle::checkDeepDaughtership(partV, quark)) continue;

          MELATopCandidate_t* cand = new MELATopCandidate_t(6*(1-2*qaq), (partV->p4+quark->p4));
          cand->setPartnerParticle(quark);
          cand->setWFermion(partV->getDaughter(0));
          cand->setWAntifermion(partV->getDaughter(1));

          alltops.push_back(cand);
        }
      }
    }
  }

  // Match reco. jets to reco. Vs
  for (MELAParticle* quark:quarkAntiquark[0][0]){
    for (MELAParticle* partV:tmpVhandle){
      if (MELAParticle::checkDeepDaughtership(partV, quark)) continue;
      //if (PDGHelpers::isAQuark(partV->getDaughter(0)) || PDGHelpers::isAQuark(partV->getDaughter(1))) continue;

      MELATopCandidate_t* cand = new MELATopCandidate_t(0, (partV->p4+quark->p4));
      cand->setPartnerParticle(quark);
      cand->setWFermion(partV->getDaughter(0));
      cand->setWAntifermion(partV->getDaughter(1));

      alltops.push_back(cand);
    }
  }

  // Disambiguate the tops that have the same daughters and pick the better ones
  for (MELATopCandidate_t* top:alltops){
    bool hasBetterCand=false;
    for (MELATopCandidate_t* tmptop:alltops){
      if (tmptop==top) continue;
      if (tmptop->id!=top->id) continue;
      if (top->hasDaughter(tmptop->getPartnerParticle()) && top->hasDaughter(tmptop->getWFermion()) && top->hasDaughter(tmptop->getWAntifermion())){
        if (top!=TopComparators::candComparator(top, tmptop, TopComparators::BestTopWMass)){ hasBetterCand=true; break; }
      }
    }
    if (!hasBetterCand) addTopCandidate(top);
    else topsToRemove.push_back(top);
  }

  // Remove unnecessary tops
  for (MELATopCandidate_t*& tmptop:topsToRemove) delete tmptop;
  topsToRemove.clear();

  for (MELAParticle*& tmpV:tmpVhandle) delete tmpV;
  if (debugVars::debugFlag) std::cout << "tmpVhandle deletion step is done." << std::endl;
  tmpVhandle.clear();

  for (std::vector<MELAParticle*>::iterator it = intermediates.begin(); it!=intermediates.end(); it++){ // Add intermediate, undecayed tops
    if (isATopQuark((*it)->id)){
      MELATopCandidate_t* cand = new MELATopCandidate_t((*it)->id, (*it)->p4);
      addTopCandidate(cand);
    }
  }
}

MELACandidate* MELAEvent::getCandidate(size_t const& index)const{
  if (candidates.size()>index) return candidates.at(index);
  else return nullptr;
}
MELATopCandidate_t* MELAEvent::getTopCandidate(size_t const& index)const{
  if ( topcandidates.size()>index) return topcandidates.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getLepton(size_t const& index)const{
  if (leptons.size()>index) return leptons.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getNeutrino(size_t const& index)const{
  if (neutrinos.size()>index) return neutrinos.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getPhoton(size_t const& index)const{
  if (photons.size()>index) return photons.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getJet(size_t const& index)const{
  if (jets.size()>index) return jets.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getMother(size_t const& index)const{
  if ( mothers.size()>index) return mothers.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getIntermediate(size_t const& index)const{
  if (intermediates.size()>index) return intermediates.at(index);
  else return nullptr;
}
MELAParticle* MELAEvent::getParticle(size_t const& index)const{
  if (particles.size()>index) return particles.at(index);
  else return nullptr;
}

TLorentzVector MELAEvent::missingP() const{
  TLorentzVector totalP(0, 0, 0, 0);
  for (unsigned int pp=0; pp<particles.size(); pp++){
    MELAParticle* part = getParticle(pp);
    if (part->passSelection) totalP = totalP + part->p4;
  }
  totalP.SetT(totalP.P());
  totalP.SetVect(-totalP.Vect());
  return totalP;
}

void MELAEvent::addVVCandidateAppendages(){
  for (std::vector<MELACandidate*>::iterator it = candidates.begin(); it!=candidates.end(); it++){
    for (std::vector<MELAParticle*>::iterator iL = leptons.begin(); iL!=leptons.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedLepton(*iL); }
    for (std::vector<MELAParticle*>::iterator iN = neutrinos.begin(); iN!=neutrinos.end(); iN++){ if ((*iN)->passSelection) (*it)->addAssociatedNeutrino(*iN); }
    for (std::vector<MELAParticle*>::iterator iP = photons.begin(); iP!=photons.end(); iP++){ if ((*iP)->passSelection) (*it)->addAssociatedPhoton(*iP); }
    for (std::vector<MELAParticle*>::iterator iJ = jets.begin(); iJ!=jets.end(); iJ++){ if ((*iJ)->passSelection) (*it)->addAssociatedJet(*iJ); }
    (*it)->addAssociatedVs();
    for (std::vector<MELATopCandidate_t*>::iterator iT = topcandidates.begin(); iT!=topcandidates.end(); iT++){ if ((*iT)->passSelection) (*it)->addAssociatedTop(*iT); }
    for (std::vector<MELAParticle*>::iterator iM = mothers.begin(); iM!=mothers.end(); iM++) (*it)->addMother(*iM);
  }
}

void MELAEvent::wipeAll(){
  leptons.clear();
  neutrinos.clear();
  photons.clear();
  jets.clear();
  wipeArray(candidates, true);
  wipeArray(topcandidates, true);
  wipeArray(intermediates, false);
  wipeArray(particles, false);
};
