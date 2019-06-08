#include <iostream>
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
    if (std::abs((*it)->id)==11 && ((*it)->pt()<=electronPTCut || std::abs((*it)->eta())>=electronEtaAcceptanceCut)) passAcceptance = false;
    else if (std::abs((*it)->id)==13 && ((*it)->pt()<=muonPTCut || std::abs((*it)->eta())>=muonEtaAcceptanceCut)) passAcceptance = false;
    else if (std::abs((*it)->id)==15) passAcceptance = false;
    for (std::vector<MELAParticle*>::iterator it2 = leptons.begin(); it2<leptons.end(); it2++){
      if ((*it2)==(*it)) continue; // Every particle is their own ghost.
      else if ((*it)->deltaR((*it2)->p4)<=ghostDeltaRCut) passAcceptance = false; // Ghost removal
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
    if ((*it)->pt()<=jetPTCut || std::abs((*it)->eta())>=jetEtaAcceptanceCut) passAcceptance = false; // ZZ4l selection and acceptance
    for (std::vector<MELAParticle*>::iterator it2 = leptons.begin(); it2<leptons.end(); it2++){ // Clean from selected leptons
      if ((*it2)->passSelection){ // If it is not selected at all, why would I care?
        if ((*it)->deltaR((*it2)->p4)<=jetDeltaR) passAcceptance = false;
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
    if ((*it)->getSortedV(0)->m()<=mV1LowCut || (*it)->getSortedV(0)->m()>=mV12HighCut){
      passAcceptance = false; (*it)->getSortedV(0)->setSelected(passAcceptance);
    } // Z1 selection
    if ((*it)->getSortedV(1)->m()<=mV2LowCut || (*it)->getSortedV(1)->m()>=mV12HighCut){
      passAcceptance = false; (*it)->getSortedV(1)->setSelected(passAcceptance);
    } // Z2 selection
    for (MELAParticle* extraV:(*it)->getSortedVs()){
      if (!isAZBoson(extraV->id)) continue;
      else{
        if (extraV->m()<=mllLowCut || extraV->m()>=mV12HighCut || (extraV->getDaughter(0)!=0 && isANeutrino(extraV->getDaughter(0)->id))) extraV->setSelected(false); // Extra Z selection, no effect on ZZ candidate
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

void MELAEvent::constructVVCandidates(int isZZ, int fstype){
  /*
  fstype  / ZZ==1 / WW==0  / Yukawa==2 / Zgam=3 / gamgam=4 / Z+nj=5
  fstype=0: 4l    / lnulnu / 2l        / 2l     / gam      / 2l
  fstype=1: 4q    / 4q     / 2q        / 2q     / -        / 2q
  fstype=2: 2l2q  / lnu2q  / -         / -      / -        / -
  fstype=3: 2l2nu / -      / -         / -      / -        / -
  fstype=4: 2q2nu / -      / -         / -      / -        / -
  fstype=5: 4nu   / -      / -         / 2nu    / -        / 2nu
  fstype=-1: Any
  fstype=-2: 2l2X
  fstype=-3: 2nu2X
  fstype=-4: 2q2X
  */

  if (
    (isZZ<=0 && fstype>2)
    ||
    (isZZ==1 && fstype>5)
    ||
    (isZZ==2 && fstype>1)
    ||
    (isZZ==3 && (fstype>1 && fstype!=5))
    ||
    (isZZ==4 && fstype>0)
    ||
    (isZZ==5 && (fstype>1 && fstype!=5))
    ||
    isZZ>5
    ||
    (fstype<-1 && isZZ>1)
    ||
    fstype<-4
    ){
    if (isZZ<0) std::cerr << "No " << "undecayed" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==0) std::cerr << "No " << "WW" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==1) std::cerr << "No " << "ZZ" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==2) std::cerr << "No " << "f-fbar" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==3) std::cerr << "No " << "Zgamma" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==4) std::cerr << "No " << "gammagamma" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else if (isZZ==5) std::cerr << "No " << "Z+(n)jets" << " candidate with final state " << fstype << " is possible!" << std::endl;
    else std::cerr << "Unknown candidate decay mode " << isZZ << " and final state " << fstype << "!" << std::endl;
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

  std::vector<MELAParticle*> tmpVhandle;

  if (isZZ==1 || isZZ==3){ // ZZ

    if (fstype<0 || (isZZ==1 && (fstype==0 || fstype==2 || fstype==3)) || (isZZ==3 && fstype==0)){ // Z->2l
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepMinusPlus[c][0]){
          for (MELAParticle* F2:lepMinusPlus[c][1]){
            TLorentzVector pV = F1->p4 + F2->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(F1);
            V->addDaughter(F2);
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || (isZZ==1 && (fstype==3 || fstype==4 || fstype==5)) || (isZZ==3 && fstype==5)){ // Z->2nu
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepNuNubar[c][0]){
          for (MELAParticle* F2:lepNuNubar[c][1]){
            TLorentzVector pV = F1->p4 + F2->p4;
            MELAParticle* V = new MELAParticle(23, pV);
            V->addDaughter(F1);
            V->addDaughter(F2);
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype<0 || (isZZ==1 && (fstype==1 || fstype==2 || fstype==4)) || (isZZ==3 && fstype==1)){ // Z->2q
      for (int c=1; c<7; c++){
        for (MELAParticle* F1:quarkAntiquark[c][0]){
          for (MELAParticle* F2:quarkAntiquark[c][1]){
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
  else if(isZZ==0){ // WW

    if (fstype<0 || fstype==0 || fstype==2){ // W->lnu
      for (int c=0; c<3; c++){
        for (int signW=0; signW<2; signW++){ // ==0: W+, ==1: W-
          for (MELAParticle* F1:lepMinusPlus[c][1-signW]){
            for (MELAParticle* F2:lepNuNubar[c][signW]){
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
  else if (isZZ==2){ // H->f fbar

    if (fstype<0 || fstype==0){ // H->2l
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepMinusPlus[c][0]){
          for (MELAParticle* F2:lepMinusPlus[c][1]){
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
  else if (isZZ==5){ // Z->f fbar

    if (fstype<0 || fstype==0){ // Z->2l
      for (int c=0; c<3; c++){
        for (MELAParticle* F1:lepMinusPlus[c][0]){
          for (MELAParticle* F2:lepMinusPlus[c][1]){
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
  else{ // Undecayed
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

  if (isZZ==3 || isZZ==4){
    for (MELAParticle* part:photons){ // Copy the photons
      MELAParticle* V = new MELAParticle(part->id, part->p4);
      V->addDaughter(part); // Photon is its own daughter!
      tmpVhandle.push_back(V);
    }
  }

  if (debugVars::debugFlag) std::cout << "Number of V/ZZ after sorting photons: " << tmpVhandle.size() << " " << getNCandidates() << std::endl;

  if (
    ((fstype<0 || fstype==1 || fstype==2 || fstype==4) && (isZZ==0 || isZZ==1)) // W/Z->2j reco.-level
    ||
    ((fstype<0 || fstype==1) && isZZ==2) // H->2j reco.-level
    ||
    ((fstype<0 || fstype==1) && isZZ==3) // H->Zgam with Z->2j
    ||
    ((fstype<0 || fstype==1) && isZZ==5) // Z+2jets with Z->2j
    ){
    for (std::vector<MELAParticle*>::iterator it1 = quarkAntiquark[0][0].begin(); it1!=quarkAntiquark[0][0].end(); it1++){
      MELAParticle* F1 = *it1;
      if (F1->id!=0) continue;
      for (std::vector<MELAParticle*>::iterator it2=it1; it2!=quarkAntiquark[0][0].end(); it2++){
        if (it1==it2) continue;
        MELAParticle* F2 = *it2;
        if (F2->id!=0) continue;

        if (isZZ==0 || isZZ==1 || isZZ==3){
          TLorentzVector pV = F1->p4 + F2->p4;
          MELAParticle* V = new MELAParticle(0, pV);
          V->addDaughter(F1);
          V->addDaughter(F2);
          tmpVhandle.push_back(V);
        }
        else if (isZZ==2){
          TLorentzVector pCand = F1->p4 + F2->p4;
          MELACandidate* cand = new MELACandidate(25, pCand, true);
          cand->addDaughter(F1);
          cand->addDaughter(F2);

          TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
          setCandidateDecayMode(TVar::CandidateDecay_ff);
          cand->sortDaughters();
          setCandidateDecayMode(defaultHDecayMode);
          addCandidate(cand);
        }
        else if (isZZ==5){
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
      if (Vi1==Vj1 || (Vi2==Vj2 && Vi2 != 0)) continue;
      bool createCandidate=true;
      if (isZZ<=1 && fstype<-1){
        unsigned int partcounter=0;
        if (fstype==-2){ // Count leptons
          if (Vi1!=0 && PDGHelpers::isALepton(Vi1->id)) partcounter++;
          if (Vi2!=0 && PDGHelpers::isALepton(Vi2->id)) partcounter++;
          if (Vj1!=0 && PDGHelpers::isALepton(Vj1->id)) partcounter++;
          if (Vj2!=0 && PDGHelpers::isALepton(Vj2->id)) partcounter++;
        }
        else if (fstype==-3){ // Count neutrinos
          if (Vi1!=0 && PDGHelpers::isANeutrino(Vi1->id)) partcounter++;
          if (Vi2!=0 && PDGHelpers::isANeutrino(Vi2->id)) partcounter++;
          if (Vj1!=0 && PDGHelpers::isANeutrino(Vj1->id)) partcounter++;
          if (Vj2!=0 && PDGHelpers::isANeutrino(Vj2->id)) partcounter++;
        }
        else if (fstype==-4){ // Count jets
          if (Vi1!=0 && PDGHelpers::isAJet(Vi1->id)) partcounter++;
          if (Vi2!=0 && PDGHelpers::isAJet(Vi2->id)) partcounter++;
          if (Vj1!=0 && PDGHelpers::isAJet(Vj1->id)) partcounter++;
          if (Vj2!=0 && PDGHelpers::isAJet(Vj2->id)) partcounter++;
        }
        if (partcounter<2) createCandidate=false;
      }
      if (!createCandidate) continue;

      if (debugVars::debugFlag){
        if (Vi1!=0) std::cout << "Vi1 not zero. Id: " << Vi1->id << std::endl;
        if (Vi2!=0) std::cout << "Vi2 not zero. Id: " << Vi2->id << std::endl;
        if (Vj1!=0) std::cout << "Vj1 not zero. Id: " << Vj1->id << std::endl;
        if (Vj2!=0) std::cout << "Vj2 not zero. Id: " << Vj2->id << std::endl;
      }

      TLorentzVector pH(0, 0, 0, 0);
      if (Vi1!=0) pH = pH + Vi1->p4;
      if (Vi2!=0) pH = pH + Vi2->p4;
      if (Vj1!=0) pH = pH + Vj1->p4;
      if (Vj2!=0) pH = pH + Vj2->p4;
      MELACandidate* cand = new MELACandidate(25, pH, true);

      if (Vi1!=0) cand->addDaughter(Vi1);
      if (Vi2!=0) cand->addDaughter(Vi2);
      if (Vj1!=0) cand->addDaughter(Vj1);
      if (Vj2!=0) cand->addDaughter(Vj2);

      TVar::CandidateDecayMode defaultHDecayMode = HDecayMode;
      if (isZZ==0) setCandidateDecayMode(TVar::CandidateDecay_WW);
      else if (isZZ==1) setCandidateDecayMode(TVar::CandidateDecay_ZZ);
      else if (isZZ==3) setCandidateDecayMode(TVar::CandidateDecay_ZG);
      else if (isZZ==4) setCandidateDecayMode(TVar::CandidateDecay_GG);
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
          if (partV->hasDaughter(quark)) continue;

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
      if (partV->hasDaughter(quark)) continue;
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

MELACandidate* MELAEvent::getCandidate(int index)const{
  if ((int)candidates.size()>index) return candidates.at(index);
  else return 0;
}
MELATopCandidate_t* MELAEvent::getTopCandidate(int index)const{
  if ((int) topcandidates.size()>index) return topcandidates.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getLepton(int index)const{
  if ((int)leptons.size()>index) return leptons.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getNeutrino(int index)const{
  if ((int)neutrinos.size()>index) return neutrinos.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getPhoton(int index)const{
  if ((int)photons.size()>index) return photons.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getJet(int index)const{
  if ((int)jets.size()>index) return jets.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getMother(int index)const{
  if ((int) mothers.size()>index) return mothers.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getIntermediate(int index)const{
  if ((int)intermediates.size()>index) return intermediates.at(index);
  else return 0;
}
MELAParticle* MELAEvent::getParticle(int index)const{
  if ((int)particles.size()>index) return particles.at(index);
  else return 0;
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
    for (std::vector<MELAParticle*>::iterator iL = leptons.begin(); iL!=leptons.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedLeptons(*iL); }
    for (std::vector<MELAParticle*>::iterator iL = neutrinos.begin(); iL!=neutrinos.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedNeutrinos(*iL); }
    for (std::vector<MELAParticle*>::iterator iP = photons.begin(); iP!=photons.end(); iP++){ if ((*iP)->passSelection) (*it)->addAssociatedPhotons(*iP); }
    for (std::vector<MELAParticle*>::iterator iJ = jets.begin(); iJ!=jets.end(); iJ++){ if ((*iJ)->passSelection) (*it)->addAssociatedJets(*iJ); }
    (*it)->addAssociatedVs();
    for (std::vector<MELATopCandidate_t*>::iterator iT = topcandidates.begin(); iT!=topcandidates.end(); iT++){ if ((*iT)->passSelection) (*it)->addAssociatedTops(*iT); }
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
