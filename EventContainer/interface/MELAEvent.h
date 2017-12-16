#ifndef EVENTBASE_H
#define EVENTBASE_H

#include <vector>
#include "TLorentzVector.h"
#include "ParticleComparators.h"
#include "MELACandidate.h"


class MELAEvent{
protected:
  std::vector<MELAParticle*> particles;
  std::vector<MELAParticle*> intermediates;
  std::vector<MELAParticle*> leptons;
  std::vector<MELAParticle*> neutrinos;
  std::vector<MELAParticle*> photons;
  std::vector<MELAParticle*> jets;
  std::vector<MELACandidate*> candidates;

  double xsec;
  double weight;
  std::vector<double> extraWeight;

public:

  // Constructors

  MELAEvent() : xsec(0), weight(0){}
  ~MELAEvent(){ wipeAll(); }

  // Data

  // Member functions
  void setWeight(double weight_){ weight=weight_; }
  void addExtraWeight(double weight_){ extraWeight.push_back(weight_); }
  int getNExtraWeights() const{ return extraWeight.size(); }
  double getWeight(int index=-1) const{ if (index==-1) return weight; else if (index<getNExtraWeights()) return extraWeight.at(index); else return 0; }
  std::vector<double>& getExtraWeights(){ return extraWeight; }
  void setXSec(double xsec_){ xsec=xsec_; }
  double getXSec() const{ return xsec; }


  void constructVVCandidates(int isZZ=1, int fstype=0);
  void applyParticleSelection();
  void addVVCandidateMother(MELAParticle* mother);
  void addVVCandidateAppendages();


  int getNCandidates() const{ return candidates.size(); }
  int getNLeptons() const{ return leptons.size(); }
  int getNNeutrinos() const{ return neutrinos.size(); }
  int getNPhotons() const{ return photons.size(); }
  int getNJets() const{ return jets.size(); }
  int getNIntermediates() const{ return intermediates.size(); }
  int getNParticles() const{ return particles.size(); }

  MELACandidate* getMELACandidate(int index) const;
  MELAParticle* getLepton(int index) const;
  MELAParticle* getNeutrino(int index) const;
  MELAParticle* getPhoton(int index) const;
  MELAParticle* getJet(int index) const;
  MELAParticle* getIntermediate(int index) const;
  MELAParticle* getParticle(int index) const;

  const std::vector<MELACandidate*>& getCandidates() const{ return candidates; }
  const std::vector<MELAParticle*>& getLeptons() const{ return leptons; }
  const std::vector<MELAParticle*>& getNeutrinos() const{ return neutrinos; }
  const std::vector<MELAParticle*>& getPhotons() const{ return photons; }
  const std::vector<MELAParticle*>& getJets() const{ return jets; }
  const std::vector<MELAParticle*>& getIntermediates() const{ return intermediates; }
  const std::vector<MELAParticle*>& getParticles() const{ return particles; }

  std::vector<MELACandidate*>& getCandidates(){ return candidates; }
  std::vector<MELAParticle*>& getLeptons(){ return leptons; }
  std::vector<MELAParticle*>& getNeutrinos(){ return neutrinos; }
  std::vector<MELAParticle*>& getPhotons(){ return photons; }
  std::vector<MELAParticle*>& getJets(){ return jets; }
  std::vector<MELAParticle*>& getIntermediates(){ return intermediates; }
  std::vector<MELAParticle*>& getParticles(){ return particles; }

  void addParticle(MELAParticle* myParticle){ particles.push_back(myParticle); }
  void addIntermediate(MELAParticle* myParticle){ intermediates.push_back(myParticle); }
  void addLepton(MELAParticle* myParticle, bool genuineParticle=true){ leptons.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  void addNeutrino(MELAParticle* myParticle, bool genuineParticle=true){ neutrinos.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  void addPhoton(MELAParticle* myParticle, bool genuineParticle=true){ photons.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  void addJet(MELAParticle* myParticle, bool genuineParticle=true){ jets.push_back(myParticle); if (genuineParticle) addParticle(myParticle); }
  TLorentzVector missingP() const;

protected:
  void addMELACandidate(MELACandidate*& myParticle); // Protected to avoid adding external MELACandidates and DELETING THEM TWICE!

  template<typename ParticleType> void wipeArray(std::vector<ParticleType*>& particleArray, bool doDelete=true){ if (doDelete){ for (unsigned int i=0; i<particleArray.size(); i++){ ParticleType* delpar = particleArray.at(i); delete delpar; } } particleArray.clear(); };
  void wipeAll(){ leptons.clear(); neutrinos.clear(); photons.clear(); jets.clear(); wipeArray(candidates, true); wipeArray(intermediates, false); wipeArray(particles, false); };

  void applyLeptonSelection();
  void applyNeutrinoSelection();
  void applyPhotonSelection();
  void applyJetSelection();
  void applyZZSelection();

};


#endif
