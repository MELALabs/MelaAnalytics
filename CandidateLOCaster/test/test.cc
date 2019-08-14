#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include "TVar.hh"
#include "PDGHelpers.h"
#include "Mela.h"
#include "MELACandidateRecaster.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


void test_VBF_Recast(){
  constexpr size_t nparts=9;
  int id_array[nparts]={ 1, 1, 12, -11, 13, -14, 1, 1, 21 };
  PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_WW);

  TLorentzVector mom_array[nparts] ={
    TLorentzVector(0 , 0 , 1377.14 , 1377.14),
    TLorentzVector(0 , 0 , -1078.06 , 1078.06),
    TLorentzVector(120.959 , -478.962 , 156.433 , 518.177),
    TLorentzVector(-3.72882 , -178.57 , 59.3964 , 188.227),
    TLorentzVector(-174.451 , 446.137 , -246.468 , 538.718),
    TLorentzVector(-30.6238 , 211.395 , -134.788 , 252.573),
    TLorentzVector(59.9463 , -20.8321 , 707.854 , 710.693),
    TLorentzVector(23.9593 , 34.5551 , -243.284 , 246.891),
    TLorentzVector(3.93967 , -13.7246 , -1375.52 , 1375.6)
  };


  std::vector<MELAParticle> aparticles; aparticles.reserve(3);
  std::vector<MELAParticle> daughters; daughters.reserve(4);
  std::vector<MELAParticle> mothers; mothers.reserve(2);
  for (size_t i=0; i<2; i++) mothers.emplace_back(id_array[i], mom_array[i]);
  for (size_t i=2; i<6; i++) daughters.emplace_back(id_array[i], mom_array[i]);
  for (size_t i=6; i<nparts; i++) aparticles.emplace_back(id_array[i], mom_array[i]);

  MELACandidate cand(25, true);
  for (MELAParticle& part:daughters){ part.setGenStatus(1); cand += &part; }
  cand.sortDaughters();
  for (MELAParticle& part:mothers){ part.setGenStatus(-1); cand.addMother(&part); }
  for (MELAParticle& part:aparticles){ part.setGenStatus(1); cand.addAssociatedJet(&part); }
  cand.addAssociatedVs();

  MELAout << "Original candidate:" << endl;
  MELAout << cand << endl;

  MELACandidateRecaster recaster(TVar::JJVBF);
  MELACandidate* candModified=nullptr;
  recaster.copyCandidate(&cand, candModified);

  MELAout << "Before recast:" << endl;
  MELAout << *candModified << endl;

  recaster.reduceJJtoQuarks(candModified);

  MELAout << "After recast:" << endl;
  MELAout << *candModified << endl;

  delete candModified;
}
