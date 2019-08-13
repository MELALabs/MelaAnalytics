#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "MELACandidateRecaster.h"
#include "MELAStreamHelpers.hh"

using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;


double MELACandidateRecaster::getQGMergeScore(MELAParticle* p1, MELAParticle* p2) const{
  if (
    p1 && p2 && p1!=p2
    &&
    !(PDGHelpers::isAQuark(p1->id) && PDGHelpers::isAQuark(p2->id))
    &&
    !(PDGHelpers::isAGluon(p1->id) && PDGHelpers::isAGluon(p2->id))
    ) return fabs(p1->dot(p2))*GeVsqunit;
  else return -1;
}
double MELACandidateRecaster::getGGMergeScore(MELAParticle* g1, MELAParticle* g2) const{
  if (
    g1 && g2 && g1!=g2
    &&
    PDGHelpers::isAGluon(g1->id) && PDGHelpers::isAGluon(g2->id)
    ){
    double sDiff = (g1->p4 - g2->p4).M2()*GeVsqunit;
    double sSum = (g1->p4 + g2->p4).M2()*GeVsqunit;
    return fabs(pow(sSum, 2)/sDiff);
  }
  else return -1;
}
double MELACandidateRecaster::getMergeScore(MELAParticle* p1, MELAParticle* p2) const{
  double res = getGGMergeScore(p1, p2);
  if (res<0.) res = getQGMergeScore(p1, p2);
  return res;
}



double MELACandidateRecaster::getVffEquivalentCoupling(int iferm, int jferm){
  const double xw = 0.23119;
  const double aR_lep = -2.*xw*(-1.);
  const double aL_lep = -2.*xw*(-1.)-1.;
  const double aR_neu = 0.;
  const double aL_neu = 1.;
  const double aR_QUp = -2.*xw*(2./3.);
  const double aL_QUp = -2.*xw*(2./3.)+1.;
  const double aR_QDn = -2.*xw*(-1./3.);
  const double aL_QDn = -2.*xw*(-1./3.)-1.;
  const double bL = sqrt(2.*(1.-xw));
  const double aSqLep = pow(aL_lep, 2)+pow(aR_lep, 2);
  const double aSqNu = pow(aL_neu, 2)+pow(aR_neu, 2);
  const double aSqUp = pow(aL_QUp, 2)+pow(aR_QUp, 2);
  const double aSqDn = pow(aL_QDn, 2)+pow(aR_QDn, 2);
  const double bSqWff = pow(bL, 2);
  const double Vqqcoupl[5][5]={
    { aSqDn, bSqWff*pow(0.974, 2), 0., bSqWff*pow(0.220, 2), 0. },
    { bSqWff*pow(0.974, 2), aSqUp, bSqWff*pow(0.225, 2), 0., bSqWff },
    { 0., bSqWff*pow(0.225, 2), aSqDn, bSqWff*pow(0.995, 2), 0. },
    { bSqWff*pow(0.220, 2), 0., bSqWff*pow(0.995, 2), aSqUp, bSqWff*pow(0.041, 2) },
    { 0., bSqWff*pow(0.004, 2), 0., bSqWff*pow(0.041, 2), aSqDn },
  };
  const double Vllcoupl[2][2]={
    { aSqLep, bSqWff },
    { bSqWff, aSqNu }
  };
  int absiferm=std::abs(iferm);
  int absjferm=std::abs(jferm);
  if (PDGHelpers::isAQuark(absiferm) && PDGHelpers::isAQuark(absjferm)){
    absiferm--;
    absjferm--;
    if (absiferm==5 || absjferm==5) return 0.;
    else return Vqqcoupl[absiferm][absjferm];
  }
  else if ((PDGHelpers::isALepton(absiferm) || PDGHelpers::isANeutrino(absiferm)) && (PDGHelpers::isALepton(absjferm) || PDGHelpers::isANeutrino(absjferm))){
    if (PDGHelpers::isALepton(absiferm)) absiferm=0;
    else if (PDGHelpers::isANeutrino(absiferm)) absiferm=1;
    if (PDGHelpers::isALepton(absjferm)) absjferm=0;
    else if (PDGHelpers::isANeutrino(absjferm)) absjferm=1;
    return Vllcoupl[absiferm][absjferm];
  }
  else return 0.;
}

void MELACandidateRecaster::readCandidate(
  MELACandidate* cand,
  SimpleParticleCollection_t& mothers,
  SimpleParticleCollection_t& daughters,
  SimpleParticleCollection_t& associated
  ){
  mothers.clear();
  daughters.clear();
  associated.clear();
  for (int ip=0; ip<cand->getNMothers(); ip++){
    MELAParticle* part = cand->getMother(ip);
    if (!part->passSelection) continue;
    mothers.push_back(
      SimpleParticle_t(part->id, part->p4)
      );
  }
  if (cand->getDecayMode()==TVar::CandidateDecay_Stable){
    daughters.push_back(
      SimpleParticle_t(cand->id, cand->p4)
    );
  }
  else{
    for (auto& part:cand->getSortedDaughters()){
      if (!part->passSelection) continue;
      daughters.push_back(
        SimpleParticle_t(part->id, part->p4)
      );
    }
  }
  for (auto& part:cand->getAssociatedLeptons()){
    if (!part->passSelection) continue;
    associated.push_back(
      SimpleParticle_t(part->id, part->p4)
      );
  }
  for (auto& part:cand->getAssociatedPhotons()){
    if (!part->passSelection) continue;
    associated.push_back(
      SimpleParticle_t(part->id, part->p4)
      );
  }
  for (auto& part:cand->getAssociatedJets()){
    if (!part->passSelection) continue;
    if (part->genStatus==1){
      associated.push_back(
        SimpleParticle_t(part->id, part->p4)
        );
    }
    else{
      mothers.push_back(
        SimpleParticle_t(part->id, part->p4)
        );
    }
  }
}

MELAParticle* MELACandidateRecaster::getBestAssociatedV(MELACandidate* cand, TVar::Production production, double* bestScore){
  MELAParticle* protectV=nullptr;
  if (bestScore!=nullptr) *bestScore=0;
  bool useZH=(production==TVar::Had_ZH || production==TVar::Lep_ZH);
  bool useWH=(production==TVar::Had_WH || production==TVar::Lep_WH);
  if (useZH || useWH){
    int nsortedvsstart=(cand->getDecayMode()!=TVar::CandidateDecay_Stable ? 2 : 0);
    double bestVLepscore=-1;
    int bestAssociatedVLepindex=-1;
    double bestVHadscore=-1;
    int bestAssociatedVHadindex=-1;
    for (int iv=nsortedvsstart; iv<cand->getNSortedVs(); iv++){
      MELAParticle* sortedV = cand->getSortedV(iv);
      if (sortedV!=nullptr && ((PDGHelpers::isAZBoson(sortedV->id) && useZH) || (PDGHelpers::isAWBoson(sortedV->id) && useWH))){
        //MELAout << "- Considering sorted V " << iv << endl;
        //MELAout << sortedV->getDaughter(0)->id << " " << sortedV->getDaughter(0)->passSelection << endl;
        //MELAout << sortedV->getDaughter(1)->id << " " << sortedV->getDaughter(1)->passSelection << endl;
        if (sortedV->getNDaughters()==2 && sortedV->getDaughter(0)->passSelection && sortedV->getDaughter(1)->passSelection){
          double s = sortedV->p4.M2()*GeVsqunit;
          double m, ga;
          if (useZH) TUtil::GetMassWidth(23, m, ga); else TUtil::GetMassWidth(24, m, ga);
          m *= GeVunit; ga *= GeVunit;
          double score = pow(s-pow(m, 2), 2)+pow(m*ga, 2)/getVffEquivalentCoupling(sortedV->getDaughter(0)->id, sortedV->getDaughter(1)->id); score = 1./score;
          if (sortedV->getDaughter(0) && (PDGHelpers::isALepton(sortedV->getDaughter(0)->id) || PDGHelpers::isANeutrino(sortedV->getDaughter(0)->id))){
            if (score>bestVLepscore || bestVLepscore<0.){
              bestVLepscore=score;
              bestAssociatedVLepindex=iv;
            }
          }
          else if (sortedV->getDaughter(0)){
            if (score>bestVHadscore || bestVHadscore<0.){
              bestVHadscore=score;
              bestAssociatedVHadindex=iv;
            }
          }
        }
      }
    }
    if (bestAssociatedVLepindex>=0){
      protectV = cand->getSortedV(bestAssociatedVLepindex);
      if (bestScore!=nullptr) *bestScore=bestVLepscore;
    }
    else if (bestAssociatedVHadindex>=0){
      protectV = cand->getSortedV(bestAssociatedVHadindex);
      if (bestScore!=nullptr) *bestScore=bestVHadscore;
    }
  }
  return protectV;
}

void MELACandidateRecaster::adjustForIncomingMomenta(
  SimpleParticleCollection_t& mothers,
  SimpleParticleCollection_t& daughters,
  SimpleParticleCollection_t& associated
  ){
  if (mothers.size()!=2) return;

  TLorentzRotation ltr; // Boost*Rotation*Boost
  TLorentzVector const& pLab1 = mothers.at(0).second;
  TLorentzVector const& pLab2 = mothers.at(1).second;
  TLorentzVector pLab12 = pLab1 + pLab2;
  if (pLab1.Pt()+pLab2.Pt()==0.) return; // If pTs are 0, don't touch the system
  // Else, boost to p=0 frame and rotate the z axis

  // First boost to p=0 frame
  TVector3 pLabBoost = -pLab12.BoostVector();
  ltr.Boost(pLabBoost);
  TLorentzVector pBoosted1 = pLab1; pBoosted1.Boost(pLabBoost);
  TLorentzVector pBoosted2 = pLab2; pBoosted2.Boost(pLabBoost);

  // Determine the z' axis such that z'.z>0
  TVector3 zprimeWU = (pBoosted1-pBoosted2).Vect();
  if (zprimeWU.Z()<0.) zprimeWU=-zprimeWU;
  if (zprimeWU.Mag()>0.){
    TVector3 zprime = zprimeWU.Unit();

    // Rotate the system to align with the z axis
    TVector3 zaxis(0, 0, 1);
    double angle = acos(zprime.Z());
    ltr.Rotate(angle, zprime.Cross(zaxis)); // This cross product - angle combination is correct to rotate back!
  }

  for (auto& part:mothers) part.second = ltr*part.second;
  for (auto& part:daughters) part.second = ltr*part.second;
  for (auto& part:associated) part.second = ltr*part.second;
}



MELACandidateRecaster::MELACandidateRecaster(const TVar::Production candScheme_) :
candScheme(candScheme_),
doDotincomingParticles(false),
protectVStrict(false),
nAssociated(0)
{
  if (candScheme==TVar::JJVBF_S || candScheme==TVar::JJVBF_TU) candScheme=TVar::JJVBF;
  else if (candScheme==TVar::Had_ZH_S || candScheme==TVar::Had_ZH_TU) candScheme=TVar::Had_ZH;
  else if (candScheme==TVar::Had_WH_S || candScheme==TVar::Had_WH_TU) candScheme=TVar::Had_WH;
  else if (candScheme==TVar::Lep_ZH_S || candScheme==TVar::Lep_ZH_TU) candScheme=TVar::Lep_ZH;
  else if (candScheme==TVar::Lep_WH_S || candScheme==TVar::Lep_WH_TU) candScheme=TVar::Lep_WH;

  if (
    candScheme==TVar::JJVBF
    ||
    candScheme==TVar::JJQCD
    ||
    candScheme==TVar::Had_ZH
    ||
    candScheme==TVar::Had_WH
    ||
    candScheme==TVar::Lep_ZH
    ||
    candScheme==TVar::Lep_WH
    ){
    doDotincomingParticles=true;
    nAssociated=2;
  }

  /*
  protectVStrict = (
  candScheme==TVar::Had_ZH
  ||
  candScheme==TVar::Had_WH
  ||
  candScheme==TVar::Lep_ZH
  ||
  candScheme==TVar::Lep_WH
  );
  */
}

MELACandidateRecaster::~MELACandidateRecaster(){ clearExtraParticles(); }

MELAParticle* MELACandidateRecaster::mergeTwoGluons(MELAParticle* glu1, MELAParticle* glu2){
  MELAParticle* mergedGluon = new MELAParticle(21);
  mergedGluon->genStatus=2;
  mergedGluon->setSelected(false);
  if (glu1->genStatus==-1) mergedGluon->p4 -= glu1->p4;
  else mergedGluon->p4 += glu1->p4;
  if (glu2->genStatus==-1) mergedGluon->p4 -= glu2->p4;
  else mergedGluon->p4 += glu2->p4;
  mergedGluon->addDaughter(glu1);
  mergedGluon->addDaughter(glu2);
  return mergedGluon;
}

bool MELACandidateRecaster::merge2Qto1G(
  MELACandidate*& cand,
  std::vector<MELAParticle*>& gluons,
  std::vector<MELAParticle*>& quarks,
  MELAParticle*& protectV
  ){
  int ifound=-1; int jfound=-1; double minqsq=-1;
  std::vector<MELAParticle*> allquarks(quarks);
  if (protectV){
    for (auto& Vdau : protectV->getDaughters()){
      if (PDGHelpers::isAQuark(Vdau->id) && !MELAParticle::checkParticleExists(Vdau, quarks)) allquarks.push_back(Vdau);
    }
  }
  //MELAout << "quarks size = " << allquarks.size() << endl;
  //for (auto& q:allquarks) MELAout << " - Quark id = " << q->id << " , E = " << q->t() << endl;
  for (unsigned int i = 0; i<allquarks.size(); i++){
    MELAParticle* quark_i = allquarks.at(i);
    double Qi = quark_i->charge();
    int id_i = quark_i->id;
    TLorentzVector p4_i=quark_i->p4;
    if (quark_i->genStatus==-1){ id_i *= -1; Qi *= -1.; p4_i = -p4_i; }

    for (unsigned int j = i+1; j<allquarks.size(); j++){
      MELAParticle* quark_j = allquarks.at(j);
      double Qj = quark_j->charge();
      int id_j = quark_j->id;
      TLorentzVector p4_j=quark_j->p4;
      if (quark_j->genStatus==-1){ id_j *= -1; Qj *= -1.; p4_j = -p4_j; }
      //MELAout << "- Checking i, j = " << i << " , " << j << endl;
      //MELAout << " - Qi, Qj = " << Qi << " , " << Qj << endl;
      //MELAout << " - id_i, id_j = " << id_i << " , " << id_j << endl;

      MELAParticle* tmpProtV = nullptr;
      double tmpVpropagator=1;
      if (protectV){
        // Temporarily set selection=false to test whether there is another suitable V to be found
        quark_i->setSelected(false);
        quark_j->setSelected(false);
        tmpProtV=getProtectedV(cand, &tmpVpropagator);
        if (tmpVpropagator!=0.) tmpVpropagator = 1./tmpVpropagator; // Need to invert propagator to find minimum score
        quark_i->setSelected(true);
        quark_j->setSelected(true);
      }
      if (
        ((tmpProtV && protectV) || !protectV)
        &&
        (((Qi+Qj)==0 && (id_i+id_j)==0) || (PDGHelpers::isAnUnknownJet(id_i) || PDGHelpers::isAnUnknownJet(id_j)))
        ){
        TLorentzVector pTotal = p4_i+p4_j;
        double qsq = fabs(pTotal.M2())*GeVsqunit;
        // If a protectV exists, multiply by the propagator as well
        if (protectV) qsq *= tmpVpropagator;
        if (minqsq<0. || minqsq>qsq){
          minqsq=qsq;
          ifound=i;
          jfound=j;
          //MELAout << "tmpProtV 1 E = " << tmpProtV->getDaughter(0)->t() << endl;
          //MELAout << "tmpProtV 2 E = " << tmpProtV->getDaughter(1)->t() << endl;
        }
      }
    }
  }
  //MELAout << "ifound, jfound = " << ifound << " " << jfound << endl;
  //MELAout << "protectV = " << (protectV ? protectV->id : -9000) << endl;
  if (ifound>=0 && jfound>=0){
    MELAParticle* quark_i = allquarks.at(ifound); quark_i->setSelected(false);
    MELAParticle* quark_j = allquarks.at(jfound); quark_j->setSelected(false);
    MELAParticle* newGluon = quark_i; quark_j->id=-9000;
    TLorentzVector newGluon_p4(0, 0, 0, 0);
    int newGluon_id=21;
    int newGluon_genStatus=1;
    if (quark_i->genStatus==-1) newGluon_p4 = newGluon_p4 - quark_i->p4;
    else newGluon_p4 = newGluon_p4 + quark_i->p4;
    if (quark_j->genStatus==-1) newGluon_p4 = newGluon_p4 - quark_j->p4;
    else newGluon_p4 = newGluon_p4 + quark_j->p4;
    newGluon->p4=newGluon_p4;
    newGluon->id=newGluon_id;
    newGluon->genStatus=newGluon_genStatus;
    newGluon->setSelected(false);
    gluons.push_back(newGluon);
    // Remove the two merged quarks from the quarks list for permutations
    allquarks.erase(allquarks.begin()+jfound);
    allquarks.erase(allquarks.begin()+ifound);
    std::swap(allquarks, quarks);
    return true;
  }
  return false;
}

std::vector<pair<MELAParticle*, std::vector<MELAParticle*>>> MELACandidateRecaster::mapGluonsToQuarks(
  const std::vector<MELAParticle*>& gluons,
  const std::vector<MELAParticle*>& quarks,
  const std::vector<int>& permutation
  ){
  assert(gluons.size()==permutation.size());

  vector<pair<MELAParticle*, vector<MELAParticle*>>> quark_gluon_collection;
  for (int iq=0; iq<(int)quarks.size(); iq++){
    MELAParticle* quark = quarks.at(iq);
    quark_gluon_collection.push_back(
      pair<MELAParticle*, vector<MELAParticle*>>(quark, vector<MELAParticle*>())
      );
    for (unsigned int ig=0; ig<gluons.size(); ig++){
      MELAParticle* gluon = gluons.at(ig);
      const int& qindex = permutation.at(ig);
      if (qindex==iq) quark_gluon_collection.back().second.push_back(gluon);
    }
  }

  return quark_gluon_collection;
}

double MELACandidateRecaster::getBestVHConfig(
  MELAParticle* protectV,
  const std::vector<pair<MELAParticle*, std::vector<MELAParticle*>>>& quarks,
  std::vector<int>* qordered,
  int* swapconfig
  ){
  const unsigned int nQreq=2;
  double bestCfgLikelihood=-1;
  double bestCfgSysE=0;
  if (quarks.size()<nQreq || protectV==nullptr) return bestCfgLikelihood;

  if (swapconfig) *swapconfig=0;
  if (qordered){
    qordered->clear();
    for (unsigned int i=0; i<nQreq; i++) qordered->push_back(-1);
  }

  vector<vector<int>> qperm;
  CombinationGenerator(quarks.size(), nQreq, qperm, 0);
  for (auto& perm:qperm){
    bool skipPermutation=false;
    int swapcfg=0;
    int id[nQreq];
    TLorentzVector p4[nQreq];

    // Assign incoming momenta
    for (unsigned int i=0; i<nQreq; i++){
      const int& qindex = perm.at(i);
      MELAParticle* quark=quarks.at(qindex).first;
      if (MELAParticle::checkParticleExists(quark, protectV->getDaughters())) skipPermutation=true;
      if (skipPermutation) break;

      int genstat = quark->genStatus;
      id[i]=quark->id;
      p4[i]=quark->p4;
      for (MELAParticle* gluon:quarks.at(qindex).second){
        if (gluon->genStatus*genstat>0) p4[i] = p4[i] + gluon->p4;
        else p4[i] = p4[i] - gluon->p4;
      }
      if (genstat>=0){
        id[i] = -id[i];
        p4[i] = -p4[i];
        swapcfg += pow(2, i);
      }
    }

    if (skipPermutation) continue;

    // Check if boosting to the pT=0 frame is still valid (MELA check for PDFs)
    TLorentzVector pTotal=p4[0]+p4[1];
    if (pTotal.M2()<0.) continue;
    pTotal.Boost(-pTotal.X()/pTotal.T(), -pTotal.Y()/pTotal.T(), 0);
    const double sysPz = pTotal.Z();
    const double sysE = pTotal.T();
    const double Ez0 = (sysE+sysPz)/2.;
    const double Ez1 = (sysE-sysPz)/2.;
    if (Ez0<0. || Ez1<0.) continue;

    // Couplings
    double couplings_allcfgs=0;
    // Obtain q1-q2 couplings
    if (PDGHelpers::getCoupledVertex(id[0], id[1])==protectV->id){
      //MELAout << "id 0, 1 = " << id[0] << " , " << id[1] << endl;
      couplings_allcfgs = getVffEquivalentCoupling(id[0], id[1]);
    }
    else continue;

    // Do not check for (bestCfgLikelihood<0. && couplings_allcfgs>0.) because we may want to let ME->0 without PDF errors
    // Also check for a greater sysE if the current combination likelihood is equivalent.
    if (bestCfgLikelihood<0. || bestCfgLikelihood<couplings_allcfgs || (bestCfgSysE<sysE && bestCfgLikelihood==couplings_allcfgs)){
      bestCfgLikelihood=couplings_allcfgs;
      bestCfgSysE=sysE;
      if (qordered){
        *qordered=perm;
        for (int iquark=0; iquark<(int)quarks.size(); iquark++){
          bool indexPresent=false;
          for (auto& jq:(*qordered)){ if (iquark==jq){ indexPresent=true; break; } }
          if (!indexPresent) qordered->push_back(iquark);
        }
      }
      if (swapconfig) *swapconfig=swapcfg;
    }
  }
  //if (bestCfgLikelihood<0.) MELAerr << "Ordering failed!" << endl;
  if (bestCfgLikelihood<0. && qordered) qordered->clear(); // Clear so that we can just check for size
  //else if (qordered){
  //  MELAout << "Final VH order:";
  //  for (auto& ord:*qordered) MELAout << ord << " (" << quarks.at(ord).first->id << ") ";
  //  MELAout << endl;
  //}
  return bestCfgLikelihood;
}

double MELACandidateRecaster::getBestVBFConfig(
  const std::vector<pair<MELAParticle*, std::vector<MELAParticle*>>>& quarks,
  std::vector<int>* qordered,
  int* swapconfig
  ){
  constexpr unsigned int nQin=2;
  constexpr unsigned int nQreq=4;
  double bestCfgLikelihood=-1;
  double bestCfgSysE=0;
  if (quarks.size()<nQreq) return bestCfgLikelihood;

  if (swapconfig) *swapconfig=0;
  if (qordered){
    qordered->clear();
    for (unsigned int i=0; i<nQreq; i++) qordered->push_back(-1);
  }

  vector<vector<int>> qperm;
  {
    vector<vector<int>> qperm2;
    CombinationGenerator(quarks.size(), nQin, qperm2, 0);
    for (auto it_in=qperm2.begin(); it_in!=qperm2.end(); it_in++){
      for (auto it_out=qperm2.begin(); it_out!=qperm2.end(); it_out++){
        if (it_in==it_out) continue;
        vector<int>& inperm = *it_in;
        vector<int>& outperm = *it_out;
        bool isUniqueComb=true;
        for (auto& ppin:inperm){
          for (auto& ppout:outperm){
            if (ppin==ppout){
              isUniqueComb=false;
              break;
            }
          }
          if (!isUniqueComb) break;
        }
        if (isUniqueComb){
          qperm.push_back(vector<int>());
          for (auto& ppin:inperm) qperm.back().push_back(ppin);
          for (auto& ppout:outperm) qperm.back().push_back(ppout);
        }
      }
    }
  }
  for (auto& perm:qperm){
    int swapcfg=0;
    int id[nQreq];
    TLorentzVector p4[nQreq];

    // Assign incoming/outgoing momenta
    for (unsigned int i=0; i<nQreq; i++){
      const int& qindex = perm.at(i);
      const MELAParticle* quark=quarks.at(qindex).first;

      int genstat = quark->genStatus;
      id[i]=quark->id;
      p4[i]=quark->p4;
      for (MELAParticle* gluon:quarks.at(qindex).second){
        if (gluon->genStatus*genstat>0) p4[i] = p4[i] + gluon->p4;
        else p4[i] = p4[i] - gluon->p4;
      }
      if ((genstat==-1 && i>=2) || (genstat>=0 && i<2)){
        id[i] = -id[i];
        p4[i] = -p4[i];
        swapcfg += pow(2, i);
      }
    }

    // Check if boosting to the pT=0 frame is still valid (MELA check for PDFs)
    TLorentzVector pTotal=p4[0]+p4[1];
    if (pTotal.M2()<0.) continue;
    pTotal.Boost(-pTotal.X()/pTotal.T(), -pTotal.Y()/pTotal.T(), 0);
    const double sysPz = pTotal.Z();
    const double sysE = pTotal.T();
    const double Ez0 = (sysE+sysPz)/2.;
    const double Ez1 = (sysE-sysPz)/2.;
    if (Ez0<0. || Ez1<0.) continue;

    // Couplings
    double couplsq13=0, couplsq24=0, couplsq23=0, couplsq14=0;
    // Obtain q1-q3 and q2-q4 couplings
    if (TMath::Sign(1, id[0])==TMath::Sign(1, id[2]) && TMath::Sign(1, id[1])==TMath::Sign(1, id[3]) && PDGHelpers::getCoupledVertex(id[0], -id[2])==PDGHelpers::getCoupledVertex(-id[1], id[3])){ // Incoming 13, outgoing 24
      couplsq13 = getVffEquivalentCoupling(id[0], id[2]);
      couplsq24 = getVffEquivalentCoupling(id[1], id[3]);
    }
    // Obtain q1-q4 and q2-q3 couplings
    if (TMath::Sign(1, id[1])==TMath::Sign(1, id[2]) && TMath::Sign(1, id[0])==TMath::Sign(1, id[3]) && PDGHelpers::getCoupledVertex(id[1], -id[2])==PDGHelpers::getCoupledVertex(-id[0], id[3])){ // Incoming 23, outgoing 14
      couplsq23 = getVffEquivalentCoupling(id[1], id[2]);
      couplsq14 = getVffEquivalentCoupling(id[0], id[3]);
    }
    // Check coupling combination and assign qordered to the best known config so far
    double couplings_13_24=couplsq13*couplsq24;
    double couplings_14_23=couplsq14*couplsq23;
    double couplings_allcfgs = couplings_13_24+couplings_14_23;
    // Do not check for (bestCfgLikelihood<0. && couplings_allcfgs>0.) because we may want to let ME->0 without PDF errors
    // Also check for a greater sysE if the current combination likelihood is equivalent.
    if (bestCfgLikelihood<0. || bestCfgLikelihood<couplings_allcfgs || (bestCfgSysE<sysE && bestCfgLikelihood==couplings_allcfgs)){
      bestCfgLikelihood=couplings_allcfgs;
      bestCfgSysE=sysE;
      //MELAout << "Perm: ";
      //for (auto& pp: perm) MELAout << pp << " ";
      //MELAout << endl;
      if (qordered) swap(*qordered, perm);
      if (swapconfig) *swapconfig=swapcfg;
      //MELAout << "New VBF config ids = " << id[0] << " , " << id[1] << " , " << id[2] << " , " << id[3] << ". ";
      //MELAout << "couplsq13 = " << couplsq13 << ", couplsq24 = " << couplsq24 << ", couplsq14 = " << couplsq14 << ", couplsq23 = " << couplsq23 << endl;
      //MELAout << "Likelihood: " << bestCfgLikelihood << endl;
    }
  }
  if (bestCfgLikelihood<0. && qordered) qordered->clear(); // Clear so that we can just check for size
  return bestCfgLikelihood;
}


double MELACandidateRecaster::getBestHJJConfig(
  const std::vector<std::pair<MELAParticle*, std::vector<MELAParticle*>>>& partons,
  std::vector<int>* order,
  int* swapconfig
){
  constexpr unsigned int nPartonsIn=2;
  constexpr unsigned int nPartonsReq=4;
  bool bestCfgLikelihoodFound = false;
  double bestCfgSysE=0;
  if (partons.size()<nPartonsReq) return -1;

  if (swapconfig) *swapconfig=0;
  if (order){
    order->clear();
    for (unsigned int i=0; i<nPartonsReq; i++) order->push_back(-1);
  }

  vector<vector<int>> parton_perm;
  {
    vector<vector<int>> parton_perm2;
    CombinationGenerator(partons.size(), nPartonsIn, parton_perm2, 0);
    for (auto it_in=parton_perm2.begin(); it_in!=parton_perm2.end(); it_in++){
      for (auto it_out=parton_perm2.begin(); it_out!=parton_perm2.end(); it_out++){
        if (it_in==it_out) continue;
        vector<int>& inperm = *it_in;
        vector<int>& outperm = *it_out;
        bool isUniqueComb=true;
        for (auto& ppin:inperm){
          for (auto& ppout:outperm){
            if (ppin==ppout){
              isUniqueComb=false;
              break;
            }
          }
          if (!isUniqueComb) break;
        }
        if (isUniqueComb){
          parton_perm.push_back(vector<int>());
          for (auto& ppin:inperm) parton_perm.back().push_back(ppin);
          for (auto& ppout:outperm) parton_perm.back().push_back(ppout);
        }
      }
    }
  }
  for (auto& perm:parton_perm){
    int swapcfg=0;
    int id[nPartonsReq];
    TLorentzVector p4[nPartonsReq];

    // Assign incoming/outgoing momenta
    for (unsigned int i=0; i<nPartonsReq; i++){
      const int& qindex = perm.at(i);
      const MELAParticle* quark=partons.at(qindex).first;

      int genstat = quark->genStatus;
      id[i]=quark->id;
      p4[i]=quark->p4;
      for (MELAParticle* gluon:partons.at(qindex).second){
        if (gluon->genStatus*genstat>0) p4[i] = p4[i] + gluon->p4;
        else p4[i] = p4[i] - gluon->p4;
      }
      if ((genstat==-1 && i>=2) || (genstat>=0 && i<2)){
        id[i] = -id[i];
        p4[i] = -p4[i];
        swapcfg += pow(2, i);
      }
    }

    // Check if boosting to the pT=0 frame is still valid (MELA check for PDFs)
    TLorentzVector pTotal=p4[0]+p4[1];
    if (pTotal.M2()<0.) continue;
    pTotal.Boost(-pTotal.X()/pTotal.T(), -pTotal.Y()/pTotal.T(), 0);
    const double sysPz = pTotal.Z();
    const double sysE = pTotal.T();
    const double Ez0 = (sysE+sysPz)/2.;
    const double Ez1 = (sysE-sysPz)/2.;
    if (Ez0<0. || Ez1<0.) continue;

    // Also check for a greater sysE if the current combination likelihood is equivalent.
    if (!bestCfgLikelihoodFound || bestCfgSysE<sysE){
      bestCfgLikelihoodFound = true;
      bestCfgSysE=sysE;
      if (order) swap(*order, perm);
      if (swapconfig) *swapconfig=swapcfg;
    }
  }
  if (!bestCfgLikelihoodFound && order) order->clear(); // Clear so that we can just check for size
  return (bestCfgLikelihoodFound ? 1 : -1);
}

double MELACandidateRecaster::getVBFLikelihood(
  const std::vector<MELAParticle*>& gluons,
  const std::vector<MELAParticle*>& quarks,
  const std::vector<int>& permutation
  ){
  assert(gluons.size()==permutation.size());
  vector<pair<MELAParticle*, vector<MELAParticle*>>> quark_gluon_collection = mapGluonsToQuarks(gluons, quarks, permutation);
  return getBestVBFConfig(quark_gluon_collection);
}

double MELACandidateRecaster::getMergeOrder_GluonsIntoQuarks(
  std::vector<MELAParticle*>& gluons,
  std::vector<MELAParticle*>& quarks,
  std::vector<int>& bestConfig,
  bool doMergeGluons,
  std::vector<MELAParticle*>* VmassMonitor
  ){
  double smallestKDP=-99;
  const int nQs = quarks.size();
  const int nGs = gluons.size();
  const bool nGsIsOdd = (nGs%2==1);
  if (nGs>nQs || nGs==0) return smallestKDP;

  vector<int> bestMergedConfig;
  double bestMergedConfigKDP=smallestKDP;
  // If there are more than one gluons and the function call requests merging the gluons,
  // check whether a pairqise merged configuration has better k.p than the possible unmerged configurations.
  if (doMergeGluons && nGs>1){
    vector<vector<int>> gluonPermutations;
    // Get nGs! permutations so that the set of pairwise merges is just equivalent to splitting of these permutations pairwise.
    PermutationGenerator(nGs, nGs, gluonPermutations, 0);
    for (auto& gperm:gluonPermutations){
      assert((int)gperm.size()==nGs);
      // Get merged gluons
      vector<MELAParticle*> mergedGluons; // Temporary collection of merged gluons to delete
      for (unsigned int ip=0; ip<gperm.size()/2; ip++){
        int ig = gperm.at(2*ip);
        int jg = gperm.at(2*ip+1);
        MELAParticle* glu1=gluons.at(ig);
        MELAParticle* glu2=gluons.at(jg);
        MELAParticle* mergedGluon = mergeTwoGluons(glu1, glu2);
        mergedGluons.push_back(mergedGluon);
      }

      // Number of merging possibilities is 2^(N merged gluons)
      int nMergeUnmerge = pow(2, (int)mergedGluons.size());
      // Iterate over the different merging possiblities of the curent gluon permutation
      // Start from 1 because 0 means all unmerged, which is not what we want to check under the doMergeGluons if-condition.
      for (int imc=1; imc<nMergeUnmerge; imc++){
        int icfg=imc;
        vector<MELAParticle*> tmpGluonList;
        double penalty=1;
        for (auto& mergedGluon:mergedGluons){ // Iterate over the gluon pairs
          if (icfg%2==0){ // Do not merge this particular gluon pair
            tmpGluonList.push_back(mergedGluon->getDaughter(0));
            tmpGluonList.push_back(mergedGluon->getDaughter(1));
          }
          else{ // Merge this gluon pair
            tmpGluonList.push_back(mergedGluon);
            penalty *= getGGMergeScore(mergedGluon->getDaughter(0), mergedGluon->getDaughter(1));
          }
          icfg = icfg >> 1;
        }
        // Remember that if there is an odd number of starting gluons, the last gluon of the permutation over the gluons is a singleton.
        if (nGsIsOdd) tmpGluonList.push_back(gluons.at(gperm.back()));

        // Find the likelihood of the merged configuration considered
        vector<int> tmpCfg;
        double tmpResult = getMergeOrder_GluonsIntoQuarks(tmpGluonList, quarks, tmpCfg, false, VmassMonitor); // Do not apply further recursion because the penalty will no longer be accurate.
        if (tmpResult>0.){
          tmpResult *= penalty;
          if (bestMergedConfigKDP<0. || bestMergedConfigKDP>tmpResult){
            bestMergedConfigKDP=tmpResult;
            bestMergedConfig.clear();
            assert(tmpCfg.size()==tmpGluonList.size());
            // We cannot simply assign bestMergedConfig=tmpCfg because the order of gluons is different and some are merged,
            // so map back to the original order.
            for (unsigned int igluon=0; igluon<gluons.size(); igluon++){ // Original order
              MELAParticle* matchGluon = gluons.at(igluon);
              for (unsigned int jtg=0; jtg<tmpGluonList.size(); jtg++){
                MELAParticle* tmpGluon=tmpGluonList.at(jtg);
                if (tmpGluon==matchGluon || MELAParticle::checkParticleExists(matchGluon, tmpGluon->getDaughters())){
                  bestMergedConfig.push_back(tmpCfg.at(jtg));
                  break;
                }
              } // Loop over tmpGluonList
            } // Loop over original gluons array
            //MELAout << "New bestmergedcfg: ";
            //for (auto& c:bestMergedConfig) MELAout << c << " ";
            //MELAout << "(k.p=" << tmpResult << ", multiplier=" << multiplier << ")" << endl;
          }

        }
      }
      for (auto& mg:mergedGluons) delete mg;
    } // Loop over gluon permutations
  } // doMergeGluons?

  vector<vector<float>> quantifiers; // Size of each quantifier should be the same as quarks, and size of the double-vector should be nGs
  for (auto& gluon:gluons){
    vector<float> theQuantifier;
    // Calculate k.p between quarks and gluons
    for (auto& quark:quarks){
      float kDp = 1e20;
      if (!(quark->genStatus==-1 && gluon->genStatus==-1) || doDotincomingParticles) kDp = getQGMergeScore(gluon, quark); // Do not dot incoming particles
      theQuantifier.push_back(kDp);
    }
    quantifiers.push_back(theQuantifier);
  }
  vector<vector<int>> permutator; PermutationGenerator(quarks.size(), gluons.size(), permutator, 0);

  int chosenPerm=-1;
  for (unsigned int iperm=0; iperm<permutator.size(); iperm++){
    vector<int>& perm = permutator.at(iperm);
    //MELAout << "Checking perm " << perm << endl;
    double theKDP=1;
    bool skipPermutation=false;
    for (unsigned int ig=0; ig<perm.size(); ig++) theKDP *= quantifiers.at(ig).at(perm.at(ig));
    if (candScheme==TVar::JJVBF){
      double VBFL= getVBFLikelihood(gluons, quarks, perm);
      if (VBFL<=0.) VBFL=1e40; // Put something huge
      theKDP *= VBFL;
    }
    if (VmassMonitor!=nullptr && !VmassMonitor->empty()){
      for (auto& part : (*VmassMonitor)){
        double m, ga;
        TUtil::GetMassWidth(part, m, ga);
        m *= GeVunit;
        ga *= GeVunit;

        TLorentzVector pV = part->p4;
        const double unmergedVE = pV.T()*GeVunit;
        const double unmergedVqsq = pV.M2()*GeVsqunit;
        for (unsigned int ig=0; ig<perm.size(); ig++){
          if (MELAParticle::checkParticleExists(quarks.at(perm.at(ig)), part->getDaughters())){
            MELAParticle* theGluon = gluons.at(ig);
            if (quarks.at(perm.at(ig))->genStatus*theGluon->genStatus>0) pV = pV + theGluon->p4;
            else{
              pV = pV - theGluon->p4;
              skipPermutation |= (theGluon->genStatus<0); // Do not flip outgoing quark sign
            }
          }
        }
        double qsq=pV.M2()*GeVsqunit;
        double VE=pV.T()*GeVunit;
        double propagator = pow(qsq-pow(m, 2), 2) + pow(m*ga, 2);
        theKDP *= propagator;
        skipPermutation |= (TMath::Sign(1., VE)!=TMath::Sign(1., unmergedVE) || VE==0. || TMath::Sign(1., qsq)!=TMath::Sign(1., unmergedVqsq) || qsq==0.); // Do not change the sign of E or Q2 of the protected V
      }
    }
    if (skipPermutation) continue;
    theKDP = fabs(theKDP);
    if (smallestKDP<0. || theKDP<smallestKDP){
      smallestKDP=theKDP;
      chosenPerm=iperm;
    }
  }

  if (chosenPerm>=0) bestConfig=permutator.at(chosenPerm);
  // If the merged config is better than split, use that one
  if (bestMergedConfigKDP>0. && bestMergedConfigKDP<smallestKDP){
    bestConfig=bestMergedConfig;
    smallestKDP=bestMergedConfigKDP;
  }
  return smallestKDP;
}


void MELACandidateRecaster::copyCandidate(MELACandidate* cand, MELACandidate*& candModified, bool adjustIncoming){
  SimpleParticleCollection_t mothersnew;
  SimpleParticleCollection_t daughtersnew;
  SimpleParticleCollection_t associatednew;
  readCandidate(cand, mothersnew, daughtersnew, associatednew);
  if (adjustIncoming) adjustForIncomingMomenta(mothersnew, daughtersnew, associatednew);

  TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
  PDGHelpers::setCandidateDecayMode(cand->getDecayMode());

  candModified = TUtil::ConvertVectorFormat(
    &daughtersnew,
    &associatednew,
    &mothersnew,
    true,
    &extraParticles, nullptr
  );

  PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
}

void MELACandidateRecaster::reduceJJtoQuarks(MELACandidate*& cand){
  if (candScheme==TVar::nProductions) return;
  if (!cand) return;
  unsigned int nQ = 2 + nAssociated; // Assume mothers are supposed to be quarks as well
  // If VH, find a V to protect.
  vector<MELAParticle*> protectedVs;
  MELAParticle* protectV = getProtectedV(cand);
  const bool hasProtectedV = (protectV!=nullptr);
  // Determine whether the protected V is leptonic.
  bool leptonicV=false;
  if (hasProtectedV && protectV->getDaughter(0)!=nullptr) leptonicV = (PDGHelpers::isALepton(protectV->getDaughter(0)->id) || PDGHelpers::isANeutrino(protectV->getDaughter(0)->id));
  if ((leptonicV || protectVStrict) && nAssociated>=2) nQ-=2; // nAssociated includes leptons, so subtract the number of outgoing leptons to get the number of quarks required.

  // Gather all quarks and gluons
  vector<MELAParticle*> gluons, quarks;
  // Gather all mother quarks and gluons
  for (int ip=0; ip<cand->getNMothers(); ip++){
    MELAParticle* part = cand->getMother(ip);
    if (PDGHelpers::isAGluon(part->id)) gluons.push_back(part);
    else quarks.push_back(part);
  }
  // Gather all associated quarks and gluons except the protected V
  for (auto& part:cand->getAssociatedJets()){
    if (PDGHelpers::isAGluon(part->id)) gluons.push_back(part);
    else{
      // Skip the quarks from the protected V (Note: Done for code protection)
      if (hasProtectedV && !leptonicV && MELAParticle::checkParticleExists(part, protectV->getDaughters())) continue;
      quarks.push_back(part);
    }
  }
  // Gather quarks from the protected V at the last step
  if (hasProtectedV && !leptonicV && !protectVStrict){ for (auto& part:cand->getAssociatedJets()){ if (MELAParticle::checkParticleExists(part, protectV->getDaughters())) quarks.push_back(part); } }

  bool hasAtLeastOneGluon = !gluons.empty();
  if (quarks.size()==nQ && !hasAtLeastOneGluon) return; // Do not continue further, you have the candidate.

  /**********************************************/
  /* DO NOT MODIFY THE OBJECTS UNTIL THIS POINT */
  /**********************************************/

  // Disable the gluons
  // FIXME: ONLY DISABLE THE GLUONS THAT GET MERGED
  for (auto& gluon:gluons) gluon->setSelected(false);

  // Merge two quarks to create a gluon
  if (quarks.size()>nQ && !hasAtLeastOneGluon){
    //TUtil::PrintCandidateSummary(cand);
    hasAtLeastOneGluon = merge2Qto1G(cand, gluons, quarks, protectV);
    // Recreate the intermediate Vs of the cand object, and re-assign protectV (although unchanged, the pointer is now invalid).
    if (hasAtLeastOneGluon){
      cand->recreateVs();
      protectV=getProtectedV(cand);
      //TUtil::PrintCandidateSummary(cand);
      //MELAout << "prptectV = " << (protectV ? protectV->id : -9000) << endl;
      if (protectVStrict && protectV){
        for (int iqq=quarks.size()-1; iqq>=0; iqq--){
          MELAParticle*& quark = quarks.at(iqq);
          if (MELAParticle::checkParticleExists(quark, protectV->getDaughters())) quarks.erase(quarks.begin()+iqq);
        }
      }
    }
    else{
      MELAerr << "Failed to merge gluons! Candidate summary:" << endl;
      TUtil::PrintCandidateSummary(cand);
      exit(1);
    }
  }

  // MERGE
  vector<int> bestMerge;
  if (hasProtectedV) protectedVs.push_back(protectV);
  // Merge the gluons into quarks
  if (getMergeOrder_GluonsIntoQuarks(gluons, quarks, bestMerge, true, &protectedVs)>=0.){
    for (unsigned int ig=0; ig<bestMerge.size(); ig++){
      unsigned int qindex, gindex;
      gindex=ig;
      qindex=bestMerge.at(ig);

      MELAParticle* theGluon = gluons.at(gindex);
      MELAParticle* theQuark = quarks.at(qindex);

      // Add gluon momentum to the chosen quark
      //MELAout << "Considering to merge gluon " << gindex << " (E=" << theGluon->t() << ", stat=" << theGluon->genStatus << ") with quark " << qindex << " (E=" << theQuark->t() << ", id=" << theQuark->id << ", stat=" << theQuark->genStatus << ")" << endl;
      if (theQuark->genStatus*theGluon->genStatus>0) (*theQuark) += theGluon->p4;
      else (*theQuark) += -theGluon->p4;
      //MELAout << "After merging: Quark " << qindex << " (E=" << theQuark->t() << ", id=" << theQuark->id << ", stat=" << theQuark->genStatus << ")" << endl;
    }
  }
  // Re-adjust V daughter momenta/ids/genStatus
  for (auto& protV:protectedVs){
    protV->p4 = TLorentzVector(0, 0, 0, 0); // Reset protected V momentum
    for (auto& Vdau:protV->getDaughters()) protV->p4 = protV->p4 + Vdau->p4;
  }
  // END MERGE

  // Reorder all particles
  int swapconfig;
  vector<int> qordered;
  vector<MELAParticle*> dummyVectorPart;
  vector<int> dummyVectorInt;
  vector<pair<MELAParticle*, vector<MELAParticle*>>> quark_gluon_collection = mapGluonsToQuarks(dummyVectorPart, quarks, dummyVectorInt);
  /*
  for (auto const& qgpair:quark_gluon_collection){
    MELAout << *(qgpair.first) << " : ";
    for (MELAParticle* p:qgpair.second) MELAout << *p << " | ";
    MELAout << endl;
  }
  */
  /*double cfglikelihood=-1;*/
  if (candScheme==TVar::JJVBF) /*cfglikelihood=*/getBestVBFConfig(quark_gluon_collection, &qordered, &swapconfig);
  else if (hasProtectedV) /*cfglikelihood=*/getBestVHConfig(protectV, quark_gluon_collection, &qordered, &swapconfig);
  //MELAout << "qordered: " << qordered << endl;
  if (!(qordered.empty() || qordered.size()==quarks.size())){
    MELAerr << "qordered.empty() ? " << qordered.empty() << endl;
    MELAerr << "qordered.size()=?quarks.size() " << qordered.size() << " ?= " << quarks.size() << endl;
    TUtil::PrintCandidateSummary(cand);
  }
  /*
  else{
  MELAout << "qordered: ";
  for (auto& ord:qordered) MELAout << quarks.at(ord)->id << " (" << ord << " , " << MELAParticle::checkParticleExists(quarks.at(ord), protectV->getDaughters()) << ") ";
  MELAout << endl;
  }
  */
  assert(qordered.empty() || qordered.size()==quarks.size());
  for (unsigned int qord=0; qord<qordered.size(); qord++){
    const int& qindex = qordered.at(qord);
    MELAParticle*& theQuark = quarks.at(qindex);
    //MELAout << "Starting quark = " << *theQuark << endl;
    //MELAout << "Starting quark status = " << theQuark->genStatus << endl;
    //MELAout << "Swap ? " << (swapconfig%2==1) << endl;
    if (swapconfig%2==1){
      theQuark->p4 = -theQuark->p4;
      theQuark->genStatus *= -1;
      theQuark->id *= -1;
    }
    swapconfig = swapconfig >> 1;
  }

  // Get the modified candidate
  SimpleParticleCollection_t associatednew;
  {
    SimpleParticleCollection_t daughtersnew;
    SimpleParticleCollection_t mothersnew;
    MELACandidate* candtmp;
    copyCandidate(cand, candtmp, true);
    readCandidate(candtmp, mothersnew, daughtersnew, associatednew);
    if (mothersnew.size()!=2){
      TUtil::PrintCandidateSummary(cand);
      TUtil::PrintCandidateSummary(candtmp);
      exit(1);
    }
    delete cand; cand=candtmp;
  }

  // Iterate if the number of associated particles is still greater than the required amount
  if (associatednew.size()>nAssociated) reduceJJtoQuarks(cand);
}

void MELACandidateRecaster::deduceLOVHTopology(MELACandidate*& cand){
  if (candScheme==TVar::nProductions) return;
  if (!cand) return;
  bool protectVStricttmp=protectVStrict;
  protectVStrict=true; // No need to modify the V in the intermediate candidate
  MELACandidate* candLookUp;
  copyCandidate(cand, candLookUp, false);
  reduceJJtoQuarks(candLookUp);
  protectVStrict=protectVStricttmp;
  MELAParticle* protectVLookUp = getProtectedV(candLookUp);
  const bool hasProtectedV = (protectVLookUp!=nullptr);
  if (!hasProtectedV) MELAerr << "MELACandidateRecaster::deduceLOVHTopology ERROR: NO PROTECTED V!" << endl;

  // Find the protected V manually
  MELAParticle* protectV = nullptr;
  {
    const double mLookUp = protectVLookUp->m();
    for (auto& part : cand->getSortedVs()){
      const double mMatch = part->m();
      if (fabs(mMatch-mLookUp)<1e-6 && protectVLookUp->id==part->id){
        protectV=part;
        break;
      }
    }
    if (!protectV) MELAerr << "MELACandidateRecaster::deduceLOVHTopology ERROR: NO MATCHED PROTECTED V!" << endl;
  }

  vector<pair<int, double>> mother_idpz;
  for (int im=0; im<candLookUp->getNMothers(); im++){
    mother_idpz.push_back(pair<int, double>(candLookUp->getMother(im)->id, candLookUp->getMother(im)->z()));
  }
  MELAParticle* mot[2] ={
    cand->getMother(0),
    cand->getMother(1)
  };
  // Count the unmodified gluons
  const bool maskMotherIdWhenNoGluons=false;
  unsigned int nGs=0;
  if (maskMotherIdWhenNoGluons){
    for (int ip=0; ip<cand->getNMothers(); ip++){ MELAParticle* part = cand->getMother(ip); if (PDGHelpers::isAGluon(part->id)) nGs++; }
    for (auto& part:cand->getAssociatedJets()){ if (PDGHelpers::isAGluon(part->id)) nGs++; }
  }

  if (mother_idpz.size()==2 && mot[0] && mot[1]){
    if (TMath::Sign(1., mot[0]->z()-mot[1]->z())==TMath::Sign(1., mother_idpz[0].second - mother_idpz[1].second)){
      mot[0]->id = mother_idpz[0].first;
      mot[1]->id = mother_idpz[1].first;
    }
    else{
      mot[1]->id = mother_idpz[0].first;
      mot[0]->id = mother_idpz[1].first;
    }
    // Sum over associated particle and daughter momenta so that mother momneta are boosted to the correct system
    // Also disable any other irrelevant particle
    TLorentzVector pAD = cand->p4 + protectV->p4;
    for (auto& part:cand->getAssociatedJets()){
      if (PDGHelpers::isAGluon(part->id) || !MELAParticle::checkParticleExists(part, protectV->getDaughters())) part->setSelected(false);
    }
    TLorentzVector pTotalVH = pAD;
    pTotalVH.Boost(-TVector3(pTotalVH.X()/pTotalVH.T(), pTotalVH.Y()/pTotalVH.T(), 0));
    double sysE = pTotalVH.T();
    double sysZ = pTotalVH.Z();
    double Ep = (sysE+sysZ)/2.;
    double Em = (sysE-sysZ)/2.;
    double pp = Ep;
    double pm = -Em;
    if (TMath::Sign(1., mot[0]->z()-mot[1]->z())==TMath::Sign(1., pp-pm)){
      mot[0]->p4 = TLorentzVector(0, 0, pp, Ep);
      mot[1]->p4 = TLorentzVector(0, 0, pm, Em);
    }
    else{
      mot[1]->p4 = TLorentzVector(0, 0, pp, Ep);
      mot[0]->p4 = TLorentzVector(0, 0, pm, Em);
    }
    TLorentzVector pMot = mot[0]->p4 + mot[1]->p4;
    for (unsigned int im=0; im<2; im++){
      mot[im]->p4.Boost(-pMot.BoostVector());
      mot[im]->p4.Boost(pAD.BoostVector());
    }
    // Lastly, set mother id to 0 if there are no gluons (ambiguity in mother ids)
    if (maskMotherIdWhenNoGluons && nGs==0){ for (unsigned int im=0; im<2; im++) mot[im]->id=0; }
  }
  delete candLookUp;
  {
    MELACandidate* candtmp;
    copyCandidate(cand, candtmp, false);
    delete cand; cand=candtmp;
  }
}

void MELACandidateRecaster::deduceLOVBFTopology(MELACandidate*& cand){ reduceJJtoQuarks(cand); }

void MELACandidateRecaster::deduceLOHJJTopology_JetMerge(MELACandidate*& cand){
  if (candScheme!=TVar::JJQCD) return;
  if (!cand) return;

  // Gather all quarks and gluons
  unsigned int nQ = 0;
  unsigned int nG = 0;
  vector<MELAParticle*> partons; // Only need a single collection
  // Gather all mother quarks and gluons
  for (int ip=0; ip<cand->getNMothers(); ip++){
    MELAParticle* part = cand->getMother(ip);
    part->p4 = -part->p4;
    if (PDGHelpers::isAGluon(part->id)) nG++;
    else nQ++;
    partons.push_back(part);
  }
  // Gather all associated quarks and gluons except the protected V
  for (auto& part:cand->getAssociatedJets()){
    if (PDGHelpers::isAGluon(part->id)) nG++;
    else nQ++;
    partons.push_back(part);
  }

  int nExtraJets = nQ+nG-4;
  if (nExtraJets<0 || nQ%2==1){
    MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: (nQ, nG) = (" << nQ << ", " << nG << ") is not valid!"  << endl;
    exit(1);
  }
  else if (nExtraJets==0) return;
  else if (nExtraJets>1){
    MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: (nQ, nG) = (" << nQ << ", " << nG << ") is valid, but the implementation only supports 1 extra jet at this moment."  << endl;
    exit(1);
  }

  vector<pair<MELAParticle*, vector<MELAParticle*>>> bestMergeConfig; bestMergeConfig.reserve(4);
  double bestScore=-99;

  vector<vector<int>> partonPermutations;
  CombinationGenerator(partons.size(), 4, partonPermutations, 0);
  for (auto& perm:partonPermutations) std::sort(perm.begin(), perm.end()); // Sort the permutation index list such that a mother always come first
  for (auto& perm:partonPermutations){
    vector<MELAParticle*> missingPartons; missingPartons.reserve(nExtraJets);
    for (unsigned int i=0; i<nQ+nG; i++){ if (std::find(perm.begin(), perm.end(), static_cast<int>(i))==perm.end()) missingPartons.push_back(partons.at(i)); }
    if ((int) missingPartons.size()!=nExtraJets){
      MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: missingPartons.size() = " << missingPartons.size() << " != nExtraJets = " << nExtraJets << endl;
      exit(1);
    }

    // Check if sum of charges of extra partons is 0.
    {
      float sumChargeQExtra=0;
      for (auto const& part:missingPartons) sumChargeQExtra += part->charge();
      if (fabs(sumChargeQExtra)>1e-5) continue;
    }

    if (nExtraJets==1){
      double score=-1; // Score of merge
      int iSelected=-1; // Index of parton to merge into
      for (auto const& ipart:perm){
        double tmpscore = getMergeScore(partons.at(ipart), missingPartons.front());
        if (tmpscore<0.) continue;
        if (score<0. || tmpscore<score){
          iSelected=ipart;
          score=tmpscore;
        }
      }
      if (iSelected>=0 && (bestScore<0. || bestScore>score)){
        bestScore = score;
        bestMergeConfig.clear();
        for (auto const& ipart:perm){
          vector<MELAParticle*> tmpMergeList; tmpMergeList.reserve(nExtraJets);
          if (ipart == iSelected) tmpMergeList.push_back(missingPartons.front());
          bestMergeConfig.emplace_back(partons.at(ipart), tmpMergeList);
        }
      }
    }
    else{
      MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: nExtraJets = " << nExtraJets << " is not implemented." << endl;
      exit(1);
    }
  }

  // MERGE
  if (!bestMergeConfig.empty()){
    for (auto& cfg:bestMergeConfig){
      MELAParticle*& targetPart = cfg.first;
      for (MELAParticle*& mergePart:cfg.second){
        *targetPart += mergePart->p4;
        mergePart->setSelected(false);
      }
      cfg.second.clear(); // Clear the list of particles to merge so that the same bestMergeConfig collection can be reused below
    }
  }
  // END MERGE
  
  // Revert the momenta of incoming particles
  for (MELAParticle*& part:partons){ if (part->genStatus<0) part->p4 = -part->p4; }

  // Reorder all particles
  int swapconfig;
  vector<int> part_order;
  getBestHJJConfig(bestMergeConfig, &part_order, &swapconfig);
  if (!(part_order.empty() || part_order.size()==bestMergeConfig.size())){
    MELAerr << "part_order.empty() ? " << part_order.empty() << endl;
    MELAerr << "part_order.size()=?bestMergeConfig.size() " << part_order.size() << " ?= " << bestMergeConfig.size() << endl;
    TUtil::PrintCandidateSummary(cand);
  }
  assert(part_order.empty() || part_order.size()==bestMergeConfig.size());
  for (unsigned int iord=0; iord<part_order.size(); iord++){
    const int& index = part_order.at(iord);
    MELAParticle*& thePart = bestMergeConfig.at(index).first;
    if (swapconfig%2==1){
      thePart->p4 = -thePart->p4;
      thePart->genStatus *= -1;
      thePart->id *= -1;
    }
    swapconfig = swapconfig >> 1;
  }

  // Get the modified candidate
  SimpleParticleCollection_t associatednew;
  {
    SimpleParticleCollection_t daughtersnew;
    SimpleParticleCollection_t mothersnew;
    MELACandidate* candtmp;
    copyCandidate(cand, candtmp, true);
    readCandidate(candtmp, mothersnew, daughtersnew, associatednew);
    if (mothersnew.size()!=2){
      TUtil::PrintCandidateSummary(cand);
      TUtil::PrintCandidateSummary(candtmp);
      exit(1);
    }
    delete cand; cand=candtmp;
  }
}

void MELACandidateRecaster::deduceLOHJJTopology_LeadingPt(MELACandidate*& cand){
  if (candScheme!=TVar::JJQCD) return;
  if (!cand) return;

  // Gather all quarks and gluons
  unsigned int nQ = 0;
  unsigned int nG = 0;
  vector<MELAParticle*> partons; // Only need a single collection
                                 // Gather all mother quarks and gluons
  for (int ip=0; ip<cand->getNMothers(); ip++){
    MELAParticle* part = cand->getMother(ip);
    if (PDGHelpers::isAGluon(part->id)) nG++;
    else nQ++;
    partons.push_back(part);
  }
  // Gather all associated quarks and gluons except the protected V
  for (auto& part:cand->getAssociatedJets()){
    if (PDGHelpers::isAGluon(part->id)) nG++;
    else nQ++;
    partons.push_back(part);
  }

  int nExtraJets = nQ+nG-4;
  if (nExtraJets<0 || nQ%2==1){
    MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: (nQ, nG) = (" << nQ << ", " << nG << ") is not valid!"  << endl;
    exit(1);
  }
  else if (nExtraJets==0) return;
  else if (nExtraJets>1){
    MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: (nQ, nG) = (" << nQ << ", " << nG << ") is valid, but the implementation only supports 1 extra jet at this moment."  << endl;
    exit(1);
  }

  vector<pair<MELAParticle*, vector<MELAParticle*>>> bestMergeConfig; bestMergeConfig.reserve(4);
  vector<MELAParticle*> bestMissingPartons; bestMissingPartons.reserve(nExtraJets);
  double bestScore=-99; // Picks smallest score

  vector<vector<int>> partonPermutations;
  CombinationGenerator(partons.size(), 4, partonPermutations, 0);
  for (auto& perm:partonPermutations) std::sort(perm.begin(), perm.end()); // Sort the permutation index list such that a mother always come first
  for (auto& perm:partonPermutations){
    if (!(perm.at(0)==0 && perm.at(1)==1)) continue;

    vector<MELAParticle*> missingPartons; missingPartons.reserve(nExtraJets);
    for (unsigned int i=0; i<nQ+nG; i++){ if (std::find(perm.begin(), perm.end(), static_cast<int>(i))==perm.end()) missingPartons.push_back(partons.at(i)); }
    if ((int) missingPartons.size()!=nExtraJets){
      MELAerr << "MELACandidateRecaster::deduceLOHJJTopology: missingPartons.size() = " << missingPartons.size() << " != nExtraJets = " << nExtraJets << endl;
      exit(1);
    }

    // Check if sum of charges of extra partons is 0.
    {
      float sumChargeQExtra=0;
      for (auto const& part:missingPartons) sumChargeQExtra += part->charge();
      if (fabs(sumChargeQExtra)>1e-5) continue;
    }

    double score=0; // Score of merge
    for (auto const& part:missingPartons) score += part->pt();
    if (bestScore<0. || bestScore>score){
      bestScore = score;
      bestMergeConfig.clear(); bestMissingPartons.clear();
      vector<MELAParticle*> tmplist;
      for (auto const& ipart:perm) bestMergeConfig.emplace_back(partons.at(ipart), tmplist);
      for (MELAParticle* mpart:missingPartons) bestMissingPartons.push_back(mpart);
    }
  }

  for (auto* mpart:bestMissingPartons) mpart->setSelected(false);

  // Reorder all particles
  int swapconfig;
  vector<int> part_order;
  getBestHJJConfig(bestMergeConfig, &part_order, &swapconfig);
  if (!(part_order.empty() || part_order.size()==bestMergeConfig.size())){
    MELAerr << "part_order.empty() ? " << part_order.empty() << endl;
    MELAerr << "part_order.size()=?bestMergeConfig.size() " << part_order.size() << " ?= " << bestMergeConfig.size() << endl;
    TUtil::PrintCandidateSummary(cand);
  }
  assert(part_order.empty() || part_order.size()==bestMergeConfig.size());
  for (unsigned int iord=0; iord<part_order.size(); iord++){
    const int& index = part_order.at(iord);
    MELAParticle*& thePart = bestMergeConfig.at(index).first;
    if (swapconfig%2==1){
      thePart->p4 = -thePart->p4;
      thePart->genStatus *= -1;
      thePart->id *= -1;
    }
    swapconfig = swapconfig >> 1;
  }

  // Get the modified candidate
  SimpleParticleCollection_t associatednew;
  {
    SimpleParticleCollection_t daughtersnew;
    SimpleParticleCollection_t mothersnew;
    MELACandidate* candtmp;
    copyCandidate(cand, candtmp, true);
    readCandidate(candtmp, mothersnew, daughtersnew, associatednew);
    if (mothersnew.size()!=2){
      TUtil::PrintCandidateSummary(cand);
      TUtil::PrintCandidateSummary(candtmp);
      exit(1);
    }
    delete cand; cand=candtmp;
  }
}
