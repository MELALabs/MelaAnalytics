#include <utility>
#include <algorithm>
#include "MELAComputation.h"

using namespace std;


MELAComputation::MELAComputation(MELAHypothesis* targetP_) :
  targetP(targetP_)
{
  setOption(targetP->getOption());
  reset();
}
MELAComputation::~MELAComputation(){
  addedP.clear();
  subtractedP.clear();
  multipliedP.clear();
  dividedP.clear();
  maximize_num.clear();
  maximize_denom.clear();
}

void MELAComputation::addContingencies(vector<MELAHypothesis*>& allHypos){
  addContingency(allHypos, opt->addedAliases, addedP/*, 0*/);
  addContingency(allHypos, opt->subtractedAliases, subtractedP/*, 0*/);
  addContingency(allHypos, opt->multipliedAliases, multipliedP/*, 0*/);
  addContingency(allHypos, opt->dividedAliases, dividedP/*, 0*/);

  addContingency(allHypos, opt->maximizationNumerators, maximize_num, 1);
  addContingency(allHypos, opt->maximizationDenominators, maximize_denom, 1);
}
void MELAComputation::addContingency(vector<MELAHypothesis*>& allHypos, vector<string>& source, vector<MELAHypothesis*>& dest, unsigned int setHypoFlag){
  for (string const& aliasToMatch:source){
    // Match by alias, not by name!
    for (MELAHypothesis*& hypo:allHypos){
      if (hypo->getOption()->isAliased() && aliasToMatch==hypo->getOption()->getAlias()){
        MELAHypothesis* matchedHypo = hypo;
        if (setHypoFlag==1) matchedHypo->setMaximizationClientStatus(true);
        dest.push_back(matchedHypo);
        break;
      }
    }

  }
}

Bool_t MELAComputation::testMaximizationCache(){
  Bool_t updateMany = ((maximize_num.size()+maximize_denom.size())>0);
  Float_t testCache = (updateMany ? 1. : -1.);
  // Always use type UseME for such comparisons, the others don't make much sense
  for (MELAHypothesis const* num:maximize_num) testCache *= num->getVal(MELAHypothesis::UseME);
  for (MELAHypothesis const* denom:maximize_denom){
    Float_t divVal = denom->getVal(MELAHypothesis::UseME);
    if (divVal!=0.) testCache /= divVal;
  }
  if (testCache>=maximizationCachedVal){
    if (updateMany) maximizationCachedVal = testCache;
    else maximizationCachedVal = 0.;
    return true;
  }
  else return false;
}
void MELAComputation::update(){
  if (contUpdate && testMaximizationCache()){
    pME = extractVal(MELAHypothesis::UseME);
    pAux = extractVal(MELAHypothesis::UsePAux);
    cMEAvg = extractVal(MELAHypothesis::UsePConstant);
    contUpdate = false;
  }
}
void MELAComputation::forceUpdate(){
  if (testMaximizationCache()){
    pME = extractVal(MELAHypothesis::UseME);
    pAux = extractVal(MELAHypothesis::UsePAux);
    cMEAvg = extractVal(MELAHypothesis::UsePConstant);
  }
}

Float_t MELAComputation::extractVal(MELAHypothesis::METype valtype){
  Float_t tmp = targetP->getVal(valtype);
  for (MELAHypothesis const* hypo:addedP){ if (hypo) tmp += hypo->getVal(valtype); }
  for (MELAHypothesis const* hypo:subtractedP){ if (hypo) tmp -= hypo->getVal(valtype); }

  Float_t factor=1;
  for (MELAHypothesis const* hypo:multipliedP){ if (hypo) factor *= hypo->getVal(valtype); }
  for (MELAHypothesis const* hypo:dividedP){
    if (hypo){
      Float_t divVal = hypo->getVal(valtype);
      if (divVal!=0.) factor /= divVal;
    }
  }

  if (!isnan(factor)) tmp *= factor;
  else tmp=0;

  return tmp;
}
Float_t MELAComputation::getVal(MELAHypothesis::METype valtype) const{
  switch (valtype){
  case MELAHypothesis::UseME:
    return pME;
  case MELAHypothesis::UsePAux:
    return pAux;
  case MELAHypothesis::UsePConstant:
    return cMEAvg;
  default:
    return 0;
  }
}

void MELAComputation::resetMaximizationCache(){
  contUpdate=true;
  maximizationCachedVal=-1.;
}
void MELAComputation::reset(){
  resetMaximizationCache();
  targetP->reset();
  pME = targetP->getVal(MELAHypothesis::UseME);
  pAux = targetP->getVal(MELAHypothesis::UsePAux);
  cMEAvg = targetP->getVal(MELAHypothesis::UsePConstant);
}

void MELAComputation::Print() const{
  cout << "**********" << endl;
  cout << "MELAComputation summary:" << endl;
  cout << "\tPrimary hypothesis name: " << opt->getName() << endl;
  cout << "\tCluster: " << opt->getCluster() << endl;
  cout << "\tAlias: " << opt->getAlias() << endl;
  if (!opt->isCopy()) cout << "\tCopy alias: " << opt->getCopyAlias() << endl;
  
  cout << "\tComputes pm4l?: " << (opt->usePM4L() ? "True" : "False") << endl;
  cout << "\tComputes associated V pmjj?: " << (opt->usePMaVJJ() ? "True" : "False") << endl;
  cout << "\tCan update? " << (contUpdate ? "True" : "False") << endl;

  cout << "\tValue formula: ";
  cout << "(Self";
  for (MELAHypothesis const* hypo:addedP){
    cout << " + ";
    cout << hypo->getOption()->getAlias();
  }
  for (MELAHypothesis const* hypo:subtractedP){
    cout << " - ";
    cout << hypo->getOption()->getAlias();
  }
  cout << ") * (";
  for (unsigned int ip=0; ip<multipliedP.size(); ip++){
    if (ip>0) cout << " * ";
    cout << multipliedP.at(ip)->getOption()->getAlias();
  }
  cout << ") / (";
  for (unsigned int ip=0; ip<dividedP.size(); ip++){
    if (ip>0) cout << " * ";
    cout << dividedP.at(ip)->getOption()->getAlias();
  }
  cout << ")" << endl;
  cout << "\tCurrent (ME, pAux, pConst) = (" << pME << ", " << pAux << ", " << cMEAvg << ")" << endl;

  cout << "\tMaximized formula: ";
  cout << "(";
  for (unsigned int ip=0; ip<maximize_num.size(); ip++){
    if (ip>0) cout << " * ";
    cout << maximize_num.at(ip)->getOption()->getAlias();
  }
  cout << ") / (";
  for (unsigned int ip=0; ip<maximize_denom.size(); ip++){
    if (ip>0) cout << " * ";
    cout << maximize_denom.at(ip)->getOption()->getAlias();
  }
  cout << ")" << endl;
  cout << "\tCached maximized value = " << maximizationCachedVal << endl;
  cout << "**********" << endl;
}


