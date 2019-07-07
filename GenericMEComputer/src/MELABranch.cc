#include "MELABranch.h"


using namespace BranchHelpers;
using namespace std;


MELABranch::MELABranch(
  TTree* theTree_, TString bname_, Float_t defVal_,
  MELAComputation* computer_
  ) :
  ExtendedBranch(theTree_, bname_, defVal_, false),
  computer(computer_)
{
  // Decide what ME to compute by the branch name because AddPAux alone does not decide what the branch does.
  // It just duplicates the branch with the same computer.
  if (bname.Contains("pAux")) valtype = MELAHypothesis::UsePAux;
  else if (bname.Contains("pConst")) valtype = MELAHypothesis::UsePConstant;
  else valtype = MELAHypothesis::UseME;
}

void MELABranch::setVal(){
  if (computer) ExtendedBranch::setValue(computer->getVal(valtype));
}

void MELABranch::Print(){
  ExtendedBranch::Print();
  cout << "\tMELABranch extras:" << endl;
  cout << "\tME type=" << valtype << endl;
  if (computer) computer->Print();
  cout << "**********" << endl;
}
