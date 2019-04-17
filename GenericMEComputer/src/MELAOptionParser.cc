#include "MELAOptionParser.h"

using namespace std;


MELAOptionParser::MELAOptionParser(string stropts) :
strCluster("Common"),
propScheme(TVar::FixedWidth),
noBranching(false),
includePAux(false),
includePConst(false),
isPM4L(false),
isPMaVJJ(false),
isPMaVJJTrue(false),
isProp(false),
isGenProb(false),
defME(0.),
hmass(-99),
h2mass(-99),
hwidth(-99),
h2width(-99)
{
  // Split all options by whitespace
  splitOptionRecursive(stropts, rawOptions, ' ');
  analyze();
}
void MELAOptionParser::analyze(){ analyze(rawOptions); }
void MELAOptionParser::analyze(const std::vector<std::string>& optcoll){
  char rawdelimiter = ':';
  for (std::string const& opt : optcoll){
    string wish, value;
    splitOption(opt, wish, value, rawdelimiter);
    interpretOption(wish, value);
  }
  // Check options
  if (strName==""){ cerr << "MELAOptionParser::analyze: No name detected. Please put a name!" << endl; assert(0); }
  if (isPM4L && isPMaVJJ){ cerr << "MELAOptionParser::analyze: Cannot be defined as both P(m4l) and P(mjj)! Choose only one" << endl; assert(0); }
  if (strAlias=="<Name>") strAlias=strName;
  if (isPMaVJJTrue) isPMaVJJ = isPMaVJJTrue;
  if (isCopy()){
    if (DEBUG_MB){
      cout
        << "MELAOptionParser::analyze: Branch " << strName
        << " will be a copy of the hypothesis with alias " << strCopyAlias
        << ". ME properties will be overridden with the ME options of the original hypothesis."
        << endl;
    }
    if (isAliased()){
      cerr << "MELAOptionParser::analyze: Branch " << strName << " cannot any aliases." << endl;
      strAlias = "";
    }
  }
}
void MELAOptionParser::splitOption(const string rawoption, string& wish, string& value, char delimiter)const{
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}
void MELAOptionParser::splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="" && !checkListVariable(splitoptions, result)) splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}
Bool_t MELAOptionParser::checkListVariable(const vector<string>& list, const string& var)const{
  for (auto const v:list){
    if (v==var) return true; // Look for exact match
  }
  return false;
}

void MELAOptionParser::interpretOption(string wish, string value){
  if (wish.empty()){
    cerr << "MELAOptionParser::interpretOption: Unknown option with value " << value << endl;
  }
  else if (wish=="Process") setProcess(value);
  else if (wish=="Production") setProduction(value);
  else if (wish=="MatrixElement") setME(value);
  else if (wish=="SuperMelaSyst") setSuperMelaSyst(value);
  else if (wish=="PropScheme") setPropagatorScheme(value);

  else if (wish=="isGen") isGenProb = Bool_t(((UShort_t)atoi(value.c_str()))>0);
  else if (wish=="isPM4L") isPM4L = Bool_t(((UShort_t)atoi(value.c_str()))>0);
  else if (wish=="isPMaVJJ") isPMaVJJ = Bool_t(((UShort_t) atoi(value.c_str()))>0);
  else if (wish=="isPMaVJJTrue") isPMaVJJTrue = Bool_t(((UShort_t) atoi(value.c_str()))>0);
  else if (wish=="isProp") isProp = Bool_t(((UShort_t) atoi(value.c_str()))>0);
  else if (wish=="NoBranch") noBranching = Bool_t(((UShort_t)atoi(value.c_str()))>0);

  else if (wish=="DefaultME") defME = (Float_t)atof(value.c_str());
  else if (wish=="hmass" || wish=="MH") hmass = (Float_t)atof(value.c_str());
  else if (wish=="h2mass" || wish=="MH2") h2mass = (Float_t)atof(value.c_str());
  else if (wish=="hwidth" || wish=="GaH") hwidth = (Float_t)atof(value.c_str());
  else if (wish=="h2width" || wish=="GaH2") h2width = (Float_t)atof(value.c_str());

  else if (wish=="Name") strName = value;
  else if (wish=="Alias") strAlias = value;
  else if (wish=="Copy" || wish=="CopyFrom") strCopyAlias = value;
  else if (wish=="Cluster") strCluster = value;

  else if (wish=="Options"){
    vector<string> vtmp;
    splitOptionRecursive(value, vtmp, ';');
    for (string const& opt:vtmp){
      if (opt.find("AddPAux")!=string::npos){
        string tmpwish, tmpvalue;
        splitOption(opt, tmpwish, tmpvalue, '=');
        includePAux=Bool_t(((UShort_t)atoi(tmpvalue.c_str()))>0);
      }
      else if (opt.find("AddPConst")!=string::npos){
        string tmpwish, tmpvalue;
        splitOption(opt, tmpwish, tmpvalue, '=');
        includePConst=Bool_t(((UShort_t)atoi(tmpvalue.c_str()))>0);
      }
      else if (opt.find("AddP")!=string::npos) setAddedAliases(opt);
      else if (opt.find("SubtractP")!=string::npos) setSubtractedAliases(opt);
      else if (opt.find("MultiplyP")!=string::npos) setMultipliedAliases(opt);
      else if (opt.find("DivideP")!=string::npos) setDividedAliases(opt);
      else if (opt.find("MaxNumerator")!=string::npos) setMaximizationNumAliases(opt);
      else if (opt.find("MaxDenominator")!=string::npos) setMaximizationDenomAliases(opt);
    }
  }

  else if (wish=="Couplings"){
    vector<string> vtmp;
    splitOptionRecursive(value, vtmp, ';');
    for (string const& opt:vtmp) extractCoupling(opt);
  }

  else if (wish=="ForceIncomingFlavors"){
    vector<string> vtmp;
    splitOptionRecursive(value, vtmp, ';');
    forcedIncomingFlavorList.reserve(forcedIncomingFlavorList.size()+vtmp.size());
    for (string const& opt:vtmp){
      string strid1, strid2;
      splitOption(opt, strid1, strid2, ',');
      Int_t id1 = atoi(strid1.c_str());
      Int_t id2 = atoi(strid2.c_str());
      if (PDGHelpers::isAJet(abs(id1)) && PDGHelpers::isAJet(abs(id2))){
        for (int iq=-5; iq<=5; iq++){
          if (!(
            PDGHelpers::isAnUnknownJet(id1)
            ||
            (PDGHelpers::isAQuark(id1) && id1==iq)
            ||
            (PDGHelpers::isAGluon(id1) && id1>0 && iq==0)
            ||
            (PDGHelpers::isAGluon(-id1) && id1<0 && iq!=0) // Special case: -21 means "not-a-gluon"
            )) continue;
          for (int jq=-5; jq<=5; jq++){
            if (!(
              PDGHelpers::isAnUnknownJet(id2)
              ||
              (PDGHelpers::isAQuark(id2) && id2==jq)
              ||
              (PDGHelpers::isAGluon(id2) && id2>0 && jq==0)
              ||
              (PDGHelpers::isAGluon(-id2) && id2<0 && jq!=0) // Special case: -21 means "not-a-gluon"
              )) continue;
            forcedIncomingFlavorList.emplace_back(iq, jq); // BE CAREFUL: 0 MEANS GLUON HERE!
          }
        }
      }
      else cerr << "MELAOptionParser::interpretOption: Incoming id pair (" << id1 << "," << id2 << ") is not supported in option " << wish << "!" << endl;
    }
  }

  else cerr << "MELAOptionParser::interpretOption: Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void MELAOptionParser::setProcess(string wish){
  for (int iproc=TVar::HSMHiggs; iproc!=TVar::nProcesses; iproc++){
    if (wish == TVar::ProcessName((TVar::Process) iproc).Data()){ proc = (TVar::Process) iproc; return; }
  }
  cerr << "MELAOptionParser::setProcess(" << wish << "): Failed to find the proper process." << endl;
}
void MELAOptionParser::setProduction(string wish){
  if (wish=="GG"){ prod = TVar::ZZGG; return; }
  else if (wish=="QQB"){ prod = TVar::ZZQQB; return; }
  else if (wish=="QQB_STU"){ prod = TVar::ZZQQB_STU; return; }
  else if (wish=="QQB_S"){ prod = TVar::ZZQQB_S; return; }
  else if (wish=="QQB_TU"){ prod = TVar::ZZQQB_TU; return; }
  else if (wish=="INDEPENDENT"){ prod = TVar::ZZINDEPENDENT; return; }
  for (int iprod=TVar::ZZGG; iprod!=TVar::nProductions; iprod++){
    if (wish == TVar::ProductionName((TVar::Production) iprod).Data()){ prod = (TVar::Production) iprod; return; }
  }
  cerr << "MELAOptionParser::setProduction(" << wish << "): Failed to find the proper production." << endl;
}
void MELAOptionParser::setME(string wish){
  if (wish=="MCFM") ME = TVar::MCFM;
  else if (wish=="JHUGen") ME = TVar::JHUGen;
  else if (wish=="Analytical" || wish=="ANALYTICAL") ME = TVar::ANALYTICAL;
  else cerr << "MELAOptionParser::setME(" << wish << "): Failed to find the proper matrix element." << endl;
}
void MELAOptionParser::setSuperMelaSyst(string wish){
  if (wish=="SMSyst_None") superSyst = TVar::SMSyst_None;
  else if (wish=="SMSyst_ScaleUp") superSyst = TVar::SMSyst_ScaleUp;
  else if (wish=="SMSyst_ResUp") superSyst = TVar::SMSyst_ResUp;
  else if (wish=="SMSyst_ScaleDown") superSyst = TVar::SMSyst_ScaleDown;
  else if (wish=="SMSyst_ResDown") superSyst = TVar::SMSyst_ResDown;
  else cerr << "MELAOptionParser::setSuperMelaSyst(" << wish << "): Failed to find the proper SuperMELA systematics case." << endl;
}
void MELAOptionParser::setPropagatorScheme(std::string wish){
  if (wish=="NoPropagator") propScheme = TVar::NoPropagator;
  else if (wish=="RunningWidth") propScheme = TVar::RunningWidth;
  else if (wish=="FixedWidth") propScheme = TVar::FixedWidth;
  else if (wish=="CPS") propScheme = TVar::CPS;
  else cerr << "MELAOptionParser::setPropagatorScheme(" << wish << "): Failed to find the proper propagator case." << endl;
}
void MELAOptionParser::extractCoupling(string opt){
  string wish, strVal, strValRe, strValIm;
  // Use double precision for couplings
  Double_t valRe=0;
  Double_t valIm=0;
  splitOption(opt, wish, strVal, '=');
  // Lambda and cz/cw couplings have no imaginary components, so do not expect to parse them with ','.
  if (
    wish.find("Lambda")==string::npos
    &&
    wish.find("q1sq")==string::npos
    &&
    wish.find("q2sq")==string::npos
    &&
    wish.find("q12sq")==string::npos
    &&
    wish.find("M_Zprime")==string::npos
    &&
    wish.find("Ga_Zprime")==string::npos
    &&
    wish.find("M_Wprime")==string::npos
    &&
    wish.find("Ga_Wprime")==string::npos
    &&
    wish!="separateWWZZcouplings"
    ){
    splitOption(strVal, strValRe, strValIm, ',');
    valRe = atof(strValRe.c_str());
    valIm = atof(strValIm.c_str());
  }
  else valRe = atof(strVal.c_str());

  // Here we go again, zillions of couplings
  if (wish=="separateWWZZcouplings") coupl_H.allow_WWZZSeparation((bool)valRe);
  // Spin-0 couplings, first resonance
  else if (wish=="kappa"){ coupl_H.Hqqcoupl[gHIGGS_KAPPA][0]=valRe; coupl_H.Hqqcoupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa_tilde"){ coupl_H.Hqqcoupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.Hqqcoupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="kappa_top"){ coupl_H.Httcoupl[gHIGGS_KAPPA][0]=valRe; coupl_H.Httcoupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa_tilde_top"){ coupl_H.Httcoupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.Httcoupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="kappa_bot"){ coupl_H.Hbbcoupl[gHIGGS_KAPPA][0]=valRe; coupl_H.Hbbcoupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa_tilde_bot"){ coupl_H.Hbbcoupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.Hbbcoupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="ghg2"){ coupl_H.Hggcoupl[gHIGGS_GG_2][0]=valRe; coupl_H.Hggcoupl[gHIGGS_GG_2][1]=valIm; }
  else if (wish=="ghg3"){ coupl_H.Hggcoupl[gHIGGS_GG_3][0]=valRe; coupl_H.Hggcoupl[gHIGGS_GG_3][1]=valIm; }
  else if (wish=="ghg4"){ coupl_H.Hggcoupl[gHIGGS_GG_4][0]=valRe; coupl_H.Hggcoupl[gHIGGS_GG_4][1]=valIm; }
  else if (wish=="kappa_4gen_top"){ coupl_H.Ht4t4coupl[gHIGGS_KAPPA][0]=valRe; coupl_H.Ht4t4coupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa_tilde_4gen_top"){ coupl_H.Ht4t4coupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.Ht4t4coupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="kappa_4gen_bot"){ coupl_H.Hb4b4coupl[gHIGGS_KAPPA][0]=valRe; coupl_H.Hb4b4coupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa_tilde_4gen_bot"){ coupl_H.Hb4b4coupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.Hb4b4coupl[gHIGGS_KAPPA_TILDE][1]=valIm; }

  else if (wish=="cz_q1sq"){ coupl_H.HzzCLambda_qsq[cLambdaHIGGS_VV_QSQ1]=(int)valRe; }
  else if (wish=="cz_q2sq"){ coupl_H.HzzCLambda_qsq[cLambdaHIGGS_VV_QSQ2]=(int)valRe; }
  else if (wish=="cz_q12sq"){ coupl_H.HzzCLambda_qsq[cLambdaHIGGS_VV_QSQ12]=(int)valRe; }
  else if (wish=="Lambda_z11"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_z21"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_z31"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_z41"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_z12"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_z22"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_z32"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_z42"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_z10"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda_z20"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda_z30"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda_z40"){ coupl_H.HzzLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="ghz1"){ coupl_H.Hzzcoupl[gHIGGS_VV_1][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1][1]=valIm; }
  else if (wish=="ghz2"){ coupl_H.Hzzcoupl[gHIGGS_VV_2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2][1]=valIm; }
  else if (wish=="ghz3"){ coupl_H.Hzzcoupl[gHIGGS_VV_3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3][1]=valIm; }
  else if (wish=="ghz4"){ coupl_H.Hzzcoupl[gHIGGS_VV_4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4][1]=valIm; }
  else if (wish=="ghzgs1_prime2"){ coupl_H.Hzzcoupl[gHIGGS_ZA_1_PRIME2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_ZA_1_PRIME2][1]=valIm; }
  else if (wish=="ghzgs2"){ coupl_H.Hzzcoupl[gHIGGS_ZA_2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_ZA_2][1]=valIm; }
  else if (wish=="ghzgs3"){ coupl_H.Hzzcoupl[gHIGGS_ZA_3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_ZA_3][1]=valIm; }
  else if (wish=="ghzgs4"){ coupl_H.Hzzcoupl[gHIGGS_ZA_4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_ZA_4][1]=valIm; }
  else if (wish=="ghgsgs2"){ coupl_H.Hzzcoupl[gHIGGS_AA_2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_AA_2][1]=valIm; }
  else if (wish=="ghgsgs3"){ coupl_H.Hzzcoupl[gHIGGS_AA_3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_AA_3][1]=valIm; }
  else if (wish=="ghgsgs4"){ coupl_H.Hzzcoupl[gHIGGS_AA_4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_AA_4][1]=valIm; }
  else if (wish=="ghz1_prime"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME][1]=valIm; }
  else if (wish=="ghz1_prime2"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME2][1]=valIm; }
  else if (wish=="ghz1_prime3"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME3][1]=valIm; }
  else if (wish=="ghz1_prime4"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME4][1]=valIm; }
  else if (wish=="ghz1_prime5"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME5][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME5][1]=valIm; }
  else if (wish=="ghz2_prime"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME][1]=valIm; }
  else if (wish=="ghz2_prime2"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME2][1]=valIm; }
  else if (wish=="ghz2_prime3"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME3][1]=valIm; }
  else if (wish=="ghz2_prime4"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME4][1]=valIm; }
  else if (wish=="ghz2_prime5"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME5][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME5][1]=valIm; }
  else if (wish=="ghz3_prime"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME][1]=valIm; }
  else if (wish=="ghz3_prime2"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME2][1]=valIm; }
  else if (wish=="ghz3_prime3"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME3][1]=valIm; }
  else if (wish=="ghz3_prime4"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME4][1]=valIm; }
  else if (wish=="ghz3_prime5"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME5][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME5][1]=valIm; }
  else if (wish=="ghz4_prime"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME][1]=valIm; }
  else if (wish=="ghz4_prime2"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME2][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME2][1]=valIm; }
  else if (wish=="ghz4_prime3"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME3][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME3][1]=valIm; }
  else if (wish=="ghz4_prime4"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME4][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME4][1]=valIm; }
  else if (wish=="ghz4_prime5"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME5][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME5][1]=valIm; }
  else if (wish=="ghz1_prime6"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME6][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME6][1]=valIm; }
  else if (wish=="ghz1_prime7"){ coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME7][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_1_PRIME7][1]=valIm; }
  else if (wish=="ghz2_prime6"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME6][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME6][1]=valIm; }
  else if (wish=="ghz2_prime7"){ coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME7][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_2_PRIME7][1]=valIm; }
  else if (wish=="ghz3_prime6"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME6][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME6][1]=valIm; }
  else if (wish=="ghz3_prime7"){ coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME7][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_3_PRIME7][1]=valIm; }
  else if (wish=="ghz4_prime6"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME6][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME6][1]=valIm; }
  else if (wish=="ghz4_prime7"){ coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME7][0]=valRe; coupl_H.Hzzcoupl[gHIGGS_VV_4_PRIME7][1]=valIm; }
  else if (wish=="cw_q1sq"){ coupl_H.HwwCLambda_qsq[cLambdaHIGGS_VV_QSQ1]=(int)valRe; }
  else if (wish=="cw_q2sq"){ coupl_H.HwwCLambda_qsq[cLambdaHIGGS_VV_QSQ2]=(int)valRe; }
  else if (wish=="cw_q12sq"){ coupl_H.HwwCLambda_qsq[cLambdaHIGGS_VV_QSQ12]=(int)valRe; }
  else if (wish=="Lambda_w11"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_w21"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_w31"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_w41"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda_w12"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_w22"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_w32"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_w42"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda_w10"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda_w20"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda_w30"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda_w40"){ coupl_H.HwwLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="ghw1"){ coupl_H.Hwwcoupl[gHIGGS_VV_1][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1][1]=valIm; }
  else if (wish=="ghw2"){ coupl_H.Hwwcoupl[gHIGGS_VV_2][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2][1]=valIm; }
  else if (wish=="ghw3"){ coupl_H.Hwwcoupl[gHIGGS_VV_3][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3][1]=valIm; }
  else if (wish=="ghw4"){ coupl_H.Hwwcoupl[gHIGGS_VV_4][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4][1]=valIm; }
  else if (wish=="ghw1_prime"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME][1]=valIm; }
  else if (wish=="ghw1_prime2"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME2][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME2][1]=valIm; }
  else if (wish=="ghw1_prime3"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME3][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME3][1]=valIm; }
  else if (wish=="ghw1_prime4"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME4][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME4][1]=valIm; }
  else if (wish=="ghw1_prime5"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME5][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME5][1]=valIm; }
  else if (wish=="ghw2_prime"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME][1]=valIm; }
  else if (wish=="ghw2_prime2"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME2][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME2][1]=valIm; }
  else if (wish=="ghw2_prime3"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME3][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME3][1]=valIm; }
  else if (wish=="ghw2_prime4"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME4][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME4][1]=valIm; }
  else if (wish=="ghw2_prime5"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME5][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME5][1]=valIm; }
  else if (wish=="ghw3_prime"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME][1]=valIm; }
  else if (wish=="ghw3_prime2"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME2][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME2][1]=valIm; }
  else if (wish=="ghw3_prime3"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME3][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME3][1]=valIm; }
  else if (wish=="ghw3_prime4"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME4][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME4][1]=valIm; }
  else if (wish=="ghw3_prime5"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME5][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME5][1]=valIm; }
  else if (wish=="ghw4_prime"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME][1]=valIm; }
  else if (wish=="ghw4_prime2"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME2][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME2][1]=valIm; }
  else if (wish=="ghw4_prime3"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME3][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME3][1]=valIm; }
  else if (wish=="ghw4_prime4"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME4][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME4][1]=valIm; }
  else if (wish=="ghw4_prime5"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME5][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME5][1]=valIm; }
  else if (wish=="ghw1_prime6"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME6][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME6][1]=valIm; }
  else if (wish=="ghw1_prime7"){ coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME7][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_1_PRIME7][1]=valIm; }
  else if (wish=="ghw2_prime6"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME6][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME6][1]=valIm; }
  else if (wish=="ghw2_prime7"){ coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME7][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_2_PRIME7][1]=valIm; }
  else if (wish=="ghw3_prime6"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME6][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME6][1]=valIm; }
  else if (wish=="ghw3_prime7"){ coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME7][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_3_PRIME7][1]=valIm; }
  else if (wish=="ghw4_prime6"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME6][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME6][1]=valIm; }
  else if (wish=="ghw4_prime7"){ coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME7][0]=valRe; coupl_H.Hwwcoupl[gHIGGS_VV_4_PRIME7][1]=valIm; }

  // Spin-0 couplings, second resonance
  else if (wish=="ghg2_4gen"){ coupl_H.Hg4g4coupl[gHIGGS_GG_2][0]=valRe; coupl_H.Hg4g4coupl[gHIGGS_GG_2][1]=valIm; }
  else if (wish=="ghg3_4gen"){ coupl_H.Hg4g4coupl[gHIGGS_GG_3][0]=valRe; coupl_H.Hg4g4coupl[gHIGGS_GG_3][1]=valIm; }
  else if (wish=="ghg4_4gen"){ coupl_H.Hg4g4coupl[gHIGGS_GG_4][0]=valRe; coupl_H.Hg4g4coupl[gHIGGS_GG_4][1]=valIm; }
  else if (wish=="kappa2"){ coupl_H.H2qqcoupl[gHIGGS_KAPPA][0]=valRe; coupl_H.H2qqcoupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa2_tilde"){ coupl_H.H2qqcoupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.H2qqcoupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="kappa2_top"){ coupl_H.H2ttcoupl[gHIGGS_KAPPA][0]=valRe; coupl_H.H2ttcoupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa2_tilde_top"){ coupl_H.H2ttcoupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.H2ttcoupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="kappa2_bot"){ coupl_H.H2bbcoupl[gHIGGS_KAPPA][0]=valRe; coupl_H.H2bbcoupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa2_tilde_bot"){ coupl_H.H2bbcoupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.H2bbcoupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="gh2g2"){ coupl_H.H2ggcoupl[gHIGGS_GG_2][0]=valRe; coupl_H.H2ggcoupl[gHIGGS_GG_2][1]=valIm; }
  else if (wish=="gh2g3"){ coupl_H.H2ggcoupl[gHIGGS_GG_3][0]=valRe; coupl_H.H2ggcoupl[gHIGGS_GG_3][1]=valIm; }
  else if (wish=="gh2g4"){ coupl_H.H2ggcoupl[gHIGGS_GG_4][0]=valRe; coupl_H.H2ggcoupl[gHIGGS_GG_4][1]=valIm; }
  else if (wish=="kappa2_4gen_top"){ coupl_H.H2t4t4coupl[gHIGGS_KAPPA][0]=valRe; coupl_H.H2t4t4coupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa2_tilde_4gen_top"){ coupl_H.H2t4t4coupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.H2t4t4coupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="kappa2_4gen_bot"){ coupl_H.H2b4b4coupl[gHIGGS_KAPPA][0]=valRe; coupl_H.H2b4b4coupl[gHIGGS_KAPPA][1]=valIm; }
  else if (wish=="kappa2_tilde_4gen_bot"){ coupl_H.H2b4b4coupl[gHIGGS_KAPPA_TILDE][0]=valRe; coupl_H.H2b4b4coupl[gHIGGS_KAPPA_TILDE][1]=valIm; }
  else if (wish=="gh2g2_4gen"){ coupl_H.H2g4g4coupl[gHIGGS_GG_2][0]=valRe; coupl_H.H2g4g4coupl[gHIGGS_GG_2][1]=valIm; }
  else if (wish=="gh2g3_4gen"){ coupl_H.H2g4g4coupl[gHIGGS_GG_3][0]=valRe; coupl_H.H2g4g4coupl[gHIGGS_GG_3][1]=valIm; }
  else if (wish=="gh2g4_4gen"){ coupl_H.H2g4g4coupl[gHIGGS_GG_4][0]=valRe; coupl_H.H2g4g4coupl[gHIGGS_GG_4][1]=valIm; }

  else if (wish=="c2z_q1sq"){ coupl_H.H2zzCLambda_qsq[cLambdaHIGGS_VV_QSQ1]=(int)valRe; }
  else if (wish=="c2z_q2sq"){ coupl_H.H2zzCLambda_qsq[cLambdaHIGGS_VV_QSQ2]=(int)valRe; }
  else if (wish=="c2z_q12sq"){ coupl_H.H2zzCLambda_qsq[cLambdaHIGGS_VV_QSQ12]=(int)valRe; }
  else if (wish=="Lambda2_z11"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_z21"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_z31"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_z41"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_z12"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_z22"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_z32"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_z42"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_z10"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda2_z20"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda2_z30"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda2_z40"){ coupl_H.H2zzLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="gh2z1"){ coupl_H.H2zzcoupl[gHIGGS_VV_1][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1][1]=valIm; }
  else if (wish=="gh2z2"){ coupl_H.H2zzcoupl[gHIGGS_VV_2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2][1]=valIm; }
  else if (wish=="gh2z3"){ coupl_H.H2zzcoupl[gHIGGS_VV_3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3][1]=valIm; }
  else if (wish=="gh2z4"){ coupl_H.H2zzcoupl[gHIGGS_VV_4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4][1]=valIm; }
  else if (wish=="gh2zgs1_prime2"){ coupl_H.H2zzcoupl[gHIGGS_ZA_1_PRIME2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_ZA_1_PRIME2][1]=valIm; }
  else if (wish=="gh2zgs2"){ coupl_H.H2zzcoupl[gHIGGS_ZA_2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_ZA_2][1]=valIm; }
  else if (wish=="gh2zgs3"){ coupl_H.H2zzcoupl[gHIGGS_ZA_3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_ZA_3][1]=valIm; }
  else if (wish=="gh2zgs4"){ coupl_H.H2zzcoupl[gHIGGS_ZA_4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_ZA_4][1]=valIm; }
  else if (wish=="gh2gsgs2"){ coupl_H.H2zzcoupl[gHIGGS_AA_2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_AA_2][1]=valIm; }
  else if (wish=="gh2gsgs3"){ coupl_H.H2zzcoupl[gHIGGS_AA_3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_AA_3][1]=valIm; }
  else if (wish=="gh2gsgs4"){ coupl_H.H2zzcoupl[gHIGGS_AA_4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_AA_4][1]=valIm; }
  else if (wish=="gh2z1_prime"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME][1]=valIm; }
  else if (wish=="gh2z1_prime2"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME2][1]=valIm; }
  else if (wish=="gh2z1_prime3"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME3][1]=valIm; }
  else if (wish=="gh2z1_prime4"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME4][1]=valIm; }
  else if (wish=="gh2z1_prime5"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME5][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME5][1]=valIm; }
  else if (wish=="gh2z2_prime"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME][1]=valIm; }
  else if (wish=="gh2z2_prime2"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME2][1]=valIm; }
  else if (wish=="gh2z2_prime3"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME3][1]=valIm; }
  else if (wish=="gh2z2_prime4"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME4][1]=valIm; }
  else if (wish=="gh2z2_prime5"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME5][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME5][1]=valIm; }
  else if (wish=="gh2z3_prime"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME][1]=valIm; }
  else if (wish=="gh2z3_prime2"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME2][1]=valIm; }
  else if (wish=="gh2z3_prime3"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME3][1]=valIm; }
  else if (wish=="gh2z3_prime4"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME4][1]=valIm; }
  else if (wish=="gh2z3_prime5"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME5][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME5][1]=valIm; }
  else if (wish=="gh2z4_prime"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME][1]=valIm; }
  else if (wish=="gh2z4_prime2"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME2][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME2][1]=valIm; }
  else if (wish=="gh2z4_prime3"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME3][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME3][1]=valIm; }
  else if (wish=="gh2z4_prime4"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME4][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME4][1]=valIm; }
  else if (wish=="gh2z4_prime5"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME5][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME5][1]=valIm; }
  else if (wish=="gh2z1_prime6"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME6][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME6][1]=valIm; }
  else if (wish=="gh2z1_prime7"){ coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME7][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_1_PRIME7][1]=valIm; }
  else if (wish=="gh2z2_prime6"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME6][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME6][1]=valIm; }
  else if (wish=="gh2z2_prime7"){ coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME7][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_2_PRIME7][1]=valIm; }
  else if (wish=="gh2z3_prime6"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME6][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME6][1]=valIm; }
  else if (wish=="gh2z3_prime7"){ coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME7][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_3_PRIME7][1]=valIm; }
  else if (wish=="gh2z4_prime6"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME6][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME6][1]=valIm; }
  else if (wish=="gh2z4_prime7"){ coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME7][0]=valRe; coupl_H.H2zzcoupl[gHIGGS_VV_4_PRIME7][1]=valIm; }
  else if (wish=="c2w_q1sq"){ coupl_H.H2wwCLambda_qsq[cLambdaHIGGS_VV_QSQ1]=(int)valRe; }
  else if (wish=="c2w_q2sq"){ coupl_H.H2wwCLambda_qsq[cLambdaHIGGS_VV_QSQ2]=(int)valRe; }
  else if (wish=="c2w_q12sq"){ coupl_H.H2wwCLambda_qsq[cLambdaHIGGS_VV_QSQ12]=(int)valRe; }
  else if (wish=="Lambda2_w11"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_w21"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_w31"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_w41"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1]=valRe; }
  else if (wish=="Lambda2_w12"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_w22"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_w32"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_w42"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2]=valRe; }
  else if (wish=="Lambda2_w10"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda2_w20"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda2_w30"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="Lambda2_w40"){ coupl_H.H2wwLambda_qsq[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12]=valRe; }
  else if (wish=="gh2w1"){ coupl_H.H2wwcoupl[gHIGGS_VV_1][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1][1]=valIm; }
  else if (wish=="gh2w2"){ coupl_H.H2wwcoupl[gHIGGS_VV_2][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2][1]=valIm; }
  else if (wish=="gh2w3"){ coupl_H.H2wwcoupl[gHIGGS_VV_3][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3][1]=valIm; }
  else if (wish=="gh2w4"){ coupl_H.H2wwcoupl[gHIGGS_VV_4][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4][1]=valIm; }
  else if (wish=="gh2w1_prime"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME][1]=valIm; }
  else if (wish=="gh2w1_prime2"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME2][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME2][1]=valIm; }
  else if (wish=="gh2w1_prime3"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME3][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME3][1]=valIm; }
  else if (wish=="gh2w1_prime4"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME4][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME4][1]=valIm; }
  else if (wish=="gh2w1_prime5"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME5][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME5][1]=valIm; }
  else if (wish=="gh2w2_prime"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME][1]=valIm; }
  else if (wish=="gh2w2_prime2"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME2][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME2][1]=valIm; }
  else if (wish=="gh2w2_prime3"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME3][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME3][1]=valIm; }
  else if (wish=="gh2w2_prime4"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME4][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME4][1]=valIm; }
  else if (wish=="gh2w2_prime5"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME5][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME5][1]=valIm; }
  else if (wish=="gh2w3_prime"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME][1]=valIm; }
  else if (wish=="gh2w3_prime2"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME2][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME2][1]=valIm; }
  else if (wish=="gh2w3_prime3"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME3][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME3][1]=valIm; }
  else if (wish=="gh2w3_prime4"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME4][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME4][1]=valIm; }
  else if (wish=="gh2w3_prime5"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME5][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME5][1]=valIm; }
  else if (wish=="gh2w4_prime"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME][1]=valIm; }
  else if (wish=="gh2w4_prime2"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME2][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME2][1]=valIm; }
  else if (wish=="gh2w4_prime3"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME3][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME3][1]=valIm; }
  else if (wish=="gh2w4_prime4"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME4][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME4][1]=valIm; }
  else if (wish=="gh2w4_prime5"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME5][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME5][1]=valIm; }
  else if (wish=="gh2w1_prime6"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME6][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME6][1]=valIm; }
  else if (wish=="gh2w1_prime7"){ coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME7][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_1_PRIME7][1]=valIm; }
  else if (wish=="gh2w2_prime6"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME6][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME6][1]=valIm; }
  else if (wish=="gh2w2_prime7"){ coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME7][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_2_PRIME7][1]=valIm; }
  else if (wish=="gh2w3_prime6"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME6][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME6][1]=valIm; }
  else if (wish=="gh2w3_prime7"){ coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME7][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_3_PRIME7][1]=valIm; }
  else if (wish=="gh2w4_prime6"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME6][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME6][1]=valIm; }
  else if (wish=="gh2w4_prime7"){ coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME7][0]=valRe; coupl_H.H2wwcoupl[gHIGGS_VV_4_PRIME7][1]=valIm; }

  // Spin-1 couplings, only for the first resonance
  else if (wish=="zprime_qq_left"){ coupl_Zp.Zqqcoupl[gZPRIME_QQ_LEFT][0]=valRe; coupl_Zp.Zqqcoupl[gZPRIME_QQ_LEFT][1]=valIm; }
  else if (wish=="zprime_qq_right"){ coupl_Zp.Zqqcoupl[gZPRIME_QQ_RIGHT][0]=valRe; coupl_Zp.Zqqcoupl[gZPRIME_QQ_RIGHT][1]=valIm; }
  else if (wish=="zprime_zz_1"){ coupl_Zp.Zvvcoupl[gZPRIME_VV_1][0]=valRe; coupl_Zp.Zvvcoupl[gZPRIME_VV_1][1]=valIm; }
  else if (wish=="zprime_zz_2"){ coupl_Zp.Zvvcoupl[gZPRIME_VV_2][0]=valRe; coupl_Zp.Zvvcoupl[gZPRIME_VV_2][1]=valIm; }

  // Spin-2 couplings, only for the first resonance
  else if (wish=="graviton_qq_left"){ coupl_X.Gqqcoupl[gGRAVITON_QQ_LEFT][0]=valRe; coupl_X.Gqqcoupl[gGRAVITON_QQ_LEFT][1]=valIm; }
  else if (wish=="graviton_qq_right"){ coupl_X.Gqqcoupl[gGRAVITON_QQ_RIGHT][0]=valRe; coupl_X.Gqqcoupl[gGRAVITON_QQ_RIGHT][1]=valIm; }
  else if (wish=="a1"){ coupl_X.Gggcoupl[gGRAVITON_GG_1][0]=valRe; coupl_X.Gggcoupl[gGRAVITON_GG_1][1]=valIm; }
  else if (wish=="a2"){ coupl_X.Gggcoupl[gGRAVITON_GG_2][0]=valRe; coupl_X.Gggcoupl[gGRAVITON_GG_2][1]=valIm; }
  else if (wish=="a3"){ coupl_X.Gggcoupl[gGRAVITON_GG_3][0]=valRe; coupl_X.Gggcoupl[gGRAVITON_GG_3][1]=valIm; }
  else if (wish=="a4"){ coupl_X.Gggcoupl[gGRAVITON_GG_4][0]=valRe; coupl_X.Gggcoupl[gGRAVITON_GG_4][1]=valIm; }
  else if (wish=="a5"){ coupl_X.Gggcoupl[gGRAVITON_GG_5][0]=valRe; coupl_X.Gggcoupl[gGRAVITON_GG_5][1]=valIm; }
  else if (wish=="b1"){ coupl_X.Gvvcoupl[gGRAVITON_VV_1][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_1][1]=valIm; }
  else if (wish=="b2"){ coupl_X.Gvvcoupl[gGRAVITON_VV_2][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_2][1]=valIm; }
  else if (wish=="b3"){ coupl_X.Gvvcoupl[gGRAVITON_VV_3][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_3][1]=valIm; }
  else if (wish=="b4"){ coupl_X.Gvvcoupl[gGRAVITON_VV_4][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_4][1]=valIm; }
  else if (wish=="b5"){ coupl_X.Gvvcoupl[gGRAVITON_VV_5][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_5][1]=valIm; }
  else if (wish=="b6"){ coupl_X.Gvvcoupl[gGRAVITON_VV_6][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_6][1]=valIm; }
  else if (wish=="b7"){ coupl_X.Gvvcoupl[gGRAVITON_VV_7][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_7][1]=valIm; }
  else if (wish=="b8"){ coupl_X.Gvvcoupl[gGRAVITON_VV_8][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_8][1]=valIm; }
  else if (wish=="b9"){ coupl_X.Gvvcoupl[gGRAVITON_VV_9][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_9][1]=valIm; }
  else if (wish=="b10"){ coupl_X.Gvvcoupl[gGRAVITON_VV_10][0]=valRe; coupl_X.Gvvcoupl[gGRAVITON_VV_10][1]=valIm; }

  // Spin-0 - V' contact terms
  else if (wish=="M_Zprime"){ coupl_Vprime.M_Zprime=valRe; }
  else if (wish=="Ga_Zprime"){ coupl_Vprime.Ga_Zprime=valRe; }
  else if (wish=="M_Wprime"){ coupl_Vprime.M_Wprime=valRe; }
  else if (wish=="Ga_Wprime"){ coupl_Vprime.Ga_Wprime=valRe; }

  else if (wish=="ezp_El_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_El_left, valRe, valIm, false); }
  else if (wish=="ezp_El_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_El_right, valRe, valIm, false); }
  else if (wish=="ezp_Mu_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Mu_left, valRe, valIm, false); }
  else if (wish=="ezp_Mu_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Mu_right, valRe, valIm, false); }
  else if (wish=="ezp_Ta_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Ta_left, valRe, valIm, false); }
  else if (wish=="ezp_Ta_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Ta_right, valRe, valIm, false); }
  else if (wish=="ezp_NuE_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_NuE_left, valRe, valIm, false); }
  else if (wish=="ezp_NuE_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_NuE_right, valRe, valIm, false); }
  else if (wish=="ezp_Dn_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Dn_left, valRe, valIm, false); }
  else if (wish=="ezp_Dn_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Dn_right, valRe, valIm, false); }
  else if (wish=="ezp_Up_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Up_left, valRe, valIm, false); }
  else if (wish=="ezp_Up_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Up_right, valRe, valIm, false); }
  else if (wish=="ezp_Str_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Str_left, valRe, valIm, false); }
  else if (wish=="ezp_Str_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Str_right, valRe, valIm, false); }
  else if (wish=="ezp_Chm_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Chm_left, valRe, valIm, false); }
  else if (wish=="ezp_Chm_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Chm_right, valRe, valIm, false); }
  else if (wish=="ezp_Bot_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Bot_left, valRe, valIm, false); }
  else if (wish=="ezp_Bot_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Bot_right, valRe, valIm, false); }
  else if (wish=="ezp_Top_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Top_left, valRe, valIm, false); }
  else if (wish=="ezp_Top_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Top_right, valRe, valIm, false); }

  else if (wish=="ewp_El_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_El_left, valRe, valIm, true); }
  else if (wish=="ewp_El_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_El_right, valRe, valIm, true); }
  else if (wish=="ewp_Mu_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Mu_left, valRe, valIm, true); }
  else if (wish=="ewp_Mu_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Mu_right, valRe, valIm, true); }
  else if (wish=="ewp_Ta_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Ta_left, valRe, valIm, true); }
  else if (wish=="ewp_Ta_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Ta_right, valRe, valIm, true); }
  //else if (wish=="ewp_NuE_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_NuE_left, valRe, valIm, true); }
  //else if (wish=="ewp_NuE_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_NuE_right, valRe, valIm, true); }
  //else if (wish=="ewp_Dn_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Dn_left, valRe, valIm, true); }
  //else if (wish=="ewp_Dn_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Dn_right, valRe, valIm, true); }
  else if (wish=="ewp_Up_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Up_left, valRe, valIm, true); }
  else if (wish=="ewp_Up_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Up_right, valRe, valIm, true); }
  //else if (wish=="ewp_Str_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Str_left, valRe, valIm, true); }
  //else if (wish=="ewp_Str_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Str_right, valRe, valIm, true); }
  else if (wish=="ewp_Chm_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Chm_left, valRe, valIm, true); }
  else if (wish=="ewp_Chm_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Chm_right, valRe, valIm, true); }
  //else if (wish=="ewp_Bot_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Bot_left, valRe, valIm, true); }
  //else if (wish=="ewp_Bot_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Bot_right, valRe, valIm, true); }
  else if (wish=="ewp_Top_left"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Top_left, valRe, valIm, true); }
  else if (wish=="ewp_Top_right"){ coupl_Vprime.SetVpffCouplings(gHIGGS_Vp_Top_right, valRe, valIm, true); }

  // aTQGC
  else if (wish=="dV_A"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dVA][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dVA][1]=valIm; }
  else if (wish=="dP_A"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dPA][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dPA][1]=valIm; }
  else if (wish=="dM_A"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dMA][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dMA][1]=valIm; }
  else if (wish=="dFour_A"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dFourA][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dFourA][1]=valIm; }

  else if (wish=="dV_Z"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dVZ][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dVZ][1]=valIm; }
  else if (wish=="dP_Z"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dPZ][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dPZ][1]=valIm; }
  else if (wish=="dM_Z"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dMZ][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dMZ][1]=valIm; }
  else if (wish=="dFour_Z"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dFourZ][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dFourZ][1]=valIm; }

  else if (wish=="dAAWpWm"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dAAWpWm][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dAAWpWm][1]=valIm; }
  else if (wish=="dZAWpWm"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dZAWpWm][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dZAWpWm][1]=valIm; }
  else if (wish=="dZZWpWm"){ coupl_aTQGC.aTQGCcoupl[gATQGC_dZZWpWm][0]=valRe; coupl_aTQGC.aTQGCcoupl[gATQGC_dZZWpWm][1]=valIm; }


  else cerr << "MELAOptionParser::extractCoupling: Coupling " << wish << " is not supported!" << endl;
}

void MELAOptionParser::setAddedAliases(string opt){
  string wish, aliases;
  splitOption(opt, wish, aliases, '=');
  splitOptionRecursive(aliases, addedAliases, ',');
}
void MELAOptionParser::setSubtractedAliases(string opt){
  string wish, aliases;
  splitOption(opt, wish, aliases, '=');
  splitOptionRecursive(aliases, subtractedAliases, ',');
}
void MELAOptionParser::setMultipliedAliases(string opt){
  string wish, aliases;
  splitOption(opt, wish, aliases, '=');
  splitOptionRecursive(aliases, multipliedAliases, ',');
}
void MELAOptionParser::setDividedAliases(string opt){
  string wish, aliases;
  splitOption(opt, wish, aliases, '=');
  splitOptionRecursive(aliases, dividedAliases, ',');
}

void MELAOptionParser::setMaximizationNumAliases(string opt){
  string wish, aliases;
  splitOption(opt, wish, aliases, '=');
  splitOptionRecursive(aliases, maximizationNumerators, ',');
}
void MELAOptionParser::setMaximizationDenomAliases(string opt){
  string wish, aliases;
  splitOption(opt, wish, aliases, '=');
  splitOptionRecursive(aliases, maximizationDenominators, ',');
}

void MELAOptionParser::pickOriginalOptions(MELAOptionParser* original_opt){
  // Check if name is the same. If so, append "_Copy".
  if (strName == original_opt->getName()) strName.append("_Copy");

  // Replace these values
  strCluster = original_opt->strCluster; // The cluster better be the same, overwrite it!

  //noBranching = original_opt->noBranching;
  includePAux = original_opt->includePAux;
  includePConst = original_opt->includePConst;
  isPM4L = original_opt->isPM4L;
  isPMaVJJ = original_opt->isPMaVJJ;
  isPMaVJJTrue = original_opt->isPMaVJJTrue;
  isProp = original_opt->isProp;
  isGenProb = original_opt->isGenProb;
  defME = original_opt->defME;

  hmass = original_opt->hmass;
  h2mass = original_opt->h2mass;
  hwidth = original_opt->hwidth;
  h2width = original_opt->h2width;

  coupl_H.copy(original_opt->coupl_H);
  coupl_Vprime.copy(original_opt->coupl_Vprime);
  coupl_Zp.copy(original_opt->coupl_Zp);
  coupl_X.copy(original_opt->coupl_X);
  coupl_aTQGC.copy(original_opt->coupl_aTQGC);

  couplingsString = original_opt->couplingsString;

  proc = original_opt->proc;
  prod = original_opt->prod;
  ME = original_opt->ME;
  superSyst = original_opt->superSyst;
  propScheme = original_opt->propScheme;

  // Append these arrays instead of replacing them
  for (unsigned int it=0; it<original_opt->addedAliases.size(); it++){ if (!checkListVariable(addedAliases, (original_opt->addedAliases).at(it))) addedAliases.push_back((original_opt->addedAliases).at(it)); }
  for (unsigned int it=0; it<original_opt->subtractedAliases.size(); it++){ if (!checkListVariable(subtractedAliases, (original_opt->subtractedAliases).at(it))) subtractedAliases.push_back((original_opt->subtractedAliases).at(it)); }
  for (unsigned int it=0; it<original_opt->multipliedAliases.size(); it++){ if (!checkListVariable(multipliedAliases, (original_opt->multipliedAliases).at(it))) multipliedAliases.push_back((original_opt->multipliedAliases).at(it)); }
  for (unsigned int it=0; it<original_opt->dividedAliases.size(); it++){ if (!checkListVariable(dividedAliases, (original_opt->dividedAliases).at(it))) dividedAliases.push_back((original_opt->dividedAliases).at(it)); }

  for (unsigned int it=0; it<original_opt->maximizationNumerators.size(); it++){ if (!checkListVariable(maximizationNumerators, (original_opt->maximizationNumerators).at(it))) maximizationNumerators.push_back((original_opt->maximizationNumerators).at(it)); }
  for (unsigned int it=0; it<original_opt->maximizationDenominators.size(); it++){ if (!checkListVariable(maximizationDenominators, (original_opt->maximizationDenominators).at(it))) maximizationDenominators.push_back((original_opt->maximizationDenominators).at(it)); }
}


