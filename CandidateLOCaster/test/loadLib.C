{
gSystem->Load("$CMSSW_BASE/src/JHUGenMELA/MELA/test/loadMELA.C");

gSystem->AddIncludePath("-I$CMSSW_BASE/src/MelaAnalytics/CandidateLOCaster/interface/");
gSystem->AddIncludePath("-I$CMSSW_BASE/src/MelaAnalytics/CandidateLOCaster/test/");
gSystem->Load("libMelaAnalyticsCandidateLOCaster.so");
}
