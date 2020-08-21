# MELA Analytics package

- The EventContainer package provides a clean and general way to organize the LHE-level particles.

- The GenericMEComputer package provides the capability to compute MEs from MELA using one-line strings.

- The CandidateLOCaster package allows the user to cast some topologies NLO in QCD into those LO in QCD by merging gluons into quarks appropriately.
	It is a pretty good approxiamtion for EW processes, but does not provide a good enough description for QCD processes themselves.

## Checkout instructions

### MELA

```
git clone https://github.com/JHUGen/JHUGenMELA.git
(cd JHUGenMELA; git checkout -b from-v231 v2.3.1; ./setup.sh -j;)
# Last line could also use './setup.sh -j standalone' in order to force standalone computation.
```

### MELA Analytics

```
git clone https://github.com/MELALabs/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v22 v2.2)
```

## Compilation

First, please follow the instructions on MELA core package to compile it. Assuming you already did that,

- if you are using CMSSW, regular scram should work.
- if you are not using CMSSW, each package contains its own makefile. You can simply do

```
make (-j #cores)
```
