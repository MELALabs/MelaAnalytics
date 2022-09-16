# MELA Analytics package

- The EventContainer package provides a clean and general way to organize the LHE-level particles.

- The GenericMEComputer package provides the capability to compute MEs from MELA using one-line strings.

- The CandidateLOCaster package allows the user to cast some topologies NLO in QCD into those LO in QCD by merging gluons into quarks appropriately.
	It is a pretty good approxiamtion for EW processes, but does not provide a good enough description for QCD processes themselves.

## Checkout and compilation instructions

### MELA

```
git clone https://github.com/JHUGen/JHUGenMELA.git
(cd JHUGenMELA; git checkout -b from-v238 v2.3.8; ./setup.sh -j;)
eval $(./JHUGenMELA/setup.sh env)
```

### MELA Analytics

```
git clone https://github.com/MELALabs/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v23 v2.3)
eval $(./MelaAnalytics/setup.sh env)
```
