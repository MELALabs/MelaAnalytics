* MELA Analytics package *

- The EventContainer package provides a clean and general way to organize the LHE-level particles.

- The GenericMEComputer package provides the capability to compute MEs from MELA using one-line strings.

- The CandidateLOCaster package allows the user to cast some topologies NLO in QCD into those LO in QCD by merging gluons into quarks appropriately.
	It is a pretty good approxiamtion for EW processes, but does not provide a good enough description for QCD processes themselves.

* Checkout instructions *

# MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v221 v2.2.1; source setup.sh -j;)

# MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v15 v1.5)
