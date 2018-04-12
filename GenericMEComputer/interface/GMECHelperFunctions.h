#ifndef GMECHELPERFUNCTIONS_H
#define GMECHELPERFUNCTIONS_H

#include "Mela.h"
#include "MELABranch.h"
#include "MELACluster.h"


namespace GMECHelperFunctions{

  void addToMELACluster(MELAComputation* me_computer, std::vector<MELACluster*>& me_clusters);

}

#endif
