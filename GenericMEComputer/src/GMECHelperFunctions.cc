#include "GMECHelperFunctions.h"


void GMECHelperFunctions::addToMELACluster(MELAComputation* me_computer, std::vector<MELACluster*>& me_clusters){
  bool isAdded=false;
  for (MELACluster* me_cluster:me_clusters){
    if (me_cluster->getName()==me_computer->getCluster()){
      me_cluster->addComputation(me_computer);
      isAdded=true;
      break;
    }
  }
  if (!isAdded){
    MELACluster* tmpcluster = new MELACluster(me_computer->getCluster());
    tmpcluster->addComputation(me_computer);
    me_clusters.push_back(tmpcluster);
  }
}
