#include "MELACluster.h"

using namespace std;


MELACluster::MELACluster(string name_) : name(name_){}
MELACluster::~MELACluster(){}
void MELACluster::addComputation(MELAComputation* comp){ computers.push_back(comp); }
void MELACluster::computeAll(){
  // Re-compute all related hypotheses
  for (MELAComputation* const& computer:computers){
    // Avoid re-computing MEs twice (could happen through copy-computations)
    if (!computer->getOption()->isCopy()) computer->getHypothesis()->computeP();
  }
}
void MELACluster::update(){ for (MELAComputation* const& computer:computers) computer->update(); }
void MELACluster::forceUpdate(){ for (MELAComputation* const& computer:computers) computer->forceUpdate(); }
void MELACluster::reset(){ for (MELAComputation* const& computer:computers) computer->reset(); }

