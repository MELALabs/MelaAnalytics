/** \class MELACluster
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*
*
*  Description:
*
*  This class contains the matrix element computations
*  that belong to a specific choice of candidates.
*  
*
*/
#ifndef MELACLUSTER_H
#define MELACLUSTER_H

#include "MELAComputation.h"


class MELACluster{

protected:

  std::string name;
  std::vector<MELAComputation*> computers;

public:

  MELACluster(std::string name_);
  virtual ~MELACluster();

  std::string const& getName() const{ return name; }
  void addComputation(MELAComputation* comp);
  std::vector<MELAComputation*> const* getComputations() const{ return &computers; }
  void computeAll();
  void update();
  void forceUpdate();
  void reset();

};


#endif
