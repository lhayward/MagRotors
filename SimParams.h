/**********************************************************************************************
* Lauren E. Hayward
***********************************************************************************************
* File:   SimParams.h
**********************************************************************************************/

#ifndef SIM_PARAMETERS_H
#define SIM_PARAMETERS_H

#include <fstream>
#include <string>
#include <vector>
#include "MersenneTwister.h"

class SimParams
{ 
  public:
    typedef unsigned int  uint;
    typedef unsigned long ulong;
  
    //params to read from file:
    std::vector<double>* TList_;  //list of temperatures
    uint                 Lx_;     //length along one chain
    double               h_;      //y-distance between the sublattices [dimensionless]
    ulong                seed_;
    uint                 numWarmUpSweeps_;
    uint                 sweepsPerBin_;
    uint                 numBins_;
    bool                 printSpins_;
  
    //variables defined based on the read params:
    MTRand               randomGen_;  //random number generator
    
  //public:
    SimParams(std::string fileName, std::string startStr);
    virtual ~SimParams();
    
    void print();
};

#endif  // SIM_PARAMETERS_H
