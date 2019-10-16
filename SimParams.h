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
    double               a_;      //half of rod length [m]
    double               Delta_;  //shortest distance between the tips of two nn rods [m]
    double               Ms_;     //saturation magnetization [A/m]
    double               rodRad_; //rod radius [m]
    ulong                seed_;
    uint                 numWarmUpSweeps_;
    uint                 sweepsPerMeas_;
    uint                 measPerBin_;
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
