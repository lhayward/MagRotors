/**********************************************************************************************
* Lauren E. Hayward
***********************************************************************************************
* File:   SimParams.cpp
**********************************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "FileReading.h"
#include "SimParams.h"

//typdefs needed because uint and ulong are return types:
typedef SimParams::uint  uint;
typedef SimParams::ulong ulong;

/********** SimParams(std::string fileName, std::string startStr) (constructor) **********/
SimParams::SimParams(std::string fileName, std::string startStr)
{
  const char EQUALS_CHAR = '=';
  const char LIST_START_CHAR = '[';
  const char LIST_END_CHAR = ']';
  
  TList_ = new std::vector<double>;
  
  std::ifstream fin;
  fin.open(fileName.c_str());
  
  if( fin.is_open() )
  {
    FileReading::readUntilFound(&fin, startStr);
    
    TList_           = FileReading::readDoubleVec(&fin, EQUALS_CHAR, LIST_START_CHAR, 
                                                  LIST_END_CHAR);
    Lx_              = FileReading::readUint     (&fin, EQUALS_CHAR);
    h_               = FileReading::readDouble   (&fin, EQUALS_CHAR);
    a_               = FileReading::readDouble   (&fin, EQUALS_CHAR);
    Delta_           = FileReading::readDouble   (&fin, EQUALS_CHAR);
    Ms_              = FileReading::readDouble   (&fin, EQUALS_CHAR);
    rodRad_          = FileReading::readDouble   (&fin, EQUALS_CHAR);
    seed_            = FileReading::readULong    (&fin, EQUALS_CHAR);
    numWarmUpSweeps_ = FileReading::readUint     (&fin, EQUALS_CHAR);
    sweepsPerMeas_   = FileReading::readUint     (&fin, EQUALS_CHAR);
    measPerBin_      = FileReading::readUint     (&fin, EQUALS_CHAR);
    numBins_         = FileReading::readUint     (&fin, EQUALS_CHAR);
    printSpins_      = FileReading::readBool     (&fin, EQUALS_CHAR);
    
    randomGen_.seed(seed_);
  }
  else
  { 
    std::cout << "ERROR in SimParams constructor: could not find file \"" << fileName 
              << "\"\n" << std::endl; 
  }
  
  fin.close();
} //SimParams constructor

/******************************* ~SimParams() (destructor) *******************************/
SimParams::~SimParams()
{
  if(TList_!=NULL)
  { delete TList_; }
  TList_ = NULL; 
}

/****************************************** print() ******************************************/
void SimParams::print()
{
  std::cout << "Simulation Parameters:\n"
            << "---------------------\n"
            << "           Temperature List = [ ";

  //print the list of temperatures:
  for( uint i=0; i<(TList_->size() - 1); i++ )
  { std::cout << (*TList_)[i] << ", "; }

  //print the last temperature element:
  if( TList_->size() > 0 )
  { std::cout << (*TList_)[TList_->size() - 1]; }
  std::cout << " ]\n";

  //print the rest of the parameters:
  std::cout << "        Length of 1D chains = " << Lx_ << "\n"
            << "                          h = " << h_ << "\n"
            << "                          a = " << a_ << "\n"
            << "                      Delta = " << Delta_ << "\n"
            << "                         Ms = " << Ms_ << "\n"
            << "                 Rod radius = " << rodRad_ << "\n"
            << "                       Seed = " << seed_ << "\n"
            << "   Number of Warm-up Sweeps = " << numWarmUpSweeps_ << "\n"
            << "     Sweeps per Measurement = " << sweepsPerMeas_ << "\n"
            << "       Measurements per Bin = " << measPerBin_ << "\n"
            << "             Number of Bins = " << numBins_ << "\n"
            << "       Print Spins Configs? = " << printSpins_ << "\n"
            << std::endl;   
  
} //print method

