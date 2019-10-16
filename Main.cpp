/**********************************************************************************************
************************************ MAGNETIC ROTORS CODE *************************************
***********************************************************************************************
* Lauren E. Hayward
***********************************************************************************************
* File: Main.cpp
* Note: MersenneTwister.h was taken from external sources and is used as the random number
*       generator
**********************************************************************************************/

//#include <cmath>
#include <ctime>
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include "FileReading.h"
#include "MersenneTwister.h"
#include "SimParams.h"

typedef unsigned long ulong;
typedef unsigned int  uint;


//function definitions:
std::string getFileSuffix(int argc, char** argv);
void        localUpdate(MTRand &randomGen);
void        printDoubleArray( double* arr, uint len );
void        print2DDoubleArray( double** arr, uint L1, uint L2 );
void        sweep(MTRand &randomGen, uint N);

//global variables:
double       d;           //distance between the centres of magnets
double       q;           //magnitude of the magnetic charge at the tip of each rod
double       m;           //magnitude of magnetic moment
double       D_dip;       //dipolar energy scale
double**     r_AA;        //distances between pairs magnets (both in A)
double**     r_BB;        //distances between pairs magnets (both in B)
double**     r_AB;        //distances between pairs magnets (one in A, one in B)
double**     r_AB_x;      //x-components of r_AB
double       r_AB_y;      //y-component of r_AB (same for all pairs)
double*      alpha_A;     //array of angles corresponding to the magnets in sublattice A
double*      alpha_B;     //array of angles corresponding to the magnets in sublattice A
double       T;           //current temperature

/**********************************************************************************************
******************************************** main *********************************************
**********************************************************************************************/
int main(int argc, char** argv) 
{
  //Constants:
  const double MU_0 = 1.0e-7;
  const double PI   = 3.14159265358979323846;
  
  //Parameters needed to run the simulation:
  SimParams*   params;      //store the simulation params read from file in SimParams object
  uint         Lx;          //stored in params but let's store again since it's used so often
  uint         N_spins;     //total number of spins in both sublattices
  uint         numT;        //number of temperatures
  time_t       sec1, sec2;  //for timing
  
  //Variables related to input/output data from/to files:
  std::string fileSuffix         = getFileSuffix( argc, argv );
  std::string paramFileName      = "params" + fileSuffix + ".txt";
  std::string simParamStr        = "SIMULATION PARAMETERS";
  std::string outFileName        = "bins" + fileSuffix + ".txt";
  std::string spinConfigFileName = "spinConfigs" + fileSuffix + ".txt";
  
  //Read in parameters from file:
  std::cout.precision(8);
  std::cout << "\nParameter File: " << paramFileName << "\n" << std::endl;
  params = new SimParams( paramFileName, simParamStr );
  Lx     = params->Lx_;
  params->print();
  
  //Define other parameters based on the ones read from file:
  N_spins = 2*Lx;
  d = (2*params->a_) + params->Delta_;
  q = PI*pow(params->rodRad_,2)*params->Ms_;
  m = 2*params->a_*q;
  D_dip = MU_0*pow(m,2)/(PI*pow(d,3));
  r_AB_y = params->h_*d;
  
  //Initialize the 2D distances arrays:
  r_AA   = new double*[Lx];
  r_BB   = new double*[Lx];
  r_AB   = new double*[Lx];
  r_AB_x = new double*[Lx];
  for( uint i=0; i<Lx; i++ )
  {
    r_AA[i]   = new double[Lx];
    r_BB[i]   = new double[Lx];
    r_AB[i]   = new double[Lx];
    r_AB_x[i] = new double[Lx];
  }
  
  //Calculate the distances r_AA and r_BB assuming PBC:
  double r1, r2;  //distance along the two lattice directions (because of PBC)
  for( uint i=0; i<Lx; i++ )
  {
    r_AA[i][i] = 0;
    r_BB[i][i] = 0;
    
    for( uint j=i; j<Lx; j++ )
    {
      r1 = abs(int(j-i))*d;
      r2 = Lx*d - r1;
      
      r_AA[i][j] = std::min(r1,r2);
      r_AA[j][i] = r_AA[i][j];
      r_BB[i][j] = r_AA[i][j];
      r_BB[j][i] = r_BB[i][j];
    } //j
  } //i
  
  //Calculate the distances r_AB and r_AB_x assuming PBC:
  double rx1, rx2; //x distance along the two lattice directions (because of PBC)
  for( uint i=0; i<Lx; i++ )
  {
    for( uint j=0; j<Lx; j++ )
    {
      rx1 = ( abs(int(j-i)) + 0.5 )*d;
      rx2 = Lx*d - rx1;
      
      r_AB_x[i][j] = std::min(rx1,rx2);
      r_AB[i][j]   = pow( pow(r_AB_x[i][j],2) + pow(r_AB_y,2) , 0.5);
    } //j
  } //i
  
  //Initialize the angles in sublattice A to be random:
  alpha_A = new double[Lx];
  alpha_B = new double[Lx];
  for( uint i=0; i<Lx; i++ )
  {
    alpha_A[i] = params->randomGen_.randDblExc( 2*PI );
    alpha_B[i] = params->randomGen_.randDblExc( 2*PI );
  }
  
  printDoubleArray(alpha_A,Lx);
  
  std::cout << "\n*** STARTING SIMULATION ***\n" << std::endl;
  sec1 = time (NULL);
  
  //loop over the different temperatures:
  numT = params->TList_->size();
  for( uint TIndex=0; TIndex<numT; TIndex++)
  {
    T = (*(params->TList_))[TIndex];
    std::cout << "******** T = " << T << " (Temperature #" << (TIndex+1) << ") ********"
              << std::endl;
    
    //equilibrate:
    for( uint i=0; i<params->numWarmUpSweeps_; i++ )
    { sweep( params->randomGen_, N_spins ); }
    
    /*
    //loop over Monte Carlo bins:
    for( uint i=0; i<params->numBins_; i++ )
    {
      model->zeroMeasurements();
      //perform the measurements for one bin:
      for( uint j=0; j<params->measPerBin_; j++ )
      {
        //perform the sweeps for one measurement:
        for( uint k=0; k<params->sweepsPerMeas_; k++ )
        { model->sweep( params->randomGen_, false ); }
        model->makeMeasurement();
      } //loop over measurements
      model->writeBin((i+1), params->measPerBin_, params->sweepsPerMeas_);
      
      if( params->printSpins_ )
      { model->writeSpins(); }
      if( (i+1)%100==0 )
      { std::cout << (i+1) << " Bins Complete" << std::endl; }
    } //loop over bins
    */
    
    std::cout << std::endl;
  } //temperature loop
  
  sec2 = time(NULL);
  std::cout << "Time: " << (sec2 - sec1) << " seconds" << std::endl;
  std::cout << "\n*** END OF SIMULATION ***\n" << std::endl;
  return 0;
} //closes main

/*************************************** getFileSuffix ****************************************/
std::string getFileSuffix(int argc, char** argv)
{
  std::string result = "";
  
  if( argc > 1 )
  { result = "_" + std::string(argv[1]); }
  
  return result;
}

/**************************************** localUpdate *****************************************/
void localUpdate(MTRand &randomGen)
{
  double r = randomGen.randDblExc();
  std::cout << alpha_A[0] << std::endl;
}

/************************************** printDoubleArray **************************************/
void printDoubleArray( double* arr, uint len )
{
  std::cout << "[ "; 
  for( uint i=0; i<(len-1); i++ )
  { std::cout << arr[i] << " , "; }
  std::cout << arr[len-1] << " ]" << std::endl;
}

/************************************* print2DDoubleArray *************************************/
void print2DDoubleArray( double** arr, uint L1, uint L2 )
{
  for( uint i=0; i<L1; i++ )
  {
    for( uint j=0; j<L2; j++ )
    {
      std::cout << arr[i][j] << " ";
    } //j
    std::cout << std::endl;
  } //i 
  std::cout << std::endl;
} //print2DDoubleArray method

/******************************************* sweep ********************************************/
void sweep(MTRand &randomGen, uint N)
{
  for( uint i=0; i<N; i++ )
  { localUpdate(randomGen); }
}
