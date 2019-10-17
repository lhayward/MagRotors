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
#include <math.h>
#include <sstream>
#include <string>
#include <typeinfo>
#include "FileReading.h"
#include "MersenneTwister.h"
#include "SimParams.h"

typedef unsigned long ulong;
typedef unsigned int  uint;


//Function definitions:
double      getEnergy();
std::string getFileSuffix(int argc, char** argv);
void        initializeDistances();
void        localUpdate(MTRand &randomGen);
void        localUpdate_A(MTRand &randomGen);
void        localUpdate_B(MTRand &randomGen);
void        printDoubleArray( double* arr, std::string name, uint len );
void        print2DDoubleArray( double** arr, uint L1, uint L2 );
void        sweep(MTRand &randomGen, uint N);

//Global constants:
const double PI   = 3.14159265358979323846;
const double MU_0 = 4*PI*(1e-7);

//Global variables:
uint         Lx;          //stored in params but let's store again since it's used so often
double       d;           //distance between the centres of magnets
double       q;           //magnitude of the magnetic charge at the tip of each rod
double       m;           //magnitude of magnetic moment
double       D_dip;       //dipolar energy scale
double**     r_AA;        //distances between pairs of magnets (both in A) [multiply by d to get m]
double**     r_BB;        //distances between pairs of magnets (both in B) [multiply by d to get m]
double**     r_AB;        //distances between pairs of magnets (first one in A, second one in B) [multiply by d to get m]
double**     r_AB_x;      //x-components of r_AB [multiply by d to get m]
double       r_AB_y;      //y-component of r_AB (same for all pairs) [multiply by d to get m]
double*      alpha_A;     //array of angles corresponding to the magnets in sublattice A
double*      alpha_B;     //array of angles corresponding to the magnets in sublattice A
double       T;           //current temperature
std::string  T_str;       //T as a string

/**********************************************************************************************
******************************************** main *********************************************
**********************************************************************************************/
int main(int argc, char** argv) 
{
  //Parameters needed to run the simulation:
  SimParams*   params;      //store the simulation params read from file in SimParams object
  uint         N_spins;     //total number of spins in both sublattices
  uint         numT;        //number of temperatures
  time_t       sec1, sec2;  //for timing
  
  //Variables related to input/output data from/to files:
  std::string        fileSuffix    = getFileSuffix( argc, argv );
  std::string        paramFileName = "params" + fileSuffix + ".txt";
  std::string        simParamStr   = "SIMULATION PARAMETERS";
  std::string        binsFileName;
  std::string        spinsFileName;
  std::ofstream      fout_bins;
  std::ofstream      fout_spins;
  std::ostringstream ss;
  
  //Read in parameters from file:
  std::cout.precision(8);
  std::cout << "\nParameter File: " << paramFileName << "\n" << std::endl;
  params = new SimParams( paramFileName, simParamStr );
  Lx     = params->Lx_;
  params->print();
  
  //Variables that store observables:
  uint   N_meas;        //number of measurements averaged over
  double e;              //current energy per spin
  double sum_e, sum_eSq; //sums of e and e^2 over measurements
  
  //Define other parameters based on the ones read from file:
  N_spins = 2*Lx;
  d = (2*params->a_) + params->Delta_;
  q = PI*pow(params->rodRad_,2)*params->Ms_;
  m = 2*params->a_*q;
  D_dip = MU_0*pow(m,2)/(PI*pow(d,3));
  r_AB_y = params->h_;
  
  std::cout << D_dip << std::endl;
  
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
  
  initializeDistances();
  
  print2DDoubleArray( r_AA, Lx, Lx );
  print2DDoubleArray( r_BB, Lx, Lx );
  print2DDoubleArray( r_AB, Lx, Lx );
  
  
  //Initialize the angles in sublattice A to be random:
  alpha_A = new double[Lx];
  alpha_B = new double[Lx];
  for( uint i=0; i<Lx; i++ )
  {
    alpha_A[i] = params->randomGen_.randDblExc( 2*PI );
    //alpha_A[i] = 0;
    alpha_B[i] = params->randomGen_.randDblExc( 2*PI );
    //alpha_B[i] = PI;
  }
  
  std::cout << "\n*** STARTING SIMULATION ***\n" << std::endl;
  sec1 = time (NULL);
  
  //loop over the different temperatures:
  numT = params->TList_->size();
  for( uint TIndex=0; TIndex<numT; TIndex++)
  {
    T = (*(params->TList_))[TIndex];
    std::cout << "******** T = " << T << " (Temperature #" << (TIndex+1) << ") ********"
              << std::endl;
    
    ss.str("");
    ss << T;
    T_str = ss.str();
    
    //Open file to write bins:
    binsFileName = "bins" + fileSuffix + "_T" + T_str + ".txt";
    fout_bins.open(binsFileName);
    fout_bins.precision(15);
    fout_bins << "# binNum \t E \t ESq" << std::endl;
    
    if( params->printSpins_ )
    {
      spinsFileName = "spins" + fileSuffix + "_T" + T_str + ".txt";
      fout_spins.open(spinsFileName);
      fout_spins.precision(15);
    }
    
    //equilibrate:
    for( uint i=0; i<params->numWarmUpSweeps_; i++ )
    { sweep( params->randomGen_, N_spins ); }
    
    
    //loop over Monte Carlo bins:
    for( uint i=0; i<params->numBins_; i++ )
    {
      sum_e = 0;
      sum_eSq = 0;
      //perform the measurements for one bin:
      for( uint j=0; j<params->measPerBin_; j++ )
      {
        //perform the sweeps for one measurement:
        for( uint k=0; k<params->sweepsPerMeas_; k++ )
        { sweep( params->randomGen_, N_spins ); }
        
        //make energy measurements:
        e = getEnergy()/(1.0*N_spins);
        sum_e   += e;
        sum_eSq += pow(e,2);
      } //loop over measurements
      
      //Write binned measurements to file:
      N_meas = params->measPerBin_;
      fout_bins << (i+1)                << '\t'
                << sum_e/(1.0*N_meas)   << '\t'
                << sum_eSq/(1.0*N_meas) << std::endl;
      
      //std::cout << "sum_e = "   << sum_e   << std::endl;
      //std::cout << "sum_eSq = " << sum_eSq << std::endl << std::endl;
      
      printDoubleArray( alpha_A, "alpha_A", Lx);
      printDoubleArray( alpha_B, "alpha_B", Lx);
      std::cout << "energy = " << getEnergy()/(1.0*N_spins) << std::endl;
      std::cout << std::endl;
      //Write current spin configuration to file:
      if( params->printSpins_ )
      {
        //Print the spins in sublattice A:
        for(uint iA=0; iA<Lx; iA++)
        { fout_spins << alpha_A[iA] << " ";  }
        
        //Print the spins in sublattice B:
        for(uint iB=0; iB<Lx; iB++)
        { fout_spins << alpha_B[iB] << " ";  }
        fout_spins << std::endl;
      }
      
      
      if( (i+1)%100==0 )
      { std::cout << (i+1) << " Bins Complete" << std::endl; }
      
    } //loop over bins
    
    std::cout << std::endl;
    fout_bins.close();
    
    if( params->printSpins_ )
    { fout_spins.close(); }
  } //temperature loop
  
  sec2 = time(NULL);
  std::cout << "Time: " << (sec2 - sec1) << " seconds" << std::endl;
  std::cout << "\n*** END OF SIMULATION ***\n" << std::endl;
  return 0;
} //closes main

/*************************************** getEnergy ****************************************/
double getEnergy()
{
  double E_AA = 0;
  double E_BB = 0;
  double E_AB = 0;
  
  //Calculate E_AA and E_BB:
  for(uint i=0; i<Lx; i++)
  {
    for(uint j=(i+1); j<Lx; j++)
    {
      E_AA += ( cos(alpha_A[i])*cos(alpha_A[j]) - 2*sin(alpha_A[i])*sin(alpha_A[j]) )/pow(r_AA[i][j],3);
      E_BB += ( cos(alpha_B[i] - alpha_B[j])  )/pow(r_BB[i][j],3);
    }
  }
  
  //Calculate E_AB:
  for(uint iA=0; iA<Lx; iA++)
  {
    for(uint jB=0; jB<Lx; jB++)
    {
      E_AB += ( cos(alpha_A[iA])*cos(alpha_B[jB])/pow(r_AB[iA][jB],3) )
              - ( 3.0*r_AB_x[iA][jB]*r_AB_y*sin(alpha_A[iA])*sin(alpha_B[jB])/pow(r_AB[iA][jB],5) );
    }
  }
  
  return D_dip/4.0*(E_AA + E_BB + E_AB);
}

/*************************************** getFileSuffix ****************************************/
std::string getFileSuffix(int argc, char** argv)
{
  std::string result = "";
  
  if( argc > 1 )
  { result = "_" + std::string(argv[1]); }
  
  return result;
}

/************************************ initializeDistances *************************************/
void initializeDistances()
{
  //Calculate the distances r_AA and r_BB assuming PBC:
  double r1, r2;  //distance along the two lattice directions (because of PBC)
  for( uint i=0; i<Lx; i++ )
  {
    r_AA[i][i] = 0;
    r_BB[i][i] = 0;
    
    for( uint j=(i+1); j<Lx; j++ )
    {
      r1 = abs(int(j-i));
      r2 = Lx - r1;
      
      r_AA[i][j] = std::min(r1,r2);
      r_AA[j][i] = r_AA[i][j];
      r_BB[i][j] = r_AA[i][j];
      r_BB[j][i] = r_BB[i][j];
    } //j
  } //i

//  //Calculate the distances r_AA and r_BB assuming OBC:
//  for( uint i=0; i<Lx; i++ )
//  {
//    r_AA[i][i] = 0;
//    r_BB[i][i] = 0;
//    
//    for( uint j=(i+1); j<Lx; j++ )
//    {
//      r_AA[i][j] = j-i;
//      r_AA[j][i] = r_AA[i][j];
//      r_BB[i][j] = r_AA[i][j];
//      r_BB[j][i] = r_BB[i][j];
//    } //j
//  } //i
  
  //Calculate the distances r_AB and r_AB_x assuming PBC:
  double rx1, rx2; //x distance along the two lattice directions (because of PBC)
  for( uint iA=0; iA<Lx; iA++ )
  {
    for( uint jB=0; jB<Lx; jB++ )
    {
      rx1 = ( abs(jB+0.5-iA) );
      rx2 = Lx - rx1;
      
      r_AB_x[iA][jB] = std::min(rx1,rx2);
      r_AB[iA][jB]   = pow( pow(r_AB_x[iA][jB],2) + pow(r_AB_y,2) , 0.5);
    } //j
  } //i
}

/**************************************** localUpdate *****************************************/
void localUpdate(MTRand &randomGen)
{
  if( randomGen.randInt(1) )
  { localUpdate_A(randomGen); }
  else
  { localUpdate_B(randomGen); }
}
  
/*************************************** localUpdate_A ****************************************/
void localUpdate_A(MTRand &randomGen)
{
  uint         site   = randomGen.randInt(Lx-1); //randomly selected spin location
  double       alphaFrac = 0.1;                  //fraction of circle on while to consider changes
  double       dAlpha = alphaFrac*2.0*PI*( randomGen.randDblExc() - 0.5 ); //proposed change in alpha
  double       E_i, E_f, dE;      //initial E, final E, and change in E of the local move
  
  //Calculate the change in energy:
  E_i = getEnergy();
  alpha_A[site] += dAlpha;
  E_f = getEnergy();
  alpha_A[site] -= dAlpha;
  dE = E_f - E_i;
  
  //Check if the move is accepted:
  if( dE<=0 || randomGen.randDblExc() < exp(-dE/T) )
  { 
    alpha_A[site] += dAlpha;
    
    //Ensure alpha: is between 0 and 2*PI:
    while( alpha_A[site] < 0 )
    { alpha_A[site] += 2*PI; }
    
    while( alpha_A[site] > 2*PI )
    { alpha_A[site] -= 2*PI; }
  }
}

/*************************************** localUpdate_B ****************************************/
void localUpdate_B(MTRand &randomGen)
{
  uint         site   = randomGen.randInt(Lx-1); //randomly selected spin location
  double       alphaFrac = 0.1;                  //fraction of circle on while to consider changes
  double       dAlpha = alphaFrac*2.0*PI*( randomGen.randDblExc() - 0.5 ); //proposed change in alpha
  double       E_i, E_f, dE;      //initial E, final E, and change in E of the local move
  
  //Calculate the change in energy:
  E_i = getEnergy();
  alpha_B[site] += dAlpha;
  E_f = getEnergy();
  alpha_B[site] -= dAlpha;
  dE = E_f - E_i;
  
  //Check if the move is accepted:
  if( dE<=0 || randomGen.randDblExc() < exp(-dE/T) )
  { 
    alpha_B[site] += dAlpha;
    
    //Ensure alpha: is between 0 and 2*PI:
    while( alpha_B[site] < 0 )
    { alpha_B[site] += 2*PI; }
    
    while( alpha_B[site] > 2*PI )
    { alpha_B[site] -= 2*PI; }
  }
}

/************************************** printDoubleArray **************************************/
void printDoubleArray( double* arr, std::string name, uint len )
{
  std::cout << name << ": [ ";
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
