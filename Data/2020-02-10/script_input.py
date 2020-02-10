################################ script_input.py ################################
#################################################################################

import math
import numpy as np
import os
import random
#from matplotlib import pyplot as plt

####################################### writeInputFile #######################################
def writeInputFile( fileDir, tag, TList, Lx, h, seed, numWarmup, sweepsPerBin, numBins, printSpins):
  if fileDir == "":
    inputFileName = "params_" + tag + ".txt"
  else:
    inputFileName = fileDir + "/params_" + tag + ".txt"
    
  #print inputFileName
  fout = open(inputFileName, 'w')
  
  fout.write("#SIMULATION PARAMETERS:\n")
  fout.write("Temperature List              = " + str(TList)         + '\n')
  fout.write("Lx                            = " + str(Lx)            + '\n')
  fout.write("h                             = " + str(h)            + '\n')
  fout.write("Random Generator Seed         = " + str(seed)          + '\n')
  fout.write("Number of Warm-up Sweeps      = " + str(numWarmup)     + '\n')
  fout.write("Number of Sweeps per Bin      = " + str(sweepsPerBin) + '\n')
  fout.write("Number of Bins                = " + str(numBins)       + '\n')
  fout.write("Print Spins?                  = " + str(printSpins)   + '\n')
  fout.write('\n')
  
  fout.close()
  

##############################################################################################
############################################ main ############################################
##############################################################################################
execProg = "mag"

bashFileName = "bash_script_input.sh"
fout = open(bashFileName, 'w')

#### SIMULATION PARAMETERS: ####
sweepsPerBin=1000
numBins=1000
numWarmup = int(0.1*sweepsPerBin*numBins)

tempList =  [0.1, 0.001, 0.0001]
Lx=10
hList = np.linspace(1.5,0.05,num=30).tolist()
#Lx=4
#hList = np.linspace(1.5,0.05,num=3).tolist()
print_spins = 0

LDir = "Lx" + str(Lx)
#Create the directory for this set of parameters:
if (LDir != "") and not(os.path.isdir(LDir)):
  os.mkdir(LDir)
fout.write( "cd " + LDir + "\n\n")
fout.write( "cp " + "../" + execProg + " . \n\n" )
    
#Loop over data sets:
for h in hList:
  
  tag = "h" + str('%1.2f'%h) + "_" + LDir
  seed = random.randint(0,1e10)
        
  writeInputFile( LDir, tag, tempList, Lx, h, seed, numWarmup, sweepsPerBin, numBins, print_spins )
    
  outFileName = "output_" + tag + ".txt"

  submitCmd = "  ./" + execProg + " " + tag + " > " + outFileName
  fout.write( "  echo " + '"' + submitCmd + '"' + "\n\n")
  fout.write(submitCmd + "\n")

#end of h loop

fout.write("cd .." + "\n\n")  #LDir

fout.close()
