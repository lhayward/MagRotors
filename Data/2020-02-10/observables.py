import matplotlib.pyplot as plt
import numpy as np
import time

MU_0 = 4*np.pi*(1e-7);

a       = 2.5e-3    #half of rod length [m]
Delta   = 4.5e-3    #shortest distance between the tips of two nn rods [m]
Ms      = 613310    #saturation magnetization [A/m]
rodRad  = 0.25e-3   #rod radius [m]

d = (2*a) + Delta         #distance between the centres of magnets
q = np.pi*(rodRad**2)*Ms  #magnitude of the magnetic charge at the tip of each rod
m = 2*a*q;                #magnitude of magnetic moment
D_dip = MU_0*m**2/(np.pi*d**3); #dipolar energy scale

#print(D_dip)

#---------

Lx = 8
N = 2*Lx
LDir = "Lx" + str(Lx)

tempList =  [0.1, 0.001, 0.0001]
T = min(tempList)
#T = 0.1

#h_list = np.linspace(1.5,0.05,num=30).tolist()
h1 = np.linspace(1.5,0.35,num=24)
h2 = np.linspace(0.3,0.01,num=30)
h_list = np.concatenate([h1,h2])

e_list      = []
Cv_list     = []

mAx_list    = []
mAz_list    = []
msAx_list   = []
msAz_list   = []

mBy_list    = []
mBz_list    = []
msBy_list   = []
msBz_list   = []
msB_list    = []

mAFpi_list  = []
mFAF_list   = []
mAFpi2_list = []

for h in h_list:

  fileName = LDir + "/bins_h" + str('%1.2f'%h) + "_" + LDir + "_T" + str(T) + ".txt"
  data = np.loadtxt( fileName, comments='#' )
  #print(data)

  e    = data[:,1]
  mAx  = data[:,2]
  mAz  = data[:,3]
  mBy  = data[:,4]
  mBz  = data[:,5]
  msAx = data[:,6]
  msAz = data[:,7]
  msBy = data[:,8]
  msBz = data[:,9]
  
  msB = np.sqrt( msBy**2 + msBz**2 )

  mAFpi  = np.abs(mAz - mBz)/2.0
  #mAFpi2 =
  mFAF   = (np.abs(mAx) + msB)/2.0
  
  e_list.append(np.mean(e))
  Cv_list.append(N**2*( np.mean(e**2) - np.mean(e)**2 )/T**2)
  
  mAx_list.append(np.mean(np.abs(mAx)))
  mAz_list.append(np.mean(np.abs(mAz)))
  msAx_list.append(np.mean(np.abs(msAx)))
  msAz_list.append(np.mean(np.abs(msAz)))
  
  mBy_list.append(np.mean(np.abs(mBy)))
  mBz_list.append(np.mean(np.abs(mBz)))
  msBy_list.append(np.mean(np.abs(msBy)))
  msBz_list.append(np.mean(np.abs(msBz)))
  
  
  msB_list.append(np.mean(msB))
  mAFpi_list.append(np.mean(mAFpi))
  mFAF_list. append(np.mean(mFAF))

#end of h loop

plt.figure()
plt.plot(h_list, e_list, 'o-')
plt.title(r'$L_x = %d$, $T = %f$' %(Lx,T))
plt.xlabel(r'$h$')
plt.ylabel(r'$E/N$')
plt.savefig('e_vs_h_' + LDir + '_T' + str(T) + '.pdf')

plt.figure()
plt.plot(h_list, Cv_list, 'o-')
plt.title(r'$L_x = %d$, $T = %f$' %(Lx,T))
plt.xlabel(r'$h$')
plt.ylabel(r'$C_v$')
plt.savefig('Cv_vs_h_' + LDir + '_T' + str(T) + '.pdf')

plt.figure()
plt.plot(h_list, mAx_list, 'o-', label=r'$m_{Ax}$')
plt.plot(h_list, mAz_list, 'o-', label=r'$m_{Az}$')
plt.plot(h_list, msAx_list, 'o-', label=r'$m_{sAx}$')
plt.plot(h_list, msAz_list, 'o-', label=r'$m_{sAz}$')
plt.title(r'$L_x = %d$, $T = %f$' %(Lx,T))
plt.xlabel(r'$h$')
plt.ylabel(r'$m_A$')
plt.legend()
plt.savefig('mA_vs_h_' + LDir + '_T' + str(T) + '.pdf')

plt.figure()
plt.plot(h_list, mBy_list, 'o-', label=r'$m_{By}$')
plt.plot(h_list, mBz_list, 'o-', label=r'$m_{Bz}$')
plt.plot(h_list, msBy_list, 'o-', label=r'$m_{sBy}$')
plt.plot(h_list, msBz_list, 'o-', label=r'$m_{sBz}$')
#plt.plot(h_list, msB_list, 'o-', label=r'$m_{sB}$')
plt.title(r'$L_x = %d$, $T = %f$' %(Lx,T))
plt.xlabel(r'$h$')
plt.ylabel(r'$m_B$')
plt.legend()
plt.savefig('mB_vs_h_' + LDir + '_T' + str(T) + '.pdf')

#plt.figure()
#plt.plot(h_list, mAFpi_list, 'o-')
#plt.title(r'$L_x = %d$, $T = %f$' %(Lx,T))
#plt.xlabel(r'$h$')
#plt.ylabel(r'$m_{zs}$')
#
#plt.figure()
#plt.plot(h_list, mFAF_list, 'o-')
#plt.title(r'$L_x = %d$, $T = %f$' %(Lx,T))
#plt.xlabel(r'$h$')
#plt.ylabel(r'$m_{F-AF}$')

plt.show()
