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

print(D_dip)
