"""
  This script collects the data from forward and backward switchings and computed the absolute free energy for each temperature.

  Usage:
    python integrate.py
"""

from numpy import *
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

# Input parameters.
N = 256 
#alpha = 1264.12 
k = array([1264.052])
rho0 = 1.04
nsteps = 29


################################################################################
# Lambda integration
################################################################################
lamb=loadtxt('../../lambda_leggauss.dat')
weight=loadtxt('../../weight_leggauss.dat')
for alpha in k:
        dE_0=array([average(loadtxt('output/phase1/alpha%ddU%d.dat'%(int(alpha),i))) for i in arange(nsteps+1)])
        savetxt('thermo_integration_alpha%d.dat'%int(alpha), transpose([lamb, dE_0/N]),
                header='lambda dE', fmt='%4f %.4f')
        W=-np.sum(np.multiply(dE_0,weight))
        print(np.multiply(dE_0/N,weight))

        ################################################################################
        # Compute free energy.
        ################################################################################

        # Define harmonic reference system free energy 
        #F_harm = -3/2*(N-1)* log(4*pi**2/alpha) # [kT].
        F_harm = -3/2*(N-1)* log(2*pi/alpha) # [kT].

        # Fixed center of mass correction [Eq.(24) in the paper].
        F_CM= log(rho0*(2*pi/(N*alpha)**1.5)) #[kT]

        # Compute absolute free energy per atom [Eq.(16) in the paper] and save data.
        F = (F_harm + W + F_CM) / N # [kT/atom].
        print (F_harm / N, W/N ,F_CM/N)
        print (average(F),std(F))
        savetxt('free_energy_%d.dat'%alpha, transpose([F, F_harm/N, W/N, F_CM/N]),
                header='F[kT/atom]  F_harm  W F_CM', fmt='%4f')

        ################################################################################
        # Plot cumulative free energy difference
        fig = plt.figure()
        ax = fig.add_subplot()
        cumFE=array([F_harm/N+F_CM/N+trapz(dE_0[:n]/N,lamb[:n]) for n in range(len(dE_0))]) 
        ax.plot(lamb,cumFE,label="hard sphere")
        ax.set_xlabel(r'$\lambda$')
        ax.legend()
        ax.set_ylabel(r'Cumulative Free Energy [kT]')
        fig.savefig("FEdiff_I%d.png"%alpha,dpi=300)
        plt.show()
        plt.close()

        ################################################################################
