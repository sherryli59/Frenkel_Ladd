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
k = array([1264.0])
rho0 = 1.28
nsteps = 30
beta = 0.5

################################################################################
# Lambda integration
################################################################################
lamb=loadtxt('../../lambda_leggauss.dat')
gamma=loadtxt('../../gamma_leggauss.dat')
weight=loadtxt('../../weight_leggauss.dat')
for alpha in k:
        dE_1=array([average(loadtxt('output/phase1/alpha%ddU%d.dat'%(int(alpha),i))) for i in arange(nsteps)])
        savetxt('thermo_integration_alpha%d.dat'%int(alpha), transpose([lamb, dE_1/N]),
                header='lambda dE', fmt='%4f %.4f')
        W1=np.sum(np.multiply(dE_1,weight))
        dE_2=array([average(loadtxt('output/phase2/alpha%ddU%d.dat'%(int(alpha),i))) for i in arange(nsteps)])
        savetxt('thermo_integration_alpha%d.dat'%int(alpha), transpose([gamma, dE_2/N]),
                header='lambda dE', fmt='%4f %.4f')
        W2=np.sum(np.multiply(dE_2,weight))
        W = W1+W2
        ################################################################################
        # Compute free energy.
        ################################################################################

        # Define harmonic reference system free energy 
        F_harm = -3/2*N* log(2*pi/alpha) # [kT].

        # Fixed center of mass correction [Eq.(24) in the paper].
        F_CM= -1.5*log(alpha*beta/(2*pi))+log(rho0)-1.5*log(N) #[kT]

        # Compute absolute free energy per atom [Eq.(16) in the paper] and save data.
        F = (F_harm - W + F_CM) / N # [kT/atom].
        print (F_harm / N, W1/N, W2/N, F_CM/N)
        print (average(F),std(F))
        savetxt('free_energy_%d.dat'%alpha, transpose([F, F_harm/N, W/N, F_CM/N]),
                header='F[kT/atom]  F_harm  W F_CM', fmt='%4f')

        ################################################################################
        # Plot cumulative free energy difference
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(lamb,dE_1/N,label="phase1")
        ax.set_xlabel(r'$\lambda$')
        ax.legend()
        ax.set_ylabel(r'Thermodynamic integration [kT]')
        fig.savefig("FEdiff_I%d.png"%alpha,dpi=300)
        plt.show()
        plt.close()

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(gamma,dE_2/N,label="phase2")
        ax.set_xlabel(r'$\gamma$')
        ax.legend()
        ax.set_ylabel(r'Thermodynamic integration [kT]')
        fig.savefig("FEdiff_II%d.png"%alpha,dpi=300)
        plt.show()
        plt.close()
        ################################################################################
