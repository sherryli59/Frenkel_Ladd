import numpy as np
x,w=np.polynomial.legendre.leggauss(30)
lamb=x/2+0.5
w=w/2
np.savetxt("lambda_leggauss.dat",lamb)
gamma=x+2
np.savetxt("gamma_leggauss.dat",gamma)
np.savetxt("weight_leggauss.dat",w)
