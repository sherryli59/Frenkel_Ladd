#include <iostream>
#ifndef CONF_H
#include "../include/conf.h"
#endif

#ifndef UTILS_H
#include "../include/utils.h"
#endif

#ifndef LJ_POTENTIAL_H
#include "../include/lj_potential.h"
#endif

#ifndef SPRING_POTENTIAL_H
#include "../include/spring_potential.h"
#endif

#ifndef MC_H
#include "../include/mc.h"
#endif

template <class T>
MC<T>::MC(T& potential0, long double& beta0){
    stepsize=0.05;
    potential=potential0;
    beta=beta0;
    nsteps=0;
    acc_rate=0;
}

template <class T>
void MC<T>::sample_step(Conf& frame, RNG& rng){
    int atom_index =rng.random_atom(frame.natoms);
    long double dx,dy,dz;
    double randn[4];
    for (int i=0; i<4; i++)
        randn[i]=rng.rand_01();
    dx = stepsize*(randn[0]*2-1);
    dy = stepsize*(randn[1]*2-1);
    dz = stepsize*(randn[2]*2-1);
    Atom proposal=frame.atoms[atom_index];
    proposal.displace(dx,dy,dz);
    long double acc_prob=potential.compute_prob_change(frame, proposal, atom_index);
    if (randn[3]<acc_prob){
        frame.update_frame(atom_index,dx,dy,dz);
        long double new_potential = frame.potential-log(acc_prob)/beta;
        frame.update_potential(new_potential);
        acc_rate= (acc_rate*nsteps+1)/(nsteps+1);
    }
    else {
        acc_rate= (acc_rate*nsteps)/(nsteps+1);
    }
    nsteps+=1;
    //if (nsteps%1000==0){adjust_stepsize();}
};

template <class T>
void MC<T>::adjust_stepsize(){
if (acc_rate<0.2){
    stepsize*=0.1;
}
else if (acc_rate<0.33){
    stepsize*=acc_rate/0.33;
}
else if (acc_rate>0.55){
    stepsize*=acc_rate/0.55;
}			
}