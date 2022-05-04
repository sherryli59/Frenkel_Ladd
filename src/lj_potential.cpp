#include <string>
#include <iostream>


#ifndef CONF_H
#include "../include/conf.h"
#endif

#ifndef UTILS_H
#include "../include/utils.h"
#endif

#ifndef BOX_H
#include "../include/box.h"
#endif

#ifndef ATOM_H
#include "../include/atom.h"
#endif


#ifndef LJ_POTENTIAL_H
#include "../include/lj_potential.h"
#endif

using std::cout;
using std::endl;

LJPotential::LJPotential(long double& sigma0, long double& epsilon0,long double& cutoff0){
    sigma=sigma0;
    epsilon=epsilon0;
    cutoff=cutoff0;
}
LJPotential::LJPotential(long double& sigma0, long double& epsilon0,long double& cutoff0, long double& beta0, Box& box0){
    sigma=sigma0;
    epsilon=epsilon0;
    cutoff=cutoff0;
    beta=beta0;
    box=box0;
}
LJPotential::LJPotential(const LJPotential& other){
    sigma = other.sigma;
    epsilon = other. epsilon;
    cutoff = other.cutoff;
    box = other.box;
    beta = other.beta;
}
long double LJPotential::conf_prob(long double& E){
    return exp(-beta*E);}       
        


long double LJPotential::pair_potential(Atom& atom1, Atom& atom2){      	
    long double r, pow_6, pow_6_shift;
    long double potent;
    r = compute_distance(atom1, atom2, box);
    pow_6 = pow(sigma/r,6);
    pow_6_shift = pow(sigma/cutoff,6);
    if (r<cutoff){
        potent = 4.0*epsilon*((pow(pow_6, 2))-pow_6)- 4.0*epsilon*((pow(pow_6_shift, 2))-pow_6_shift); 
    }
    else{
        potent=0.0;
    }	
    return(potent);
    }   

long double LJPotential::compute_potential(Conf& frame){
    long double E = 0.0;
    for(int atom1 = 0;atom1<frame.natoms-1;atom1++){
        for(int atom2 = atom1+1;atom2<frame.natoms;atom2++){
            if (compute_distance(frame.atoms[atom1], frame.atoms[atom2],box)==0){
                cout<<"atom1:"<<atom1<<"atom2:"<<atom2<<endl;
            }
            E += pair_potential(frame.atoms[atom1], frame.atoms[atom2]);
            }
        }
return(E);
}

long double LJPotential::compute_pair_dUdr(Atom& atom1, Atom& atom2){
    long double r, pow_6, pow_6_shift;
    r = compute_distance(atom1, atom2, box);
    pow_6 = pow(sigma/r,6);
    pow_6_shift = pow(sigma/cutoff,6);
    if (r<cutoff){
		return -24.0*epsilon/r*(2*(pow(pow_6, 2))-pow_6)+ 24.0*epsilon/cutoff*(2*(pow(pow_6_shift, 2))-pow_6_shift); //Currently I'm doing cut and shifted (should I?)
	}
	else{
		return 0.0;
	}	
}

long double LJPotential::compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=0;
    for (int i = 0; i < frame.natoms; ++i){
        if(i!=atom_index){
            dE+= pair_potential(frame.atoms[i], new_atom)-pair_potential(frame.atoms[i], frame.atoms[atom_index]);}	
    }
    return dE;
}

long double LJPotential::compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=compute_potential_change(frame, new_atom, atom_index);
    return conf_prob(dE);}



