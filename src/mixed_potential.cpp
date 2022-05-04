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

#ifndef HS_POTENTIAL_H
#include "../include/hs_potential.h"
#endif

#ifndef SPRING_POTENTIAL_H
#include "../include/spring_potential.h"
#endif

#ifndef MIXED_POTENTIAL_H
#include "../include/mixed_potential.h"
#endif

using std::cout;
using std::endl;

template <class T1, class T2>
MixedPotential<T1,T2>::MixedPotential(T1& pot10, T2& pot20,
        long double& weight10,long double& weight20,long double& beta0)
    {
        pot1=pot10;
        pot2=pot20;
        weight1=weight10;
        weight2=weight20;
        beta = beta0;
    }
template <class T1, class T2>
MixedPotential<T1,T2>::MixedPotential(const MixedPotential& other){
    pot1=other.pot1;
    pot2=other.pot2;
    weight1=other.weight1;
    weight2=other.weight2;
    beta=other.beta;
};
template <class T1, class T2>
long double MixedPotential<T1,T2>::conf_prob(long double& E){
    return exp(-beta*E);}  

template <class T1, class T2>
long double MixedPotential<T1,T2>::compute_potential(Conf& frame){
    long double E=0;
    E += weight1 * pot1.compute_potential(frame);
    E += weight2 * pot2.compute_potential(frame);
    return E;
}

template <class T1, class T2>
long double MixedPotential<T1,T2>::compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=0;
    dE += weight1 * pot1.compute_potential_change(frame,new_atom, atom_index);
    dE += weight2 * pot2.compute_potential_change(frame, new_atom, atom_index);
    return dE;
}

template <class T1, class T2>
long double MixedPotential<T1,T2>::compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=compute_potential_change(frame, new_atom, atom_index);
    return conf_prob(dE);}



