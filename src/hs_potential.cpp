#include <string>
#include <iostream>
#include <limits>

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


#ifndef HS_POTENTIAL_H
#include "../include/hs_potential.h"
#endif

#ifndef UTILS_H
#include "../include/utils_potential.h"
#endif


using std::cout;
using std::endl;

HSPotential::HSPotential(long double& radius0, Box& box0){
    radius=radius0;
    box=box0;
    inf=std::numeric_limits<long double>::infinity();
}
long double HSPotential::conf_prob(long double& E){
    if (E==0){
        return 1;
    }
    else return 0;
}        


long double HSPotential::pair_potential(Atom& atom1, Atom& atom2){      	
    long double r;
    r = compute_distance(atom1, atom2, box);
    if (r<2*radius){
        // cout<<r<<endl;
        // print_coord(atom1);
        // print_coord(atom2);
        return inf;
    }
    else return(0);
    }   

long double HSPotential::compute_potential(Conf& frame){
    long double E = 0.0;
    for(int atom1 = 0;atom1<frame.natoms-1;atom1++){
        for(int atom2 = atom1+1;atom2<frame.natoms;atom2++){
            E += pair_potential(frame.atoms[atom1], frame.atoms[atom2]);
            }
        }
return(E);
}


long double HSPotential::compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=0;
    for (int i = 0; i < frame.natoms; ++i){
        if(i!=atom_index){
            dE+= pair_potential(frame.atoms[i], new_atom)-pair_potential(frame.atoms[i], frame.atoms[atom_index]);}	
    }
    return dE;
}

long double HSPotential::compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=compute_potential_change(frame, new_atom, atom_index);
    return conf_prob(dE);}



