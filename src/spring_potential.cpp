#include <string>
#include <iostream>
#include <vector>

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

#ifndef SPRING_POTENTIAL_H
#include "../include/spring_potential.h"
#endif


using std::cout;
using std::endl;
using std::vector;

SpringPotential::SpringPotential(vector<Atom>& eq_pos0, long double& alpha0){
    eq_pos=eq_pos0;
    alpha = alpha0;
}
SpringPotential::SpringPotential(vector<Atom>& eq_pos0, long double& alpha0,
 long double& beta0, Box& box0){
    eq_pos=eq_pos0;
    alpha = alpha0;
    beta=beta0;
    box=box0;
    com_adjusted=true;
}
SpringPotential::SpringPotential(vector<Atom>& eq_pos0, long double& alpha0,
 long double& beta0, Box& box0, bool& com_adjusted0){
    eq_pos=eq_pos0;
    alpha = alpha0;
    beta=beta0;
    box=box0;
    com_adjusted = com_adjusted0;
}
SpringPotential::SpringPotential(const SpringPotential& other){
    alpha = other.alpha;
    eq_pos = other. eq_pos;
    box = other.box;
    beta = other.beta;
    com_adjusted = other.com_adjusted;
}
long double SpringPotential::conf_prob(long double& E){
    return exp(-beta*E);
}

long double SpringPotential::atom_potential(Atom& atom, Atom& atom_eq, vector<long double>& com_adjust){
    long double displacement[3] = {atom.x-atom_eq.x-com_adjust[0],
    atom.y-atom_eq.y-com_adjust[1],
    atom.z-atom_eq.z-com_adjust[2]};
	//cout<<"coordinate:("<<atom.x-com_adjust[0]<<","<<atom.y-com_adjust[1]<<","<<atom.z-com_adjust[2]<<")"<<endl;
    //print_coord(atom_eq);
	return(0.5 * alpha * (pow(displacement[0],2)+pow(displacement[1],2)+pow(displacement[2],2)));
}
long double SpringPotential::atom_potential(Atom& atom, Atom& atom_eq){
    long double displacement[3] = {atom.x-atom_eq.x,
    atom.y-atom_eq.y,
    atom.z-atom_eq.z};
	return(0.5 * alpha * (pow(displacement[0],2)+pow(displacement[1],2)+pow(displacement[2],2)));
}
long double SpringPotential::compute_potential(Conf& frame){
    long double E = 0.0;
    vector<long double> com_adjust;
    if (com_adjusted){ 
        com_adjust={frame.com[0],
            frame.com[1],frame.com[2]};
            //cout<<"com"<<frame.com[0]<<endl;
            }
       
    for(int i = 0; i<frame.natoms;i++){
        if (com_adjusted){ 
            E += atom_potential(frame.atoms[i],eq_pos[i],com_adjust);
            //cout<<E<<endl;
            }
        else {
            E += atom_potential(frame.atoms[i],eq_pos[i]);}
        }
    return(E);
}

long double SpringPotential::compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index){
    Atom old= frame.atoms[atom_index];
    float natoms = (float) frame.natoms;
    Atom eq = eq_pos[atom_index];
    if (com_adjusted){
        vector<long double> com_adjust={frame.com[0],
            frame.com[1],frame.com[2]};       
        long double displacement[3] = {new_atom.x- old.x, new_atom.y- old.y,new_atom.z- old.z};
        long double sq_distance = pow(displacement[0],2)+pow(displacement[1],2)+pow(displacement[2],2);
	    long double dRCM[3]= {com_adjust[0] + displacement[0] / natoms,
         com_adjust[1] + displacement[1] / natoms,com_adjust[2] + displacement[2] / natoms};
	    long double deltar[3] = {old.x- eq.x-com_adjust[0], old.y- eq.y-com_adjust[1],old.z-eq.z-com_adjust[2]};
	    return(0.5* alpha * (2*(displacement[0]*deltar[0]
        +displacement[1]*deltar[1]+displacement[2]*deltar[2])+(natoms-1)/natoms*sq_distance));
    }
    else {
        return atom_potential(new_atom,eq)-atom_potential(old,eq);
    }
}

long double SpringPotential::compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index){
    long double dE=compute_potential_change(frame, new_atom, atom_index);
    return conf_prob(dE);
}



