#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H
#include <string>

#ifndef CONF_H
#include "conf.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef BOX_H
#include "box.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif



class LJPotential
{
    public:
        long double sigma, epsilon,cutoff,beta;
        Box box;
        LJPotential(){};
        LJPotential(long double& sigma0, long double& epsilon0,long double& cutoff0);
        LJPotential(long double& sigma0, long double& epsilon0,long double& cutoff0, long double& beta0,Box& box0);
		LJPotential(const LJPotential& other);
        long double conf_prob(long double& E);
        long double compute_potential(Conf& frame);
        long double compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index);
        long double pair_potential(Atom& atom1, Atom& atom2);
        long double compute_pair_dUdr(Atom& atom1, Atom& atom2);
        long double compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index);

};
      
#endif 
