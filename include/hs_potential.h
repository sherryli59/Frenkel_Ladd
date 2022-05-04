#ifndef HS_POTENTIAL_H
#define HS_POTENTIAL_H
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

#ifndef UTILS_H
#include "../include/utils_potential.h"
#endif


class HSPotential
{
    public:
        long double radius;
        Box box;
        long double inf;
        HSPotential(){};
        HSPotential(long double& radius, Box& box);
		HSPotential(const HSPotential& other);
        long double conf_prob(long double& E);
        long double compute_potential(Conf& frame);
        long double compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index);
        long double pair_potential(Atom& atom1, Atom& atom2);
        long double compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index);

};
      
#endif 
