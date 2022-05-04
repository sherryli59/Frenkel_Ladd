#ifndef SPRING_POTENTIAL_H
#define SPRING_POTENTIAL_H

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


using std::cout;
using std::endl;
using std::vector;
class SpringPotential
{
    public:
        long double alpha, beta;
        vector<Atom> eq_pos;
        Box box;
        bool com_adjusted;
        SpringPotential(){};
        SpringPotential(vector<Atom>& eq_pos0 , long double& alpha0);
        SpringPotential(vector<Atom>& eq_pos0, long double& alpha0,
        long double& beta0, Box& box0);
        SpringPotential(vector<Atom>& eq_pos0, long double& alpha0,
            long double& beta0, Box& box0, bool& com_adjusted0);
        SpringPotential(const SpringPotential& other);
        long double conf_prob(long double& E);
        long double atom_potential(Atom& atom, Atom& atom_eq, vector<long double>& com_adjust);
        long double atom_potential(Atom& atom, Atom& atom_eq);
        long double compute_potential(Conf& frame);

        long double compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index);

        long double compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index);


};
      
#endif 
