#ifndef MIXED_POTENTIAL_H
#define MIXED_POTENTIAL_H
#include <string>
#include <vector>

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

#ifndef LJ_POTENTIAL_H
#include "lj_potential.h"
#endif

#ifndef SPRING_POTENTIAL_H
#include "spring_potential.h"
#endif

template <class T1, class T2>
class MixedPotential
{
    public:
        long double weight1,weight2,beta;
        Box box;
        T1 pot1;
        T2 pot2;
        MixedPotential(){};
        MixedPotential(T1& pot10, T2& pot20,
        long double& weight10,long double& weight20,long double& beta0);
		MixedPotential(const MixedPotential& other);
        long double conf_prob(long double& E);
        long double compute_potential(Conf& frame);
        long double compute_potential_change(Conf& frame, Atom& new_atom, int& atom_index);
        long double compute_prob_change(Conf& frame, Atom& new_atom, int& atom_index);

};
      
#endif 
