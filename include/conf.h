#ifndef CONF_H
#define CONF_H

#include <string.h>
#include <vector>
#include <iostream>


#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef BOX_H
#include "box.h"
#endif

using std::vector;

class Conf
{
    public:
        int natoms;
        Atom * atoms;
        vector<long double> com;
        Box box;
        long double mass;
        long double potential;

        Conf();

        Conf(int& na, const Box& box0);

        ~Conf();

        Conf(const Conf& other);

        void set_lattice(const std::string& type);
        void update_mass();
        vector<long double> compute_com();
        void update_com();
        void update_frame(int& atom_index, long double& dx, long double& dy, long double& dz);
        void shift_com(int& atom_index, long double& dx, long double& dy, long double& dz);
        void update_potential(long double& new_potential);
        void wrap_atom_coordinates(Atom& atom, Box& box);
};
#endif
