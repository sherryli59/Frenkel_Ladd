#ifndef UTILS_H
#define UTILS_H

#ifndef CONF_H
#include "conf.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef BOX_H
#include "box.h"
#endif

#include <iostream>
#include <random>
#include<bits/stdc++.h>
#include <string>
#include <vector>

using namespace std;

class RNG
{
    private:
        std::random_device rd;
        int max_num;
        typedef std::mt19937 MyRng;
        MyRng rng;
        std::uniform_int_distribution<int> randint;
        std::uniform_real_distribution<double> randreal;
    public:
        RNG();
        RNG(int& natoms);
        int random_atom();
        int random_atom(int& natoms);
        double rand_01();
};



long double compute_distance(Atom& atom1, Atom&atom2, Box& box);

void output_coord(Conf& frame, std::string& coutname);
void output_number(long double& number, std::string& coutname);
void print_coord(Atom& atom);

#endif