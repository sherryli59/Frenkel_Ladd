#include <iostream>
#include <random>
#include <string>
#include<bits/stdc++.h>

#ifndef CONF_H
#include "../include/conf.h"
#endif


#ifndef ATOM_H
#include "../include/atom.h"
#endif

#ifndef BOX_H
#include "../include/box.h"
#endif

#ifndef LJ_POTENTIAL_H
#include "../include/lj_potential.h"
#endif

#ifndef UTILS_H
#include "../include/utils.h"
#endif

using namespace std;

RNG::RNG():max_num(1), rng(rd()), randreal(std::uniform_real_distribution<>(0, 1)){
rng.seed(::time(NULL));
}
RNG::RNG(int& natoms) : max_num(natoms-1), rng(rd()), randint(std::uniform_int_distribution<>(0, natoms-1)),randreal(std::uniform_real_distribution<>(0.0, 1.0)){
rng.seed(::time(NULL));
}
int RNG::random_atom() {
    return randint(rng);
}
int RNG::random_atom(int& natoms) {
    std::uniform_int_distribution<int> distr(0, natoms-1);
    return distr(rng);
}
double RNG::rand_01(){
    return randreal(rng);
}


long double compute_distance(Atom& atom1, Atom&atom2, Box& box){
	long double dx, dy, dz, r, rs;
	long double x1, y1, z1, x2, y2, z2;
	x1 = atom1.x, y1=atom1.y, z1=atom1.z;
	x2 = atom2.x, y2=atom2.y, z2=atom2.z;
	dx = x1-x2;	
	dy = y1-y2;
	dz = z1-z2;
    if (box.min_image){
        if(abs(dx)> box.lx/2.0){
		    dx = box.lx - abs(dx);
	    }
	    if(abs(dy)>box.ly/2.0){
		    dy = box.ly - abs(dy);
	    }
	    if(abs(dz)>box.lz/2.0){
		    dz = box.lz - abs(dz);
	    }   
    }
	rs= (pow(dx, 2))+(pow(dy, 2))+(pow(dz, 2));
	r = sqrt(rs);
    return r;
}

void output_coord(Conf& frame, std::string& coutname){
	std::ofstream coutfile(coutname.c_str(),ios::app);
	coutfile << frame.natoms << endl;
	coutfile <<"U: "<<frame.potential<<endl;
	coutfile << setiosflags(ios::fixed) << setprecision(10);
	for(int i = 0; i<frame.natoms; i++){
		coutfile<<"1 "<<frame.atoms[i].x<<" "<<frame.atoms[i].y<<" "<<frame.atoms[i].z<<endl; 
	}
}

void output_number(long double& number, std::string& coutname){
	std::ofstream coutfile(coutname.c_str(),ios::app);
	coutfile << number << endl;
}

void print_coord(Atom& atom){
	cout<<"coordinate:("<<atom.x<<","<<atom.y<<","<<atom.z<<")"<<endl;
}
