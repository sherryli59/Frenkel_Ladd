#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <time.h>
#include <sstream>
#include <stdexcept>
#include<vector>


using namespace std;

#include "include/conf.h"
#include "include/box.h"
#include "include/lj_potential.h"
#include "include/spring_potential.h"
#include "include/mixed_potential.h"
#include "src/mixed_potential.cpp"
#include "include/mc.h"
#include "src/mc.cpp"
#include "include/utils.h"


int main(int argc, char const *argv[]){
    int na=256;
    int ncell=4;
    long double rho=1.28;
    bool min_image = true;
    std::string lattice_type="FCC";

    long double sigma = 1.0;
    long double epsilon = 1.0;
    long double cutoff = 2.7;
    long double beta = 0.5;
    long double alpha = 1264.0;
    long double lambda = atof(argv[1]);
    int lambda_index = atoi(argv[2]);
    string s_index = to_string(lambda_index);

    
    int niterations = 100000;
    int nequil = 100000;

    int out_freq = 500;
	string s_alpha = to_string((int)alpha);
    string s_folder = argv[3];
	string dir_name = "output/"+s_folder+"/";
    std::string coord_out=dir_name+"alpha"+s_alpha+"coord"+s_index+".xyz";
    std::string init_coord="output/coord.xyz";
    std::string potential_out=dir_name+"alpha"+s_alpha+"dU"+s_index+".dat";
    std::ofstream coord_output(coord_out.c_str());
    std::ofstream test_output(init_coord.c_str());
	std::ofstream potential_output(potential_out.c_str());

    Box box(na,ncell,ncell,ncell,rho, min_image);
    Conf frame(na,box);
    frame.set_lattice(lattice_type);
    output_coord(frame,init_coord);
    std::vector<Atom> lattice_site;
    lattice_site.reserve(na);
    for (int i = 0; i < na; ++i){
            lattice_site.push_back(frame.atoms[i]);
        } 
    
    SpringPotential spring_pot(lattice_site,alpha, beta, box);
    LJPotential lj_pot(sigma, epsilon, cutoff, beta, box);
    long double default_weight=1;
    MixedPotential<LJPotential,SpringPotential> pot(lj_pot,spring_pot,default_weight,lambda,beta);
    long double potential= pot.compute_potential(frame);
    frame.update_potential(potential);

    RNG rng(frame.natoms);
    MC<MixedPotential<LJPotential,SpringPotential>> mc(pot,beta);

    for (int n = 0; n<niterations+nequil; ++n){
        mc.sample_step(frame,rng);
        if( (n>nequil) && ( ((n-nequil)%out_freq)==0 ) ){
            output_coord(frame,coord_out);
            long double dU=spring_pot.compute_potential(frame);
            output_number(dU, potential_out);
        }
    }
    cout<<"final rate:"<<mc.acc_rate<<endl;
    cout<<"final stepsize:"<<mc.stepsize<<endl;
    lattice_site.clear();
    return 0;
}


