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
#include "include/hs_potential.h"
#include "include/spring_potential.h"
#include "include/mixed_potential.h"
#include "src/mixed_potential.cpp"
#include "include/mc.h"
#include "src/mc.cpp"
#include "include/utils.h"


int main(int argc, char const *argv[]){
    int na=256;
    int ncell=4;
    long double rho=1.04;
    bool min_image = true;
    std::string lattice_type="FCC";

    long double radius = 0.5;
    long double beta = 1;
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
    std::string coord_test="output/coord.xyz";
    std::string potential_out=dir_name+"alpha"+s_alpha+"dU"+s_index+".dat";
    std::ofstream coord_output(coord_out.c_str());
    std::ofstream test_output(coord_test.c_str());
	std::ofstream potential_output(potential_out.c_str());

    Box box(na,ncell,ncell,ncell,rho, min_image);
    Conf frame(na,box);
    cout<<box.lcellx<<endl;
    frame.set_lattice(lattice_type);
    std::vector<Atom> lattice_site;
    lattice_site.reserve(na);
    for (int i = 0; i < na; ++i){
            lattice_site.push_back(frame.atoms[i]);
        } 
    
    SpringPotential spring_pot(lattice_site,alpha, beta, box);
    HSPotential hs_pot(radius, box);
    long double default_weight=1;
    MixedPotential<HSPotential,SpringPotential> pot(hs_pot,spring_pot,default_weight,lambda,beta);
    
    output_coord(frame,coord_test);
    RNG rng(frame.natoms);

    MC<MixedPotential<HSPotential,SpringPotential>> mc(pot,beta);
    long double potential= pot.compute_potential(frame);
    frame.update_potential(potential);
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


