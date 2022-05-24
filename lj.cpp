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
    int na=32;
    int ncell=2;
    long double rho=1.28;
    bool min_image = true;
    std::string lattice_type="FCC";

    long double sigma = 1.0;
    long double epsilon = 1.0;
    long double cutoff = 2.7;
    long double beta = 0.5;



    
    int niterations = 5000000;
    int nequil = 100000;

    int out_freq = 500;
    string s_folder = argv[1];
	string dir_name = "output/"+s_folder+"/";
    std::string coord_out=dir_name+"lj.xyz";
    std::string potential_out=dir_name+"U_lj.dat";
    std::string init_coord="output/ref.xyz";
    std::ofstream coord_output(coord_out.c_str());
	std::ofstream potential_output(potential_out.c_str());

    Box box(na,ncell,ncell,ncell,rho, min_image);
    Conf frame(na,box);
    frame.set_lattice(lattice_type);

    output_coord(frame,init_coord);
    
    LJPotential pot(sigma, epsilon, cutoff, beta, box);
    long double potential= pot.compute_potential(frame);
    frame.update_potential(potential);

    RNG rng(frame.natoms);
    MC<LJPotential> mc(pot,beta);

    for (int n = 0; n<niterations+nequil; ++n){
        mc.sample_step(frame,rng);
        if( (n>nequil) && ( ((n-nequil)%out_freq)==0 ) ){
            frame.recenter();
            output_coord(frame,coord_out);
            long double U=pot.compute_potential(frame);
            output_number(U, potential_out);
        }
    }
    cout<<"final rate:"<<mc.acc_rate<<endl;
    cout<<"final stepsize:"<<mc.stepsize<<endl;
    return 0;
}


