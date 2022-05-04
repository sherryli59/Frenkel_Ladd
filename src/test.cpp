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

#include "../include/conf.h"
#include "../include/box.h"
#include "../include/mc.h"
#include "../include/potential.h"

int main(int argc, char const *argv[]){
    int na=256;
    int ncell=4;
    long double rho=1.28;
    bool min_image = true;
    Box box(na,ncell,ncell,ncell,rho, min_image);
    Conf frame(na,box);
    long double sigma = 1.0;
    long double epsilon = 1.0;
    long double cutoff = 2.7;
    long double beta = 0.5;
    LJPotential lj_pot(sigma, epsilon, cutoff, box);
    MC mc(lj_pot,beta);
    std::string type="FCC";
    frame.set_lattice(type);
    const int niterations = 10 ;
    RNG rng(frame.natoms);
    cout<<frame.natoms<<endl;
    for (int n = 0; n<niterations; n++){
        mc.sample_step(frame,rng);
        //cout<<frame.atoms[0].x<<endl;
        //cout<<lj_pot.compute_potential(frame)<<endl;
    }
    return 0;
}


