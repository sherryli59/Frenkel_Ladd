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

long double dUdgamma(LJPotential& pot, Conf& c, Box& box, vector<Atom>& lattice,Box& box0);
long double PairdUdgamma(LJPotential& pot, Atom& one,Atom& two, Box& box,
    Atom& one_ref, Atom& two_ref,Box& box0);
long double min_image(long double& x,long double& y,long double& boxlen);
int main(int argc, char const *argv[]){
    int na=256;
    int ncell=4;
    bool min_image = true;
    std::string lattice_type="FCC";
    long double rho0 = 1.28;
    long double gamma = atof(argv[1]);
    long double rho =rho0/pow(gamma,3);
    int gamma_index = atoi(argv[2]);
    string s_index = to_string(gamma_index);

    long double sigma = 1.0;
    long double epsilon = 1.0;
    long double cutoff = 2.7;
    long double beta = 1;
    long double alpha = 1264;


    
    int niterations = 100000;
    int nequil = 100000;

    int out_freq = 500;
	string s_alpha = to_string((int)alpha);
    string s_folder = argv[3];
	string dir_name = "output/"+s_folder+"/";
    std::string coord_out=dir_name+"alpha"+s_alpha+"coord"+s_index+".xyz";
    std::string coord_test="output/coord.xyz";
    std::string potential_out=dir_name+"alpha"+s_alpha+"dU"+s_index+".dat";
    std::string lj_out=dir_name+"alpha"+s_alpha+"lj"+s_index+".dat";
    std::string spring_out=dir_name+"alpha"+s_alpha+"spring"+s_index+".dat";
    std::ofstream coord_output(coord_out.c_str());
    std::ofstream test_output(coord_test.c_str());
	std::ofstream potential_output(potential_out.c_str());
	std::ofstream lj_output(lj_out.c_str());
	std::ofstream spring_output(spring_out.c_str());

 
    Box box(na,ncell,ncell,ncell,rho, min_image);
    Conf frame(na,box);
    frame.set_lattice(lattice_type);
    std::vector<Atom> lattice_site;
    lattice_site.reserve(na);
    for (int i = 0; i < na; ++i){
            lattice_site.push_back(frame.atoms[i]);
        } 
    
    Box box0(na,ncell,ncell,ncell,rho0, min_image);
    Conf init_frame(na,box0);
    init_frame.set_lattice(lattice_type);
    std::vector<Atom> init_lattice;
    output_coord(init_frame,coord_test);
    init_lattice.reserve(na);
    for (int i = 0; i < na; ++i){
            init_lattice.push_back(init_frame.atoms[i]);
        } 
    
    SpringPotential spring_pot(lattice_site,alpha, beta, box);
    LJPotential lj_pot(sigma, epsilon, cutoff, beta, box);
    long double default_weight=1;
    MixedPotential<LJPotential,SpringPotential> pot(lj_pot,spring_pot,default_weight,default_weight,beta);
    long double potential= pot.compute_potential(frame);
    frame.update_potential(potential);

    RNG rng(frame.natoms);
    MC<MixedPotential<LJPotential,SpringPotential>> mc(pot,beta);

    for (int n = 0; n<niterations+nequil; ++n){
        mc.sample_step(frame,rng);
        if( (n>nequil) && ( ((n-nequil)%out_freq)==0 ) ){
            output_coord(frame,coord_out);
            long double dU=dUdgamma(lj_pot,frame,box,init_lattice,box0);
            output_number(dU, potential_out);
            
            long double spring=spring_pot.compute_potential(frame);
            output_number(spring, spring_out);
            long double lj=lj_pot.compute_potential(frame);
            output_number(lj, lj_out);
        }
    }

    
    lattice_site.clear();
    return 0;
}

long double dUdgamma(LJPotential& pot, Conf& c, Box& box, vector<Atom>& lattice,Box& box0){
	
	long double E= 0.0;
	unsigned int natoms = c.natoms;
		for(unsigned int atom1 = 0;atom1<natoms-1;++atom1){

			for(unsigned int atom2 = atom1+1;atom2<natoms;++atom2){
				
				long double pairE= PairdUdgamma(pot, c.atoms[atom1], 
                c.atoms[atom2],box, lattice[atom1],lattice[atom2],box0);
				E += pairE;	
			}
		}
	return(E);
    }
    long double PairdUdgamma(LJPotential& pot, Atom& one, Atom& two, Box& box,
    Atom& one_ref, Atom& two_ref,Box& box0){
	
	    long double dU=pot.compute_pair_dUdr(one,two);
	    if (dU!=0){
        cout<<"dUdr:"<<dU<<endl;
        long double r=compute_distance(one,two,box);
        cout<<"final separation:"<<r<<endl;
        // cout<<one.x<<","<<one.y<<","<<one.z<<endl;
        // cout<<two.x<<","<<two.y<<","<<two.z<<endl;
        long double rvec[3];
        rvec[0]=min_image(one.x,two.x,box.lx);
        rvec[1]=min_image(one.y,two.y,box.ly);
        rvec[2]=min_image(one.z,two.z,box.lz);
        long double rvec_ref[3];
        rvec_ref[0]=min_image(one_ref.x,two_ref.x,box0.lx);
        rvec_ref[1]=min_image(one_ref.y,two_ref.y,box0.ly);
        rvec_ref[2]=min_image(one_ref.z,two_ref.z,box0.lz);
		long double drdgamma=1.0 / r *(rvec[0]*rvec_ref[0]
		+rvec[1]*rvec_ref[1]+rvec[2]*rvec_ref[2]);
        cout<<"drdgamma:"<<drdgamma<<endl;
        dU*=drdgamma;
	}
        long double r0=compute_distance(one_ref,two_ref,box);
        cout<<"initial separation:"<<r0<<endl;
	return(dU);
}

long double min_image(long double& x,long double& y,long double& boxlen){
    long double d=x-y;
    if (d>0.5*boxlen){
        d-=0.5*boxlen;
    }
    else if (d<-0.5*boxlen){
        d+=0.5*boxlen;
    }
    return d;
}
