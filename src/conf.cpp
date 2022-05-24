#include <string.h>
#include <iostream>
#include <vector>

#include "../include/utils.h"
#include "../include/atom.h"
#include "../include/box.h"
#include "../include/conf.h"

using std::vector;

Conf::Conf(){
    atoms = new Atom[natoms];
} 

Conf::Conf(int& na, const Box& box0){
    natoms = na;
    atoms = new Atom[natoms];
    box = box0;

}

Conf::~Conf()
{
    delete [] atoms;
}

Conf::Conf(const Conf& other){
    natoms = other.natoms;
    box = other.box;
    com = other.com;
    mass = other.mass;
    atoms = new Atom[natoms];
    for (int i = 0; i < natoms; ++i){
            atoms[i] = other.atoms[i];
        }
}

void Conf::set_lattice(const std::string& type){
    int nbasis;
    std::vector< std::vector<double> > basis ;
    if (type.compare("FCC")==0) {
        nbasis=4;
        basis={{0,0,0},{0,0.5,0.5},{0.5,0,0.5},{0.5,0.5,0}}; 
    }
    for (int x = 0; x < box.ncellx; ++x)
    {
        for (int y = 0; y < box.ncelly; ++y)
        {
            for (int z = 0; z < box.ncellz; ++z)
            {
                int cell_index = z*box.ncellx*box.ncelly+y*box.ncellx+x;
                for (int k = 0; k < nbasis; ++k)
                {
                    int atom_index = cell_index * 4 + k;
                    atoms[atom_index].x = (x+basis[k][0]+0.25) * box.lcellx - box.lx/2;
                    atoms[atom_index].y = (y+basis[k][1]+0.25) * box.lcelly - box.ly/2;
                    atoms[atom_index].z = (z+basis[k][2]+0.25) * box.lcellz - box.lz/2;
                    atoms[atom_index].mass = 1;
                }
            }
        }
    }
    natoms = box.ncellx*box.ncelly*box.ncellz*nbasis;
    update_mass();
    update_com();
}
vector<long double> Conf::compute_com(){
    long double M = 0.0;
    long double sumrx=0.0, sumry=0.0, sumrz=0.0;
    for(unsigned int i = 0;i<natoms;i++){
        
        long double mass = atoms[i].mass;
        sumrx += atoms[i].x*mass;
        sumry += atoms[i].y*mass;
        sumrz += atoms[i].z*mass;
        M += mass;	
    }
    vector<long double> com= {sumrx/M,sumry/M,sumrz/M};
    return com;
}
void Conf::update_mass(){			
    mass = 0.0;			
    for(unsigned int i = 0;i<natoms;++i){	
        mass += atoms[i].mass;
    }
    return;
} 
void Conf::update_com(){			
    com=compute_com();
    return;
}
void Conf::update_frame(int& atom_index, long double& dx, long double& dy, long double& dz){
    atoms[atom_index].x+=dx;
    atoms[atom_index].y+=dy;
    atoms[atom_index].z+=dz;
    shift_com(atom_index, dx, dy,dz);
    if (box.min_image){
        wrap_atom_coordinates(atoms[atom_index],box);
    }
}
void Conf::shift_com(int& atom_index, long double& dx, long double& dy, long double& dz){
    com[0]+=dx*atoms[atom_index].mass/mass;
    com[1]+=dy*atoms[atom_index].mass/mass;
    com[2]+=dz*atoms[atom_index].mass/mass;
}
void Conf::update_potential(long double& new_potential){
    potential=new_potential;
}

void Conf::wrap_atom_coordinates(Atom& atom, Box& box){
        long double x,y,z;
        if (atom.x-com[0]<box.x0){
            atom.x+=box.lx;
        }
        else if (atom.x-com[0]>box.x1){
            atom.x-=box.lx;
        }
        if (atom.y-com[1]<box.y0){
            atom.y+=box.ly;
        }
        else if (atom.y-com[1]>box.y1){
            atom.y-=box.ly;
        }
        if (atom.z-com[2]<box.z0){
            atom.z+=box.lz;
        }
        else if (atom.z-com[2]>box.z1){
            atom.z-=box.lz;
        }
    }
void Conf::recenter(){
    for (unsigned int i = 0; i < natoms; ++i)
		{
			atoms[i].x -= com[0];
			atoms[i].y -= com[1];
			atoms[i].z -= com[2];
		}
    for (int i=0; i<3; ++i){
        com[i]=0;
    }
}
