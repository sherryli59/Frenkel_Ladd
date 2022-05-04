#include "../include/atom.h"


Atom::Atom(){
	x = 0;
	y = 0;
	z = 0;
	mass = 0;
}

Atom::Atom(long double& x_coord, long double& y_coord, long double& z_coord): mass(1.0){
	x = x_coord;
	y = y_coord;
	z = z_coord;
}

Atom::Atom(const Atom& other){
	x=other.x;
	y=other.y;
	z=other.z;
	mass=other.mass;
}

void Atom::displace(long double& dx, long double& dy, long double& dz){
	x+=dx;
	y+=dy;
	z+=dz;
}

