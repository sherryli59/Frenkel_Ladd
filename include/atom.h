#ifndef ATOM_H
#define ATOM_H

class Atom
{
    public:
        long double x, y, z;
		long double mass;

        Atom();

        Atom(long double& x_coord, long double& y_coord, long double& z_coord);
		Atom(const Atom& other);
		
		void displace(long double& dx, long double& dy, long double& dz);
};
#endif