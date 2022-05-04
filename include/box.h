#ifndef BOX_H
#define BOX_H
#include <stdlib.h>
#include <math.h>

class Box
{
    public:
        long double ncellx, ncelly, ncellz;
		long double lcellx, lcelly, lcellz;
        long double lx, ly, lz;
        long double x0, x1, y0, y1, z0, z1;
        bool min_image;

        Box();

        Box(int& natoms, int& ncellx0, int& ncelly0,int& ncellz0,long double& rho, bool& min_image0);

		Box(const Box& other);

};
#endif