#include <stdlib.h>
#include <math.h>
#ifndef BOX_H
#include "../include/box.h"
#endif

Box::Box(){
    ncellx = 1;
    ncelly = 1;
    ncellz = 1;
    //Setting up cubic cells 
    lcellx = 1;
    lcelly = 1;
    lcellz = 1;
    lx = ncellx * lcellx;
    x0 = -0.5 * lx;
    x1 = 0.5 * lx;
    ly = ncelly * lcelly;
    y0 = -0.5 * ly;
    y1 = 0.5 * ly;
    lz = ncellz * lcellz;
    z0 = -0.5 * lz;
    z1 = 0.5 * lz;
    min_image = true;
}

Box::Box(int& natoms, int& ncellx0, int& ncelly0,int& ncellz0,long double& rho, bool& min_image0){
    ncellx = ncellx0;
    ncelly = ncelly0;
    ncellz = ncellz0;
    //Setting up cubic cells 
    lcellx = pow(natoms/(rho*ncellx*ncelly*ncellz),1.0/3);
    lcelly = lcellx;
    lcellz = lcellx;
    lx = ncellx * lcellx;
    x0 = -0.5 * lx;
    x1 = 0.5 * lx;
    ly = ncelly * lcelly;
    y0 = -0.5 * ly;
    y1 = 0.5 * ly;
    lz = ncellz * lcellz;
    z0 = -0.5 * lz;
    z1 = 0.5 * lz;
    min_image = min_image0;
}

Box::Box(const Box& other){
    ncellx = other.ncellx;
    ncelly = other.ncelly;
    ncellz = other.ncellz;
    lcellx = other.lcellx;
    lcelly = other.lcelly;
    lcellz = other.lcellz;
    lx = other.lx;
    ly = other.ly;
    lz = other.lz;
    x0 = other.x0;
    x1 = other.x1;
    y0 = other.y0;
    y1 = other.y1;
    z0 = other.z0;
    z1 = other.z1;
    min_image = other.min_image;
}

