#ifndef MC_H
#define MC_H

#ifndef CONF_H
#include "conf.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef LJ_POTENTIAL_H
#include "lj_potential.h"
#endif

#ifndef SPRING_POTENTIAL_H
#include "spring_potential.h"
#endif

#include <iostream>
template <class T>
class MC{
    public: 
        long double stepsize,beta, acc_rate;
        int nsteps;
        T potential;
        MC(T& potential0, long double& beta0);
        void sample_step(Conf& frame, RNG& rng);
        void adjust_stepsize();
};
#endif