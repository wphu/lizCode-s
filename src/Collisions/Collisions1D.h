/*
Collisions1D class
*/

#ifndef COLLISIONS1D_H
#define COLLISIONS1D_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions.h"


class Collisions1D : public Collisions
{

public:
    //! Constructor for Collisions between two species
    Collisions1D(PicParams &params)
    : Collisions(params)
    {
        npairsRem.resize(params.n_space[0]);
        for(int i = 0; i < npairsRem.size(); i++)
        {
            npairsRem[i] = 0;
        }
    };
    virtual ~Collisions1D(){};
private:
    double dx_inv_;
    unsigned int index_domain_begin;

};


#endif
