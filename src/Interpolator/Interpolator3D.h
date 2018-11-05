#ifndef INTERPOLATOR3D_H
#define INTERPOLATOR3D_H

#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 3D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D : public Interpolator
{
public:
    Interpolator3D(PicParams&params): Interpolator(params) {};

    virtual ~Interpolator3D() {};

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    double dy_inv_;
    double dz_inv_;
};

#endif
