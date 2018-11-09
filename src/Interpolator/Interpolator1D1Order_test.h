#ifndef INTERPOLATOR1D1ORDER_TEST_H
#define INTERPOLATOR1D1ORDER_TEST_H

#include "Field1D.h"
#include "ElectroMagn.h"
#include "Particles.h"
//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D1Order_test
{

public:
    Interpolator1D1Order_test(PicParams&);
    ~Interpolator1D1Order_test(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);

    inline double compute( double* coeff, Field1D* f, int idx) {
    	double interp_res =  coeff[0] * (*f)(idx)   + coeff[1] * (*f)(idx+1);
    	return interp_res;
    };

private:
    double dx_inv_;
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Interpolation coefficient on Prim grid
    double coeffp_[2];
    // Interpolation coefficient on Dual grid
    double coeffd_[2];


};//END class

#endif
