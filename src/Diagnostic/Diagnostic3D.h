#ifndef DIAGNOSTIC3D_H
#define DIAGNOSTIC3D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "PSI3D.h"
#include "Grid3D.h"
#include "Particles.h"
#include "Field3D.h"


class Field;
class PSI;

class Diagnostic3D : public Diagnostic {

public :

    Diagnostic3D(PicParams& params, Grid* grid, ElectroMagn* EMfields, vector<PSI*>& vecPSI);
    virtual ~Diagnostic3D() {};

    //! Runs the diag for all patches for local diags.
    virtual void run( Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* EMfields, vector<PSI*>& vecPSI, int itime );

    int n_species;
    vector<unsigned int> dims_global;

    vector<Field3D*> particleFlux;
    vector<Field3D*> heatFlux;

    vector<Field3D*> particleFlux_global;
    vector<Field3D*> heatFlux_global;

    

    // calculate velocity and temperature of each species
    // not implemented
	void calVT(vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime){};

protected :
    double dx, dy, dz;
    double dx_inv_, dy_inv_, dz_inv_;
    int i_domain_begin, j_domain_begin, k_domain_begin;
    int dim_global;


};

#endif
