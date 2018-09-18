/*
PSI1D class
*/

#ifndef PSI1D_H
#define PSI1D_H

#include <vector>
#include "PSI.h"
#include "Diagnostic1D.h"

class PSI1D : public PSI
{

public:
    //! Constructor for PSI between two species
    PSI1D(PicParams& params) : PSI(params)
    {
        nPartEmit_rem = 0.0;
    };
    virtual ~PSI1D(){};

    // sputtered particle number
    int nPartEmit;
    double nPartEmit_rem;

    //! Method called in the main smilei loop to apply collisions at each timestep
    virtual void performPSI(PicParams& params, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime){};

private:

};


#endif
