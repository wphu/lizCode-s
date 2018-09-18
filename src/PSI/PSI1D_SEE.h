
#ifndef PSI1D_SEE_H
#define PSI1D_SEE_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"

using namespace std;

class PSI1D_SEE : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_SEE(
        PicParams& params,
        unsigned int psi_species1,
        unsigned int psi_species2,
        bool psi_is_self_consistent,
        string psiPosition,
        double emitTemperature,
        double SEEYield
    );
    ~PSI1D_SEE();

    // secondary electron emission yield
    double SEEYield;

    Particles new_particles;

    //! Method called in the main smilei loop to apply collisions at each timestep
    void performPSI(PicParams& params, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime);

    void emit(PicParams& params, vector<Species*>& vecSpecies, unsigned int species_emit);
private:


};


#endif
