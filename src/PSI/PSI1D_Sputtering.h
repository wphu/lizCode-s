#ifndef PSI1D_SPUTTERING_H
#define PSI1D_SPUTTERING_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI1D.h"
#include "PhysicalSputtering_EmpiricalFormula.h"

class PSI1D_Sputtering : public PSI1D
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Sputtering(
        PicParams& params,
        vector<Species*>& vecSpecies,
        unsigned int psi_species1,
        unsigned int psi_species2,
        bool psi_is_self_consistent,
        string psiPosition,
        double emitTemperature
    );

    ~PSI1D_Sputtering();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams& params, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime);

    // emit particles
    void emit(PicParams&, vector<Species*>&);

    void init(std::vector<Species*>&);

    Sputtering *sputtering;

private:


};


#endif
