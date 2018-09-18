// PSI class for Reflection and Deposition

#ifndef PSI2D_RELDEP_H
#define PSI2D_RELDEP_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "PSI2D.h"
#include "Backscatterin_EmpiricalFormula.h"

class PSI2D_RefDep : public PSI2D
{

public:
    //! Constructor
    PSI2D_RefDep(        
        PicParams& params,
        vector<Species*>& vecSpecies,
        unsigned int psi_species1,
        unsigned int psi_species2,
        bool psi_is_self_consistent,
        string psiPosition,
        double emitTemperature);
    ~PSI2D_RefDep();

    //! Method called in the main smilei loop to apply PSI at each timestep
    void performPSI(PicParams& params, Grid* grid, vector<Species*>& vecSpecies, ElectroMagn* fields, Diagnostic* diag, int itime);

    void init(std::vector<Species*>&);
    
    Backscattering *backscattering;
private:


};


#endif
