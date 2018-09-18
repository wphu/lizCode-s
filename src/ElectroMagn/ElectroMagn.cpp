#include "ElectroMagn.h"

#include <limits>
#include <iostream>

#include "PicParams.h"
#include "Species.h"
#include "Projector.h"
#include "Field.h"
#include "Profile.h"
#include "SolverFactory.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(PicParams &params, InputData &input_data) :
timestep(params.timestep),
cell_length(params.cell_length),
n_species(params.species_param.size()),
nDim_field(params.nDim_field),
cell_volume(params.cell_volume),
n_space(params.n_space),
n_space_global(params.n_space_global),
oversize(params.oversize),
avg_step(params.ntime_step_avg),
dump_step(params.dump_step)
{
    species_param = params.species_param;

    // initialize poynting vector
    poynting[0].resize(nDim_field,0.0);
    poynting[1].resize(nDim_field,0.0);
    poynting_inst[0].resize(nDim_field,0.0);
    poynting_inst[1].resize(nDim_field,0.0);

    // initialize charge vector: the charge is scalar, not field
    emitCharge[0].resize(nDim_field, 0.0);
    emitCharge[1].resize(nDim_field, 0.0);
    depCharge[0].resize(nDim_field, 0.0);
    depCharge[1].resize(nDim_field, 0.0);
    totCharge[0].resize(nDim_field, 0.0);
    totCharge[1].resize(nDim_field, 0.0);



    // take useful things from params
    for (unsigned int i=0; i<3; i++) {
        DEBUG("____________________ OVERSIZE: " <<i << " " << oversize[i]);
    }

    if (n_space.size() != 3) ERROR("this should not happend");

    Ex_=NULL;
    Ey_=NULL;
    Ez_=NULL;
    Bx_=NULL;
    By_=NULL;
    Bz_=NULL;
    Bx_m=NULL;
    By_m=NULL;
    Bz_m=NULL;
    Jx_=NULL;
    Jy_=NULL;
    Jz_=NULL;
    rho_=NULL;

    Ex_avg=NULL;
    Ey_avg=NULL;
    Ez_avg=NULL;
    Bx_avg=NULL;
    By_avg=NULL;
    Bz_avg=NULL;

    // Species charge currents and density
    Jx_s.resize(n_species);
    Jy_s.resize(n_species);
    Jz_s.resize(n_species);
    rho_s.resize(n_species);
    rho_s_avg.resize(n_species);

    Vx_s.resize(n_species);
    Vx_s_avg.resize(n_species);


    Vy_s.resize(n_species);
    Vy_s_avg.resize(n_species);


    Vz_s.resize(n_species);
    Vz_s_avg.resize(n_species);


    Vp_s.resize(n_species);
    Vp_s_avg.resize(n_species);


    T_s.resize(n_species);
    T_s_avg.resize(n_species);

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]             = NULL;
        Jy_s[ispec]             = NULL;
        Jz_s[ispec]             = NULL;
        rho_s[ispec]            = NULL;
        rho_s_avg[ispec]        = NULL;

        Vx_s[ispec]             = NULL;
        Vx_s_avg[ispec]         = NULL;

        Vy_s[ispec]             = NULL;
        Vy_s_avg[ispec]         = NULL;

        Vz_s[ispec]             = NULL;
        Vz_s_avg[ispec]         = NULL;

        Vp_s[ispec]             = NULL;
        Vp_s_avg[ispec]         = NULL;

        T_s[ispec]              = NULL;
        T_s_avg[ispec]          = NULL;
    }

    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<2; j++) {
            istart[i][j]=0;
            bufsize[i][j]=0;
        }
    }

    //emBoundCond = ElectroMagnBC_Factory::create(params, laser_params);

    //MaxwellFaradaySolver_ = SolverFactory::create(params);

}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::~ElectroMagn()
{
    delete Ex_;
    delete Ey_;
    delete Ez_;
    delete Bx_;
    delete By_;
    delete Bz_;
    delete Bx_m;
    delete By_m;
    delete Bz_m;
    delete Jx_;
    delete Jy_;
    delete Jz_;
    delete rho_;

    if (Ex_avg!=NULL) {
        delete Ex_avg;
        delete Ey_avg;
        delete Ez_avg;
        delete Bx_avg;
        delete By_avg;
        delete Bz_avg;
    }

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
      delete Jx_s[ispec];
      delete Jy_s[ispec];
      delete Jz_s[ispec];
      delete rho_s[ispec];
    }


}//END Destructer



// ---------------------------------------------------------------------------------------------------------------------
// Method used to initialize the total charge density
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::initRhoJ(vector<Species*>& vecSpecies, Projector* Proj)
{
    //! \todo Check that one uses only none-test particles
    // number of (none-test) used in the simulation
    //! \todo fix this: n_species is already a member of electromagn, is it this confusing? what happens if n_species grows (i.e. with ionization)?
    unsigned int n_species = vecSpecies.size();

    //loop on all (none-test) Species
    for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ ) {
        Particles &cuParticles = vecSpecies[iSpec]->getParticlesList();
        unsigned int n_particles = vecSpecies[iSpec]->getNbrOfParticles();

        DEBUG(n_particles<<" species "<<iSpec);
	if (!cuParticles.isTestParticles) {
	    for (unsigned int iPart=0 ; iPart<n_particles; iPart++ ) {
		// project charge & current densities
		(*Proj)(Jx_s[iSpec], Jy_s[iSpec], Jz_s[iSpec], rho_s[iSpec], cuParticles, iPart,
			cuParticles.lor_fac(iPart));
	    }
	}

    }//iSpec
    DEBUG("before computeTotalRhoJ");
    computeTotalRhoJ();
    DEBUG("projection done for initRhoJ");

}

