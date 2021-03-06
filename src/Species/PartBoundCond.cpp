#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond(PicParams& params, int ispec)
{
    // number of dimensions for the particle
    //!\todo (MG to JD) isn't it always 3?
    nDim_particle = params.nDim_particle;

    // Absolute global values
    double x_min_global = 0;
    double x_max_global = params.cell_length[0]*(params.n_space[0]);
    double y_min_global = 0;
    double y_max_global = params.cell_length[1]*(params.n_space[1]);
    double z_min_global = 0;
    double z_max_global = params.cell_length[2]*(params.n_space[2]);

    // by default apply no bcs
    bc_west   = NULL;
    bc_east   = NULL;
    bc_south  = NULL;
    bc_north  = NULL;
    bc_bottom = NULL;
    bc_up     = NULL;

    // -----------------------------
    // Define limits of local domain
    // -----------------------------

    // 1d3v or 2d3v or 3d3v
    x_min = 0.0;
    x_max = params.cell_length[0]*(params.n_space[0]);

    // 2d3v or 3d3v
    if ( nDim_particle > 1 ) {
        if ( (params.bc_em_type_y[0]=="periodic") || (params.bc_em_type_y[1]=="periodic") ) {
            y_min = 0.0;
            y_max = params.sim_length[1];
        }
        else {
            y_min = 0.0;
            y_max = params.sim_length[1];
        }
    }

    // 3d3v
    if ( nDim_particle > 2 ) {
        if ( (params.bc_em_type_z[0]=="periodic") || (params.bc_em_type_z[1]=="periodic") ) {
            z_min = 0.0;
            z_max = params.sim_length[2];
        }
        else {
            z_min = 0.0;
            z_max = params.sim_length[2];
        }
    }

    // Check for inconsistencies between EM and particle BCs
    if ( ((params.bc_em_type_x[0]=="periodic")&&(params.species_param[ispec].bc_part_type_west!="periodic"))
     ||  ((params.bc_em_type_x[1]=="periodic")&&(params.species_param[ispec].bc_part_type_east!="periodic")) ) {
        WARNING("For species #" << ispec << ", periodic EM boundary conditions require x particle BCs to be periodic.");
    }
    if ( nDim_particle > 1 ) {
        if ( ((params.bc_em_type_y[0]=="periodic")&&(params.species_param[ispec].bc_part_type_south!="periodic"))
         ||  ((params.bc_em_type_y[1]=="periodic")&&(params.species_param[ispec].bc_part_type_north!="periodic")) ) {
            WARNING("For species #" << ispec << ", periodic EM boundary conditions require y particle BCs to be periodic.");
        }
        if ( nDim_particle > 2 ) {
            if ( ((params.bc_em_type_z[0]=="periodic")&&(params.species_param[ispec].bc_part_type_bottom!="periodic"))
             ||  ((params.bc_em_type_z[1]=="periodic")&&(params.species_param[ispec].bc_part_type_up!="periodic"    )) ) {
                WARNING("For species #" << ispec << ", periodic EM boundary conditions require z particle BCs to be periodic.");
            }
        }
    }

    // ----------------------------------------------
    // Define the kind of applied boundary conditions
    // ----------------------------------------------

    bool thermCond = false;

    // West
    if ( params.species_param[ispec].bc_part_type_west == "refl" ) {
        bc_west = &refl_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "supp" ) {
        bc_west = &supp_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "stop" ) {
        bc_west = &stop_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "periodic" ) {
        bc_west = &periodic_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "thermalize" ) {
        thermCond = true;
        bc_west = &thermalize_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "none" ) {
        WARNING( "West boundary condition for species " << ispec << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "West boundary condition undefined" );
    }

    // East
    if ( params.species_param[ispec].bc_part_type_east == "refl" ) {
        bc_east = &refl_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "supp" ) {
        bc_east = &supp_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "stop" ) {
        bc_east = &stop_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "periodic" ) {
        bc_east = &periodic_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "thermalize" ) {
        thermCond = true;
        bc_east = &thermalize_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "none" ) {
        WARNING( "East boundary condition for species " << ispec << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "East boundary condition undefined" );
    }


    if ( nDim_particle > 1 ) {
        // South
        if ( params.species_param[ispec].bc_part_type_south == "refl" ) {
            bc_south = &refl_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "supp" ) {
            bc_south = &supp_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "stop" ) {
            bc_south = &stop_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "periodic" ) {
            bc_south = &periodic_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "thermalize" ) {
            thermCond = true;
            bc_south = &thermalize_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "none" ) {
            WARNING( "South boundary condition for species " << ispec << " is 'none', which means the same as fields");
        }
        else {
            ERROR( "South boundary condition undefined : " << params.species_param[ispec].bc_part_type_south  );
        }

        // North
        if ( params.species_param[ispec].bc_part_type_north == "refl" ) {
            bc_north = &refl_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "supp" ) {
            bc_north = &supp_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "stop" ) {
            bc_north = &stop_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "periodic" ) {
            bc_north = &periodic_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "thermalize" ) {
            thermCond = true;
            bc_north = &thermalize_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "none" ) {
            WARNING( "North boundary condition for species " << ispec << " is 'none', which means the same as fields");
        }
        else {
            ERROR( "North boundary condition undefined : " << params.species_param[ispec].bc_part_type_north  );
        }


        if ( nDim_particle > 2 ) {
            if ( params.species_param[ispec].bc_part_type_bottom == "refl" ) {
                if (z_min==z_min_global) bc_bottom = &refl_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_bottom == "supp" ) {
                if (z_min==z_min_global) bc_bottom = &supp_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_bottom == "stop" ) {
                if (z_min==z_min_global) bc_bottom = &stop_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_bottom == "periodic" ) {
                if (z_min==z_min_global) bc_bottom = &periodic_particle;
            }

            if ( params.species_param[ispec].bc_part_type_up == "refl" ) {
                if (z_min==z_min_global) bc_up = &refl_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_up == "supp" )  {
                if (z_min==z_min_global) bc_up = &supp_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_up == "stop" ) {
                if (z_min==z_min_global) bc_up = &stop_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_up == "periodic" ) {
                if (z_min==z_min_global) bc_up = &periodic_particle;
            }

        }//nDim_particle>2

    }//nDim_particle>1

    /* NOT USED ANYMORE AS WE USE THE ERFINV FCT FROM TOOLS/USERFUNCTIONS
    // ---------------------------------------------------------------------
    // Compute the tabulated inverse error function used in thermalizing bcs
    // ---------------------------------------------------------------------
    if ( thermCond ) {
        erfinv::instance().prepare();
    }//thermCond
     */

}


PartBoundCond::~PartBoundCond()
{
}

