#include "ElectroMagn3D.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "Field3D.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::ElectroMagn3D(PicParams &params, InputData &input_data) :
ElectroMagn(params, input_data)
{
    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = cell_length[0];
    dt_ov_dx = timestep/dx;
    dx_ov_dt = 1.0/dt_ov_dx;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = cell_length[1];
    dt_ov_dy = timestep/dy;
    dy_ov_dt = 1.0/dt_ov_dy;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dz       = cell_length[2];
    dt_ov_dz = timestep/dz;
    dz_ov_dt = 1.0/dt_ov_dz;

    // ----------------------
    // Electromagnetic fields
    // ----------------------
    //! \todo Homogenize 3D/3D dimPrim/dimDual or nx_p/nx_d/ny_p/ny_d

    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );

    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<nDim_field ; i++) 
    {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = n_space[0]+2+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = n_space[1]+2+2*oversize[1];
    // number of nodes of the primal and dual grid in the y-direction
    nz_p = n_space[2]+1+2*oversize[2];
    nz_d = n_space[2]+2+2*oversize[2];


    // Allocation of the EM fields
    Ex_  = new Field3D(dimPrim, "Ex" );
    Ey_  = new Field3D(dimPrim, "Ey" );
    Ez_  = new Field3D(dimPrim, "Ez" );
    Bx_  = new Field3D(dimPrim, "Bx" );
    By_  = new Field3D(dimPrim, "By" );
    Bz_  = new Field3D(dimPrim, "Bz" );
    Bx_m = new Field3D(dimPrim, "Bx_m" );
    By_m = new Field3D(dimPrim, "By_m" );
    Bz_m = new Field3D(dimPrim, "Bz_m" );

    // Total charge currents and densities
    Jx_   = new Field3D(dimPrim, "Jx" );
    Jy_   = new Field3D(dimPrim, "Jy" );
    Jz_   = new Field3D(dimPrim, "Jz" );
    rho_  = new Field3D(dimPrim, "Rho" );
    rho_avg  = new Field3D(dimPrim, "Rho_avg" );

    // Allocation of time-averaged EM fields
    phi_avg = new Field3D(dimPrim, "Phi_avg" );
    Ex_avg  = new Field3D(dimPrim, "Ex_avg");
    Ey_avg  = new Field3D(dimPrim, "Ey_avg");
    Ez_avg  = new Field3D(dimPrim, "Ez_avg");
    Bx_avg  = new Field3D(dimPrim, "Bx_avg");
    By_avg  = new Field3D(dimPrim, "By_avg");
    Bz_avg  = new Field3D(dimPrim, "Bz_avg");

    Ex_->put_to(0.0);
    Ey_->put_to(0.0);
    Ez_->put_to(0.0);
    Bx_->put_to(params.externB[0]);
    By_->put_to(params.externB[1]);
    Bz_->put_to(params.externB[2]);
    Bx_m->put_to(0.0);
    By_m->put_to(0.0);
    Bz_m->put_to(0.0);
    rho_->put_to(0.0);



    // Allocation of the time-averaged EM fields
    Ex_avg  = new Field3D(dimPrim, "Ex_avg" );
    Ey_avg  = new Field3D(dimPrim, "Ey_avg" );
    Ez_avg  = new Field3D(dimPrim, "Ez_avg" );
    Bx_avg  = new Field3D(dimPrim, "Bx_avg" );
    By_avg  = new Field3D(dimPrim, "By_avg" );
    Bz_avg  = new Field3D(dimPrim, "Bz_avg" );

    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) 
    {
        Jx_s[ispec]  = new Field3D(dimPrim, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[ispec]  = new Field3D(dimPrim, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[ispec]  = new Field3D(dimPrim, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[ispec] = new Field3D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
        rho_s_avg[ispec]        = new Field3D(dimPrim, ("Rho_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vx_s[ispec]             = new Field3D(dimPrim, ("Vx_"+params.species_param[ispec].species_type).c_str());
        Vx_s_avg[ispec]         = new Field3D(dimPrim, ("Vx_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vy_s[ispec]             = new Field3D(dimPrim, ("Vy_"+params.species_param[ispec].species_type).c_str());
        Vy_s_avg[ispec]         = new Field3D(dimPrim, ("Vy_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vz_s[ispec]             = new Field3D(dimPrim, ("Vz_"+params.species_param[ispec].species_type).c_str());
        Vz_s_avg[ispec]         = new Field3D(dimPrim, ("Vz_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vp_s[ispec]             = new Field3D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type).c_str());
        Vp_s_avg[ispec]         = new Field3D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type+"_avg").c_str());

        T_s[ispec]              = new Field3D(dimPrim, ("T_"+params.species_param[ispec].species_type).c_str());
        T_s_avg[ispec]          = new Field3D(dimPrim, ("T_"+params.species_param[ispec].species_type+"_avg").c_str());

    }



    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (unsigned int i=0 ; i<nDim_field ; i++) 
    {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    /*
     MESSAGE("index_bc_min / index_bc_max / nx_p / nx_d" << index_bc_min[0]
            << " " << index_bc_max[0] << " " << nx_p<< " " << nx_d);
     */


    // Define limits of non duplicated elements
    // (by construction 1 (prim) or 2 (dual) elements shared between per MPI process)
    // istart
    for (unsigned int i=0 ; i<3 ; i++)
    {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
        {
            istart[i][isDual] = 0;
        }
    }
        
    for (unsigned int i=0 ; i<nDim_field ; i++) 
    {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++) 
        {
            istart[i][isDual] = oversize[i];
        }
    }

    // bufsize = nelements
    for (unsigned int i=0 ; i<3 ; i++)
    {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
        {
            bufsize[i][isDual] = 1;
        }
    }
        

    for (unsigned int i=0 ; i<nDim_field ; i++) 
    {
        for (int isDual=0 ; isDual<2 ; isDual++)
        {
            bufsize[i][isDual] = n_space[i] + 1;
        }

        for (int isDual=0 ; isDual<2 ; isDual++) 
        {
            bufsize[i][isDual] += isDual;
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field

}//END constructor Electromagn3D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::~ElectroMagn3D()
{
}//END ElectroMagn3D


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::saveMagneticFields()
{
    // Static cast of the fields
    Field3D* Bx3D   = static_cast<Field3D*>(Bx_);
    Field3D* By3D   = static_cast<Field3D*>(By_);
    Field3D* Bz3D   = static_cast<Field3D*>(Bz_);
    Field3D* Bx3D_m = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m = static_cast<Field3D*>(Bz_m);


}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Ampere equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::solveMaxwellAmpere()
{
    // Static-cast of the fields
    Field3D* Ex3D = static_cast<Field3D*>(Ex_);
    Field3D* Ey3D = static_cast<Field3D*>(Ey_);
    Field3D* Ez3D = static_cast<Field3D*>(Ez_);
    Field3D* Bx3D = static_cast<Field3D*>(Bx_);
    Field3D* By3D = static_cast<Field3D*>(By_);
    Field3D* Bz3D = static_cast<Field3D*>(Bz_);
    Field3D* Jx3D = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D = static_cast<Field3D*>(Jz_);



}//END solveMaxwellAmpere


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::centerMagneticFields()
{
    // Static cast of the fields
    Field3D* Bx3D   = static_cast<Field3D*>(Bx_);
    Field3D* By3D   = static_cast<Field3D*>(By_);
    Field3D* Bz3D   = static_cast<Field3D*>(Bz_);
    Field3D* Bx3D_m = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m = static_cast<Field3D*>(Bz_m);


}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::incrementAvgFields(unsigned int time_step)
{
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1) % dump_step == 0 )
    {
        for (unsigned int ispec=0; ispec<n_species; ispec++)
        {
            rho_s_avg[ispec]->put_to(0.0);
        }
    }

    // Calculate the sum values for global rho phi Ex and Ey
    if( (time_step % dump_step) > (dump_step - avg_step) || (time_step % dump_step) == 0 )
    {
        // Calculate the sum values for density of each species
        for (unsigned int ispec=0; ispec<n_species; ispec++) 
        {
            // all fields are defined on the primal grid
            for (unsigned int ix=0 ; ix<dimPrim[0]*dimPrim[1]*dimPrim[2] ; ix++) 
            {
                (*rho_s_avg[ispec])(ix) += (*rho_s[ispec])(ix);
            }
        }//END loop on species ispec
    }

    // calculate the averaged values
    if ( time_step % dump_step == 0 )
    {
        for (unsigned int ispec=0; ispec<n_species; ispec++) 
        {
            for (unsigned int ix=0 ; ix<dimPrim[0]*dimPrim[1]*dimPrim[2] ; ix++) 
            {
                (*rho_s_avg[ispec])(ix) /= avg_step;
            }
        }//END loop on species ispec
    }

}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge densities and currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::restartRhoJ()
{
    // --------------------------
    // Total currents and density
    // --------------------------

    // static cast of the total currents and densities
    Field3D* Jx3D    = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D    = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D    = static_cast<Field3D*>(Jz_);
    Field3D* rho3D   = static_cast<Field3D*>(rho_);

    // Charge density rho^(p,p) to 0

    rho3D->put_to(0.0);
    //Jx3D->put_to(0.0);
    //Jy3D->put_to(0.0);
    //Jz3D->put_to(0.0);

}//END restartRhoJ


void ElectroMagn3D::restartRhoJs(int ispec, bool currents)
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    Field3D* Jx3D_s  = static_cast<Field3D*>(Jx_s[ispec]);
    Field3D* Jy3D_s  = static_cast<Field3D*>(Jy_s[ispec]);
    Field3D* Jz3D_s  = static_cast<Field3D*>(Jz_s[ispec]);
    Field3D* rho3D_s = static_cast<Field3D*>(rho_s[ispec]);

    rho3D_s->put_to(0.0);

    if (currents)
    {
        Jx3D_s->put_to(0.0);
        Jy3D_s->put_to(0.0);
        Jz3D_s->put_to(0.0);
    }
}//END restartRhoJs




// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computeTotalRhoJ()
{

    // static cast of the total currents and densities
    Field3D* Jx3D    = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D    = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D    = static_cast<Field3D*>(Jz_);
    Field3D* rho3D   = static_cast<Field3D*>(rho_);


    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) 
    {
        Field3D* Jx3D_s  = static_cast<Field3D*>(Jx_s[ispec]);
        Field3D* Jy3D_s  = static_cast<Field3D*>(Jy_s[ispec]);
        Field3D* Jz3D_s  = static_cast<Field3D*>(Jz_s[ispec]);
        Field3D* rho3D_s = static_cast<Field3D*>(rho_s[ispec]);

        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) 
        {
            for (unsigned int j=0 ; j<ny_p ; j++) 
            {
                for (unsigned int k=0 ; k<nz_p ; k++) 
                {
                    (*rho3D)(i,j,k) += ( (*rho3D_s)(i,j,k) * species_param[ispec].charge );
                }
            }
        }

    }

}//END computeTotalRhoJ

