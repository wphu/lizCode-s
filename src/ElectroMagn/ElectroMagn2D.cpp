#include "ElectroMagn2D.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "Field2D.h"


#include "Profile.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::ElectroMagn2D(PicParams &params, InputData &input_data) :
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

    // ----------------------
    // Electromagnetic fields
    // ----------------------
    //! \todo Homogenize 2D/2D dimPrim/dimDual or nx_p/nx_d/ny_p/ny_d

    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    dim_global.resize( nDim_field );

    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];

        dim_global[i] = n_space_global[i] + 1;
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = n_space[0]+2+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = n_space[1]+2+2*oversize[1];

    // Allocation of the EM fields
    Ex_  = new Field2D(dimPrim, "Ex" );
    Ey_  = new Field2D(dimPrim, "Ey" );
    Ez_  = new Field2D(dimPrim, "Ez" );
    Bx_  = new Field2D(dimPrim, "Bx" );
    By_  = new Field2D(dimPrim, "By" );
    Bz_  = new Field2D(dimPrim, "Bz" );
    Bx_m = new Field2D(dimPrim, "Bx_m" );
    By_m = new Field2D(dimPrim, "By_m" );
    Bz_m = new Field2D(dimPrim, "Bz_m" );

    // Total charge currents and densities
    Jx_   = new Field2D(dimPrim, "Jx" );
    Jy_   = new Field2D(dimPrim, "Jy" );
    Jz_   = new Field2D(dimPrim, "Jz" );
    rho_  = new Field2D(dimPrim, "Rho" );
    rho_avg  = new Field2D(dimPrim, "Rho_avg" );

    // Allocation of time-averaged EM fields
    phi_avg = new Field2D(dimPrim, "Phi_avg");
    Ex_avg  = new Field2D(dimPrim, "Ex_avg");
    Ey_avg  = new Field2D(dimPrim, "Ey_avg");
    Ez_avg  = new Field2D(dimPrim, "Ez_avg");
    Bx_avg  = new Field2D(dimPrim, "Bx_avg");
    By_avg  = new Field2D(dimPrim, "By_avg");
    Bz_avg  = new Field2D(dimPrim, "Bz_avg");

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
    Ex_avg  = new Field2D(dimPrim, "Ex_avg" );
    Ey_avg  = new Field2D(dimPrim, "Ey_avg" );
    Ez_avg  = new Field2D(dimPrim, "Ez_avg" );
    Bx_avg  = new Field2D(dimPrim, "Bx_avg" );
    By_avg  = new Field2D(dimPrim, "By_avg" );
    Bz_avg  = new Field2D(dimPrim, "Bz_avg" );

    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field2D(dimPrim, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[ispec]  = new Field2D(dimPrim, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[ispec]  = new Field2D(dimPrim, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[ispec] = new Field2D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
        rho_s_avg[ispec]        = new Field2D(dimPrim, ("Rho_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vx_s[ispec]             = new Field2D(dimPrim, ("Vx_"+params.species_param[ispec].species_type).c_str());
        Vx_s_avg[ispec]         = new Field2D(dimPrim, ("Vx_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vy_s[ispec]             = new Field2D(dimPrim, ("Vy_"+params.species_param[ispec].species_type).c_str());
        Vy_s_avg[ispec]         = new Field2D(dimPrim, ("Vy_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vz_s[ispec]             = new Field2D(dimPrim, ("Vz_"+params.species_param[ispec].species_type).c_str());
        Vz_s_avg[ispec]         = new Field2D(dimPrim, ("Vz_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vp_s[ispec]             = new Field2D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type).c_str());
        Vp_s_avg[ispec]         = new Field2D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type+"_avg").c_str());

        T_s[ispec]              = new Field2D(dimPrim, ("T_"+params.species_param[ispec].species_type).c_str());
        T_s_avg[ispec]          = new Field2D(dimPrim, ("T_"+params.species_param[ispec].species_type+"_avg").c_str());


    }

    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (unsigned int i=0 ; i<nDim_field ; i++) {
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
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            istart[i][isDual] = 0;
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++) {
            istart[i][isDual] = oversize[i];
        }
    }

    // bufsize = nelements
    for (unsigned int i=0 ; i<3 ; i++)
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            bufsize[i][isDual] = 1;

    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (int isDual=0 ; isDual<2 ; isDual++)
            bufsize[i][isDual] = n_space[i] + 1;

        for (int isDual=0 ; isDual<2 ; isDual++) {
            bufsize[i][isDual] += isDual;
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field

}//END constructor Electromagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::~ElectroMagn2D()
{
}//END ElectroMagn2D


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::saveMagneticFields()
{
    // Static cast of the fields
    Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
    Field2D* By2D   = static_cast<Field2D*>(By_);
    Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
    Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);

    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bx2D_m)(i,j)=(*Bx2D)(i,j);
        }

    // Magnetic field By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j)=(*By2D)(i,j);
        }

    // Magnetic field Bz^(d,d)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bz2D_m)(i,j)=(*Bz2D)(i,j);
        }
    }// end for j
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(nx_p,j)=(*By2D)(nx_p,j);
        }
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bz2D_m)(nx_p,j)=(*Bz2D)(nx_p,j);
        }

}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Ampere equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwellAmpere()
{
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(Bx_);
    Field2D* By2D = static_cast<Field2D*>(By_);
    Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    Field2D* Jx2D = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D = static_cast<Field2D*>(Jz_);

    // Electric field Ex^(d,p)

//{
    for (unsigned int i=0 ; i<nx_p ; i++) {
//    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex2D)(i,j) += -timestep*(*Jx2D)(i,j) + dt_ov_dy * ( (*Bz2D)(i,j+1) - (*Bz2D)(i,j) );
        }// end for j
//    }// end for i

    // Electric field Ey^(p,d)
//    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ey2D)(i,j) += -timestep*(*Jy2D)(i,j) - dt_ov_dx * ( (*Bz2D)(i+1,j) - (*Bz2D)(i,j) );
        }// end for j
    //} // end for i

    // Electric field Ez^(p,p)
    //for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D)(i,j) += -timestep*(*Jz2D)(i,j)
            +               dt_ov_dx * ( (*By2D)(i+1,j) - (*By2D)(i,j) )
            -               dt_ov_dy * ( (*Bx2D)(i,j+1) - (*Bx2D)(i,j) );
        } // end for j
    }// end for i
//} // end parallel

    for (unsigned int j=0 ; j<ny_p ; j++) {
        (*Ex2D)(nx_p,j) += -timestep*(*Jx2D)(nx_p,j) + dt_ov_dy * ( (*Bz2D)(nx_p,j+1) - (*Bz2D)(nx_p,j) );
    }

//    }
}//END solveMaxwellAmpere


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::centerMagneticFields()
{
    // Static cast of the fields
    Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
    Field2D* By2D   = static_cast<Field2D*>(By_);
    Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
    Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);

    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bx2D_m)(i,j) = ( (*Bx2D)(i,j) + (*Bx2D_m)(i,j) )*0.5;
        }
//    }

    // Magnetic field By^(d,p)
//    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j) = ( (*By2D)(i,j) + (*By2D_m)(i,j) )*0.5;
        }
//    }

    // Magnetic field Bz^(d,d)
//    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Bz2D_m)(i,j) = ( (*Bz2D)(i,j) + (*Bz2D_m)(i,j) )*0.5;
        } // end for j
      } // end for i

    for (unsigned int j=0 ; j<ny_p ; j++) {
        (*By2D_m)(nx_p,j) = ( (*By2D)(nx_p,j) + (*By2D_m)(nx_p,j) )*0.5;
    }
    for (unsigned int j=0 ; j<ny_p ; j++) {
        (*Bz2D_m)(nx_p,j) = ( (*Bz2D)(nx_p,j) + (*Bz2D_m)(nx_p,j) )*0.5;
    } // end for j


}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::incrementAvgFields(unsigned int time_step)
{
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1) % dump_step == 0 )
    {
        for (unsigned int ispec=0; ispec<n_species; ispec++) 
        {
            rho_s_avg[ispec]->put_to(0.0);
        }//END loop on species ispec
    }

    // Calculate the sum values for global rho phi Ex and Ey
    if( (time_step % dump_step) > (dump_step - avg_step) || (time_step % dump_step) == 0 )
    {
        // Calculate the sum values for density of each species
        for (unsigned int ispec=0; ispec<n_species; ispec++) 
        {
            // all fields are defined on the primal grid
            for (unsigned int ix=0 ; ix<dimPrim[0]*dimPrim[1] ; ix++) 
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
            for (unsigned int ix=0 ; ix<dimPrim[0]*dimPrim[1] ; ix++) 
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
void ElectroMagn2D::restartRhoJ()
{
    // --------------------------
    // Total currents and density
    // --------------------------

    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);

    // Charge density rho^(p,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*rho2D)(i,j) = 0.0;
        }
    }

    // Current Jx^(d,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jx2D)(i,j) = 0.0;
        }
    }

    // Current Jy^(p,d) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jy2D)(i,j) = 0.0;
        }
    }

    // Current Jz^(p,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jz2D)(i,j) = 0.0;
        }
    }
}//END restartRhoJ


void ElectroMagn2D::restartRhoJs(int ispec, bool currents)
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
    Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
    Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
    Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);

    // Charge density rho^(p,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++)
    {
        for (unsigned int j=0 ; j<ny_p ; j++)
        {
            (*rho2D_s)(i,j) = 0.0;
        }
    }
    if (currents){
        // Current Jx^(d,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++)
        {
            for (unsigned int j=0 ; j<ny_p ; j++)
            {
                (*Jx2D_s)(i,j) = 0.0;
            }
        }

        // Current Jy^(p,d) to 0
        for (unsigned int i=0 ; i<nx_p ; i++)
        {
            for (unsigned int j=0 ; j<ny_p ; j++)
            {
                (*Jy2D_s)(i,j) = 0.0;
            }
        }

        // Current Jz^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++)
        {
            for (unsigned int j=0 ; j<ny_p ; j++)
            {
                (*Jz2D_s)(i,j) = 0.0;
            }
        }
    }
}//END restartRhoJs




// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::computeTotalRhoJ()
{

    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);


    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
        Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
        Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
        Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);

        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*rho2D)(i,j) += ( (*rho2D_s)(i,j) * species_param[ispec].charge );
            }
        }
/*
        // Current Jx^(d,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jx2D)(i,j) += (*Jx2D_s)(i,j);
            }
        }

        // Current Jy^(p,d) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jy2D)(i,j) += (*Jy2D_s)(i,j);
            }
        }

        // Current Jz^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jz2D)(i,j) += (*Jz2D_s)(i,j);
            }
        }
*/
    }//END loop on species ispec

}//END computeTotalRhoJ


