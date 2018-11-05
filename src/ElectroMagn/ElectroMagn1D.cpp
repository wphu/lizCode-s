#include "ElectroMagn1D.h"

#include <cmath>

#include <sstream>
#include <string>
#include <iostream>

#include "PicParams.h"
#include "Field1D.h"
#include "MF_Solver1D_Yee.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::ElectroMagn1D(PicParams &params, InputData &input_data):
ElectroMagn(params, input_data)
{
    oversize_ = oversize[0];

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step
    dx       = cell_length[0];
    dt_ov_dx = timestep/cell_length[0];
    dx_ov_dt = 1.0/dt_ov_dx;

    // Electromagnetic fields
    // ----------------------
    // number of nodes of the primal-grid
    nx_p = n_space[0]+1 + 2*oversize[0];
    // number of nodes of the dual-grid
    nx_d = n_space[0]+2 + 2*oversize[0];
    // dimPrim/dimDual = nx_p/nx_d
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+1;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }

    // Allocation of the EM fields
    phi_ = new Field1D(dimPrim, "Phi" );
    Ex_  = new Field1D(dimPrim, "Ex");
    Ey_  = new Field1D(dimPrim, "Ey");
    Ez_  = new Field1D(dimPrim, "Ez");
    Bx_  = new Field1D(dimPrim, "Bx");
    By_  = new Field1D(dimPrim, "By");
    Bz_  = new Field1D(dimPrim, "Bz");
    Bx_m = new Field1D(dimPrim, "Bx_m");
    By_m = new Field1D(dimPrim, "By_m");
    Bz_m = new Field1D(dimPrim, "Bz_m");

    // for (unsigned int i=0 ; i<nx_d ; i++) {
    //         double x = ( (double)(smpi1D->getCellStartingGlobalIndex(0)+i-0.5) )*params.cell_length[0];
    //         (*By_)(i) = 0.001 * sin(x * 2.0*M_PI/params.sim_length[0] * 40.0);
    //     }
    //     smpi1D->exchangeField(By_);
    //     for (unsigned int i=0 ; i<nx_d ; i++) {
    // //        double x = ( (double)(smpi1D->getCellStartingGlobalIndex(0)+i-0.5) )*params.cell_length[0];
    //         (*By_m)(i) = (*By_)(i);
    //     }
    //
    // Allocation of time-averaged EM fields
    phi_avg = new Field1D(dimPrim, "Phi_avg" );
    Ex_avg  = new Field1D(dimPrim, "Ex_avg");
    Ey_avg  = new Field1D(dimPrim, "Ey_avg");
    Ez_avg  = new Field1D(dimPrim, "Ez_avg");
    Bx_avg  = new Field1D(dimPrim, "Bx_avg");
    By_avg  = new Field1D(dimPrim, "By_avg");
    Bz_avg  = new Field1D(dimPrim, "Bz_avg");

    // Total charge currents and densities
    Jx_   = new Field1D(dimPrim, "Jx");
    Jy_   = new Field1D(dimPrim, "Jy");
    Jz_   = new Field1D(dimPrim, "Jz");
    rho_  = new Field1D(dimPrim, "Rho" );
    rho_avg  = new Field1D(dimPrim, "Rho_avg" );

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

    MESSAGE("B "<<params.externB[0]<<" "<<params.externB[1]<<" "<<params.externB[2]);

    // Charge currents currents and density for each species

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]             = new Field1D(dimPrim, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[ispec]             = new Field1D(dimPrim, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[ispec]             = new Field1D(dimPrim, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[ispec]            = new Field1D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
        rho_s_avg[ispec]        = new Field1D(dimPrim, ("Rho_"+params.species_param[ispec].species_type+"_avg").c_str());


        Vx_s[ispec]             = new Field1D(dimPrim, ("Vx_"+params.species_param[ispec].species_type).c_str());
        Vx_s_avg[ispec]         = new Field1D(dimPrim, ("Vx_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vy_s[ispec]             = new Field1D(dimPrim, ("Vy_"+params.species_param[ispec].species_type).c_str());
        Vy_s_avg[ispec]         = new Field1D(dimPrim, ("Vy_"+params.species_param[ispec].species_type+"_avg").c_str());

        Vz_s[ispec]             = new Field1D(dimPrim, ("Vz_"+params.species_param[ispec].species_type).c_str());
        Vz_s_avg[ispec]         = new Field1D(dimPrim, ("Vz_"+params.species_param[ispec].species_type+"_avg").c_str());


        Vp_s[ispec]             = new Field1D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type).c_str());
        Vp_s_avg[ispec]         = new Field1D(dimPrim, ("Vparallel_"+params.species_param[ispec].species_type+"_avg").c_str());

        T_s[ispec]              = new Field1D(dimPrim, ("T_"+params.species_param[ispec].species_type).c_str());
        T_s_avg[ispec]          = new Field1D(dimPrim, ("T_"+params.species_param[ispec].species_type+"_avg").c_str());
    }


    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (size_t i=0 ; i<nDim_field ; i++) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }


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

}//END constructor Electromagn1D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::~ElectroMagn1D()
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::saveMagneticFields()
{
    // Static cast of the fields
    Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
    Field1D* By1D   = static_cast<Field1D*>(By_);
    Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
    Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
    Field1D* By1D_m = static_cast<Field1D*>(By_m);
    Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);

    // for Bx^(p)
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*Bx1D_m)(i)=(*Bx1D)(i);
    }
    //for By^(d) & Bz^(d)
    for (unsigned int i=0 ; i<dimDual[0] ; i++) {
        (*By1D_m)(i) = (*By1D)(i);
        (*Bz1D_m)(i) = (*Bz1D)(i);
    }

}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwellAmpere()
{

    Field1D* Ex1D = static_cast<Field1D*>(Ex_);
    Field1D* Ey1D = static_cast<Field1D*>(Ey_);
    Field1D* Ez1D = static_cast<Field1D*>(Ez_);
    Field1D* By1D = static_cast<Field1D*>(By_);
    Field1D* Bz1D = static_cast<Field1D*>(Bz_);
    Field1D* Jx1D = static_cast<Field1D*>(Jx_);
    Field1D* Jy1D = static_cast<Field1D*>(Jy_);
    Field1D* Jz1D = static_cast<Field1D*>(Jz_);

    // --------------------
    // Solve Maxwell-Ampere
    // --------------------
    // Calculate the electrostatic field ex on the dual grid
    //for (unsigned int ix=0 ; ix<nx_d ; ix++){
    for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
        (*Ex1D)(ix)= (*Ex1D)(ix) - timestep* (*Jx1D)(ix) ;
    }
    // Transverse fields ey, ez  are defined on the primal grid
    //for (unsigned int ix=0 ; ix<nx_p ; ix++) {
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        (*Ey1D)(ix)= (*Ey1D)(ix) - dt_ov_dx * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - timestep * (*Jy1D)(ix) ;
        (*Ez1D)(ix)= (*Ez1D)(ix) + dt_ov_dx * ( (*By1D)(ix+1) - (*By1D)(ix)) - timestep * (*Jz1D)(ix) ;
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::centerMagneticFields()
{
    // Static cast of the fields
    Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
    Field1D* By1D   = static_cast<Field1D*>(By_);
    Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
    Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
    Field1D* By1D_m = static_cast<Field1D*>(By_m);
    Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);

    // for Bx^(p)
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*Bx1D_m)(i) = ( (*Bx1D)(i)+ (*Bx1D_m)(i))*0.5 ;
    }

    // for By^(d) & Bz^(d)
    for (unsigned int i=0 ; i<dimDual[0] ; i++) {
        (*By1D_m)(i)= ((*By1D)(i)+(*By1D_m)(i))*0.5 ;
        (*Bz1D_m)(i)= ((*Bz1D)(i)+(*Bz1D_m)(i))*0.5 ;
    }

}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::incrementAvgFields(unsigned int time_step)
{
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1) % dump_step == 0 ){
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            rho_s_avg[ispec]->put_to(0.0);
        }//END loop on species ispec
    }

    // Calculate the sum values for global rho phi Ex and Ey
    if( (time_step % dump_step) > (dump_step - avg_step) || (time_step % dump_step) == 0 )
    {
        // Calculate the sum values for density of each species
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            // all fields are defined on the primal grid
            for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
                (*rho_s_avg[ispec])(ix) += (*rho_s[ispec])(ix);
            }
        }//END loop on species ispec
    }

    // calculate the averaged values
    if ( time_step % dump_step == 0 )
    {
        for (unsigned int ispec=0; ispec<n_species; ispec++) 
        {
            for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) 
            {
                (*rho_s_avg[ispec])(ix) /= avg_step;
            }
        }//END loop on species ispec
    }

}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::restartRhoJ()
{
    //Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
    //Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
    //Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
    Field1D* rho1D   = static_cast<Field1D*>(rho_);

    // --------------------------
    // Total currents and density
    // --------------------------

    // put longitudinal current to zero on the dual grid
    ///for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
    //    (*Jx1D)(ix)    = 0.0;
    //}

    // all fields are defined on the primal grid
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        (*rho1D)(ix)   = 0.0;
        //(*Jy1D)(ix)    = 0.0;
        //(*Jz1D)(ix)    = 0.0;
    }
}
void ElectroMagn1D::restartRhoJs(int ispec, bool currents)
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
    Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
    Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
    Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);

    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++)
    {
        (*rho1D_s)(ix) = 0.0;
    }
    if (currents){
        // put longitudinal current to zero on the dual grid
        for (unsigned int ix=0 ; ix<dimDual[0] ; ix++)
        {
            (*Jx1D_s)(ix)  = 0.0;
        }
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++)
        {
            // all fields are defined on the primal grid
            (*Jy1D_s)(ix)  = 0.0;
            (*Jz1D_s)(ix)  = 0.0;
        }
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::computeTotalRhoJ()
{
    //Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
    //Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
    //Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
    Field1D* rho1D   = static_cast<Field1D*>(rho_);

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        //Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
        //Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
        //Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
        Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);

        // put longitudinal current to zero on the dual grid
        //for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
        //    (*Jx1D)(ix)  += (*Jx1D_s)(ix);
        //}

        // all fields are defined on the primal grid
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
            //(*Jy1D)(ix)  += (*Jy1D_s)(ix);
            //(*Jz1D)(ix)  += (*Jz1D_s)(ix);
            (*rho1D)(ix) += ( (*rho1D_s)(ix) * species_param[ispec].charge );
        }
    }//END loop on species ispec
}

// --------------------------------------------------------------------------
// Compute Poynting (return the electromagnetic energy injected at the border
// --------------------------------------------------------------------------
void ElectroMagn1D::computePoynting() {

}
