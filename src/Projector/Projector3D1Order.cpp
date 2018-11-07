#include "Projector3D1Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D1Order
// ---------------------------------------------------------------------------------------------------------------------
Projector3D1Order::Projector3D1Order (PicParams& params) : Projector3D(params)
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    dz_inv_   = 1.0/params.cell_length[2];
    dz_ov_dt  = params.cell_length[2] / params.timestep;

    one_third = 1.0/3.0;

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D1Order
// ---------------------------------------------------------------------------------------------------------------------
Projector3D1Order::~Projector3D1Order()
{
}


//! Below, in this order :
//!   Project global current densities (EMfields->Jx_/Jy_/Jz_), not used
//!   Projection by species
//!   Project global current charge
//!   Project local current densities (sort)
//!   Project global current densities (ionize)




// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (EMfields->Jx_/Jy_/Jz_), not used
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D1Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double gf)
{

} // END Project global current densities, not used


// ---------------------------------------------------------------------------------------------------------------------
//!   Projection by species
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D1Order::operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf)
{

}//END Projection by species



// ---------------------------------------------------------------------------------------------------------------------
//! Project global current charge
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D1Order::operator() (Field* rho, Particles &particles, int ipart, double weight)
{

    //Static cast of the total charge density
    Field3D* rho3D  = static_cast<Field3D*>(rho);

    //Declaration of local variables
    double delta, delta2;
    double rho_p = weight;   // charge density of the macro-particle
    double Sx[2], Sy[2], Sz[2];             // projection coefficient arrays

    //Locate particle on the primal grid & calculate the projection coefficients
    double       xpn = particles.position(0, ipart) * dx_inv_;  // normalized distance to the first node
    int ic  = floor(xpn);                   // index of the central node
    delta  = xpn - (double)ic;                       // normalized distance to the nearest grid point
    Sx[0]  = 1.0-delta;
    Sx[1]  = delta;

    double       ypn = particles.position(1, ipart) * dy_inv_;  // normalized distance to the first node
    int jc   = floor(ypn);                  // index of the central node
    delta  = ypn - (double)jc;                       // normalized distance to the nearest grid point
    Sy[0]  = 1.0-delta;
    Sy[1]  = delta;

    double       zpn = particles.position(2, ipart) * dz_inv_;  // normalized distance to the first node
    int kc   = floor(zpn);                  // index of the central node
    delta  = zpn - (double)kc;                       // normalized distance to the nearest grid point
    Sz[0]  = 1.0-delta;
    Sz[1]  = delta;

    //cout << "Pos = " << particles.position(0, ipart) << " - i global = " << i << " - i local = " << i-index_domain_begin <<endl;
    int i = ic;
    int j = jc;
    int k = kc;

    // 1nd order projection for the total charge density
    for (unsigned int iloc=0 ; iloc<2 ; iloc++) 
    {
        for (unsigned int jloc=0 ; jloc<2 ; jloc++) 
        {
            for(unsigned int kloc=0 ; kloc<2 ; kloc++)
            {
                DEBUGEXEC
                (
                    if(i >= rho3D->dims_[0] || j >= rho3D->dims_[1] || k >= rho3D->dims_[2])
                    {
                        ERROR("Project Error, positions are: "<<particles.position(0, ipart)<<" "<<particles.position(1, ipart)<<" "<<particles.position(2, ipart));
                        //ERROR("Project Error, ic, jc, kc are: "<<ic<<" "<<jc<<" "<<kc);
                        //ERROR("Project Error, i, j, k are: "<<i<<" "<<j<<" "<<k);
                        //ERROR("Project Error, i_domain_begin are: "<<i_domain_begin<<" "<<j_domain_begin<<" "<<k_domain_begin);
                    }
                )
                (*rho3D)(i+iloc, j+jloc, k+kloc) += Sx[iloc]*Sy[jloc]*Sz[kloc]*rho_p;
            }
            
        }
    }

} // END Project global current charge



// ---------------------------------------------------------------------------------------------------------------------
//! Project global current charge
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D1Order::operator() (Field* rho, Particles &particles, int ipart)
{

} // END Project global current charge


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D1Order::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim1)
{

} // END Project local current densities (sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D1Order::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
    ERROR("Projection of ionization current not yet defined for 3D 2nd order");

} // END Project global current densities (ionize)
