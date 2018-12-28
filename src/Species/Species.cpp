#include "Species.h"

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>

#include <iostream>

#include <omp.h>

#include "PusherFactory.h"

#include "PartBoundCond.h"
//#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "Interpolator1D1Order_test.h"
#include "Grid.h"
#include "Grid2D.h"
#include "Projector.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

#include "Timer.h"

using namespace std;

Species::Species(PicParams &params, int ispec) : speciesNumber(ispec),
                                                 oversize(params.oversize),
                                                 cell_length(params.cell_length),
                                                 species_param(params.species_param[ispec]),
                                                 ndim(params.nDim_particle)
{
    electron_species = NULL;

    n_cell = params.n_space[0];
    
    q = species_param.charge;
    m = species_param.mass;
    q_over_m = q / m;
    weight = species_param.density / species_param.n_part_per_cell;


    particles.resize(n_cell);
    if (!params.restart)
    {
        createParticles(params);
    }

    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond(params, ispec);

} //END Species creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete Push;
    if (partBoundCond)
        delete partBoundCond;
    DEBUG(10, "Species deleted ");
}

int Species::createParticles(PicParams &params)
{
    for (int i_cell = 0; i_cell < n_cell; i++)
    {
        particles[i_cell].reserve(species_param.n_part_per_cell * species_param.c_part_max);
        particles[i_cell].resize(species_param.n_part_per_cell);
        n_particle_in_cell[i_cell] = species_param.n_part_per_cell;
        initPosition(i_cell, species_param.initPosition_type);
        initMomentum(i_cell, species_param.initMomentum_type);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// for all particles in a mesh initialize their position
// local cordinate and normalization by cell_length is used, so x = 0 ~ 1
// ---------------------------------------------------------------------------------------------------------------------
void Species::initPosition(int i_cell, string initPosition_type)
{
    for (int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
    {
        if (initPosition_type == "regular")
        {
            particles[i_cell][i_particle].x = i_particle / n_part_per_cell[i_cell];
        }
        else if (initPosition_type == "random")
        {
            particles[i_cell][i_particle].x = (double)rand() / RAND_MAX;
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles in a mesh initialize their velocity
// ---------------------------------------------------------------------------------------------------------------------
void Species::initVelocity(int i_cell, string initMomentum_type, PicParams &params)
{
    if (initMomentum_type == "cold")
    {
        for (int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        {
            particles[i_cell][i_particle].vx = 0.0;
            particles[i_cell][i_particle].vy = 0.0;
            particles[i_cell][i_particle].vz = 0.0;
        }
    }
    else if (initMomentum_type == "maxwell")
    {
        // initialize using the Maxwell distribution function
        for (int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        {
            double vt = sqrt(2.0 * species_param.temperature * params.const_e / species_param.mass);
            double x1;
            double x2;

            do
            {
                x1 = (double)rand() / RAND_MAX;
            } while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles[i_cell][i_particle].vx = vt * sqrt(-log(x1)) * sin(2.0 * M_PI * x2);

            do
            {
                x1 = (double)rand() / RAND_MAX;
            } while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles[i_cell][i_particle].vy = vt * sqrt(-log(x1)) * sin(2.0 * M_PI * x2);

            do
            {
                x1 = (double)rand() / RAND_MAX;
            } while (x1 == 0.0);
            x2 = (double)rand() / RAND_MAX;
            particles[i_cell][i_particle].vz = vt * sqrt(-log(x1)) * sin(2.0 * M_PI * x2);
        }
    }
    else if (initMomentum_type == "rectangular")
    {
        for (int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        {
            particles[i_cell][i_particle].vx = 2.0 * (2. * (double)rand() / RAND_MAX - 1.) * sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            particles[i_cell][i_particle].vy = (0.00001 * 2. * (double)rand() / RAND_MAX - 1.) * sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
            particles[i_cell][i_particle].vz = (0.00001 * 2. * (double)rand() / RAND_MAX - 1.) * sqrt(2.0 * temp[0] * params.const_e / species_param.mass);
        }
    }


    for (int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
    {
        particles[i_cell][i_particle].vx += species_param.mean_velocity[0];
        particles[i_cell][i_particle].vy += species_param.mean_velocity[1];
        particles[i_cell][i_particle].vz += species_param.mean_velocity[2];
    }


}

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void Species::heat(PicParams &params)
{


}


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics(double time_dual, unsigned int ispec, ElectroMagn *EMfields, Interpolator *Interp,
                       Projector *Proj, PicParams &params)
{
    LocalFields E_local, B_local;

    Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* rho1D  = static_cast<Field1D*>(EMfields->rho_s[species_number]);

    int iDirection = -1;

    clearExchList();


    psi_particles.clear();


    for(int i_cell = 0; i_cell < n_space; i_cell++)
    {
        for(int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        { 
            interpolator1D1Order(Ex1D,i_cell, i_particle, E_local);
            
            E_local.x = 0.0;
            //pusherBoris(i_cell, i_particle, E_local);
            pusherBoris(i_cell, i_particle, E_local, B_local);
        }
    }

    sort_particles();

    for(int i_cell = 0; i_cell < n_space; i_cell++)
    {
        for(int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        { 
            projector1D1Order(rho1D, i_cell, i_particle);
        }
    }

    // Apply boundary condition on the particles
    // Boundary Condition may be physical or due to domain decomposition
    // apply returns 0 if iPart is no more in the domain local
    // if omp, create a list per thread
    if (!partBoundCond->apply(particles, iPart, params.species_param[ispec], ener_iPart, iDirection))
    {
        addPartInExchList(iPart);
        if (iDirection >= 0)
        {
            addPartInPsiList(iDirection, iPart);
        }
    }


    // copy PSI particles to psi_particles, because after MPi particle exchanging
    // the PSI particles will be erased
    for (int iD = 0; iD < indexes_of_particles_to_perform_psi.size(); iD++)
    {
        for (int iPart = 0; iPart < indexes_of_particles_to_perform_psi[iD].size(); iPart++)
        {
            int iPart_psi = indexes_of_particles_to_perform_psi[iD][iPart];
            particles.cp_particle(iPart_psi, psi_particles);
        }
    }

    erase_particles_from_bins(indexes_of_particles_to_exchange);

}

void Species::interpolator1D1Order(Field1D* Ex1D, int i_cell, int i_particle, LocalFields &E_local)
{
    double xjn;

    xjn = particles[i_cell][i_particle].x;
    E_local.x = (*Ex1D)(i_cell) * (1.0 - xjn) + (*Ex1D)(i_cell) * xjn;
}

void Species::projector1D1Order(Field1D* rho1D, int i_cell, int i_particle)
{
    double xjn;

    xjn = particles[i_cell][i_particle].x;
    (*rho1D)(i_cell) += (1.0 - xjn);
    (*rho1D)(i_cell + 1) += xjn;
} 


void Species::pusherBoris(int i_cell, int i_particle, LocalFields &E_local)
{
    particles[i_cell][i_particle].vx += charge_over_mass * E_local.x * dt;
    particles[i_cell][i_particle].x  += particles[i_cell][i_particle].vx * dt_over_dx;
}

void Species::pusherBoris(int i_cell, int i_particle, LocalFields &E_local, LocalFields &B_local)
{
    // for pushBoris with magnetic field
    double umx, umy, umz, upx, upy, upz, pxdot, pydot, pzdot;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2, Sx, Sy, Sz;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    double dl;

    // --------------------------------------
    // SOLVE THE PARTICLE EQUATION OF MOTIONS
    // --------------------------------------

    // Half-acceleration in the electric field
    umx = particles[i_cell][i_particle].vx + charge_over_mass * E_local.x * dts2;
    umy = particles[i_cell][i_particle].vy + charge_over_mass * E_local.y * dts2;
    umz = particles[i_cell][i_particle].vz + charge_over_mass * E_local.z * dts2;


    // Rotation in the magnetic field
    alpha = charge_over_mass * dts2;
    Tx    = alpha * B_local.x;
    Ty    = alpha * B_local.y;
    Tz    = alpha * B_local.z;
    Tx2   = Tx*Tx;
    Ty2   = Ty*Ty;
    Tz2   = Tz*Tz;
    inv_det_T = 2.0/(1.0+Tx2+Ty2+Tz2);
    Sx    = Tx * inv_det_T;
    Sy    = Ty * inv_det_T;
    Sz    = Tz * inv_det_T;


    pxdot = umx + umy * Tz - umz * Ty;
    pydot = umy + umz * Tx - umx * Tz;
    pzdot = umz + umx * Ty - umy * Tx;

    upx = umx + pydot * Sz - pzdot * Sy;
    upy = umy + pzdot * Sx - pxdot * Sz;
    upz = umz + pxdot * Sy - pydot * Sx;


    // Half-acceleration in the electric field
    pxsm = upx + charge_over_mass * E_local.x * dts2;
    pysm = upy + charge_over_mass * E_local.y * dts2;
    pzsm = upz + charge_over_mass * E_local.z * dts2;

    particles[i_cell][i_particle].vx = pxsm;
    particles[i_cell][i_particle].vy = pysm;
    particles[i_cell][i_particle].vz = pzsm;

    // Move the particle
    particles[i_cell][i_particle].x += dt_over_dx * particles[i_cell][i_particle].vx;
}




void Species::sort_particles()
{
    int i_cell_change;
    int n_space_check_boundary;

    n_space_check_boundary = 5;


    //cout<<"aaa"<<endl;
    for(int i_cell = 0; i_cell < n_space_check_boundary; i_cell++)
    {
        for(int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        {
            i_cell_change =  floor(particles[i_cell][i_particle].x);

            if(i_cell + i_cell_change < 0)
            {
                particles[i_cell][i_particle].x = -particles[i_cell][i_particle].x - 2.0 * i_cell;
            }
        }
    }
    //cout<<"bbb"<<endl;
    for(int i_cell = n_space - n_space_check_boundary; i_cell < n_space; i_cell++)
    {
        for(int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        {
            i_cell_change =  floor(particles[i_cell][i_particle].x);

            if(i_cell + i_cell_change >= particles.size())
            {
                particles[i_cell][i_particle].x = -particles[i_cell][i_particle].x - 2.0 * (i_cell - n_space);
            }
        }
    }

    //cout<<"ccc"<<endl;
    for(int i_cell = 0; i_cell < n_space; i_cell++)
    {
        for(int i_particle = 0; i_particle < n_particle_in_cell[i_cell]; i_particle++)
        {
            i_cell_change =  floor(particles[i_cell][i_particle].x);

            if(i_cell_change != 0)
            {
                particles[i_cell][i_particle].x = abs(particles[i_cell][i_particle].x - i_cell_change);
                particles[i_cell + i_cell_change].push_back(particles[i_cell][i_particle]);
                list_removed_particles[i_cell].push_back(i_particle);
            }
        }
    }
    //cout<<"ddd"<<endl;
    for(int i_cell = 0; i_cell < n_space; i_cell++)
    {
        int i_particle_last = particles[i_cell].size() - 1;
        for(int i_list = 0; i_list < list_removed_particles[i_cell].size(); i_list++)
        {
            int i_particle = list_removed_particles[i_cell][list_removed_particles[i_cell].size() - 1 - i_list];
            if(i_particle == i_particle_last)
            {
                i_particle_last--;
            }
            else
            {
                swap(particles[i_cell][i_particle], particles[i_cell][i_particle_last]);
                i_particle_last--;
            }
        }
        particles[i_cell].resize(particles[i_cell].size() - list_removed_particles[i_cell].size());
        n_particle_in_cell[i_cell] =  particles[i_cell].size();
        list_removed_particles[i_cell].clear();
    } 
    //cout<<"eee"<<endl;
}


void Species::insert_particles_to_bins(Particles &insert_Particles, std::vector<int> &count_in_bins)
{
    int n_part_insert;
    int begin_id = 0;
    for (int ibin = 0; ibin < bmax.size(); ibin++)
    {
        n_part_insert = count_in_bins[ibin];
        insert_Particles.cp_particles(begin_id, n_part_insert, particles, bmax[ibin]);
        bmax[ibin] += n_part_insert;
        for (int ibin_inLoop = ibin + 1; ibin_inLoop < bmax.size(); ibin_inLoop++)
        {
            bmax[ibin_inLoop] += n_part_insert;
            bmin[ibin_inLoop] = bmax[ibin_inLoop - 1];
        }
        begin_id += n_part_insert;
    }
}

void Species::insert_particles(Particles &insert_particles)
{
    int n_part_insert;
    int begin_id = 0;
    int iDirection = -1;
    double ener_iPart = 0.0;
    double limit_min, limit_max;
    Particles insert_particles_temp;
    vector<int> count_in_bins;
    count_in_bins.resize(bmax.size());

    clearExchList();
    for (int ibin = 0; ibin < bmax.size(); ibin++)
    {
        limit_min = min_loc + ibin * cell_length[0];
        limit_max = min_loc + (ibin + 1) * cell_length[0];
        for (int iPart = 0; iPart < insert_particles.size(); iPart++)
        {
            if (insert_particles.position(0, iPart) >= limit_min && insert_particles.position(0, iPart) < limit_max)
            {
                count_in_bins[ibin]++;
                insert_particles.cp_particle(iPart, insert_particles_temp);
            }
        }
    }

    insert_particles_to_bins(insert_particles_temp, count_in_bins);

    // determine if particle is in other MPI region
    for (int ibin = 0; ibin < bmax.size(); ibin++)
    {
        for (int iPart = bmax[ibin] - count_in_bins[ibin]; iPart < bmax[ibin]; iPart++)
        {
            if (!partBoundCond->apply(particles, iPart, species_param, ener_iPart, iDirection))
            {
                addPartInExchList(iPart);
            }
        }
    }
}

void Species::erase_particles_from_bins(std::vector<int> &indexs_to_erase)
{
    int ii, iPart;
    for (unsigned int ibin = 0; ibin < bmax.size(); ibin++)
    {
        //        DEBUG(ibin << " bounds " << bmin[ibin] << " " << bmax[ibin]);
        ii = indexs_to_erase.size() - 1;
        if (ii >= 0)
        { // Push lost particles to the end of the bin
            iPart = indexs_to_erase[ii];
            while (iPart >= bmax[ibin] && ii > 0)
            {
                ii--;
                iPart = indexs_to_erase[ii];
            }
            while (iPart == bmax[ibin] - 1 && iPart >= bmin[ibin] && ii > 0)
            {
                bmax[ibin]--;
                ii--;
                iPart = indexs_to_erase[ii];
            }
            while (iPart >= bmin[ibin] && ii > 0)
            {
                particles.overwrite_part(bmax[ibin] - 1, iPart);
                bmax[ibin]--;
                ii--;
                iPart = indexs_to_erase[ii];
            }
            if (iPart >= bmin[ibin] && iPart < bmax[ibin])
            { //On traite la dernière particule (qui peut aussi etre la premiere)
                particles.overwrite_part(bmax[ibin] - 1, iPart);
                bmax[ibin]--;
            }
        }
    }
    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for (int unsigned ibin = 1; ibin < bmax.size(); ibin++)
    {                                             //First bin don't need to be shifted
        ii = bmin[ibin] - bmax[ibin - 1];         // Shift the bin in memory by ii slots.
        iPart = min(ii, bmax[ibin] - bmin[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
        if (iPart > 0)
            particles.overwrite_part(bmax[ibin] - iPart, bmax[ibin - 1], iPart);
        bmax[ibin] -= ii;
        bmin[ibin] = bmax[ibin - 1];
    }

    // Delete useless Particles
    //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
    //Nevertheless, you might want to free memory and have the actual number of particles
    //really equal to the size of the vector. So we do:
    particles.erase_particle_trail(bmax.back());
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}