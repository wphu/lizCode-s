#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>

#include "Grid.h"
#include "Particles.h"
#include "PicParams.h"
#include "Pusher.h"
#include "PicParams.h"
#include "Pusher.h"
#include "ElectroMagn.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class Field3D;

struct Particle
{
    double x;
    double vx;
    double vy;
    double vz;
};

//! class Species
class Species
{
  public:
    //! Species creator
    Species(PicParams &, int);

    //! Species destructor
    virtual ~Species();

    //! Species index
    int species_number;

    //! Method returning the Particle list for the considered Species
    inline Particles getParticlesList() const
    {
        return particles;
    }
    inline Particles &getParticlesList()
    {
        return particles;
    }

    //! Method returning the effective number of Particles for the considered Species
    inline unsigned int getNbrOfParticles() const
    {
        return particles.size();
    }
    // capacity() = vect ever oversize
    // TO do defince particles.capacity = min.capacity
    inline unsigned int getParticlesCapacity() const
    {
        return particles.capacity();
    }

    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    // only caculate the number density, no electric currents
    virtual void dynamics(double time, unsigned int ispec, ElectroMagn *EMfields, Interpolator *interp,
                          Projector *proj, PicParams &params);

    virtual void dynamics_imp_firstPush(double time, unsigned int ispec, ElectroMagn *EMfields, Interpolator *interp,
                                        Projector *proj, PicParams &params);

    virtual void dynamics_imp_secondPush(double time, unsigned int ispec, ElectroMagn *EMfields, Interpolator *interp,
                                         Projector *proj, PicParams &params);
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    // Caculate the number density and the electric currents
    virtual void dynamics_EM(double time, unsigned int ispec, ElectroMagn *EMfields, Interpolator *interp,
                             Projector *proj, PicParams &params);

    virtual void Project(double time, unsigned int ispec, ElectroMagn *EMfields, Projector *Proj, PicParams &params);

    // absort particles according to the Grid for 2-dimension
    virtual void absorb2D(double time, unsigned int ispec, Grid *grid, PicParams &params);

    //! Method used to initialize the Particle position in a given cell
    void initPosition(unsigned int, unsigned int, double *, unsigned int, std::vector<double>, std::string);

    //! Method used to initialize the Particle 3d momentum in a given cell
    void initMomentum(unsigned int, unsigned int, double *, double *, std::string, std::vector<double> &, PicParams &);

    // Method used to initialize the Particle acceleration for implicit method
    void initAcceleration_imp(unsigned int nPart, unsigned int iPart);

    //! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
    void initWeight(unsigned int, unsigned int, unsigned int, double);

    // Method used to initialize the Particles weight by a constant
    void initWeight_constant(unsigned int nPart, unsigned int ispec, unsigned int iPart, double weight_const);

    //! Method used to initialize the Particle charge
    void initCharge(unsigned int, unsigned int, unsigned int, double);

    //! Method used to heat particles in a given cell
    // Zoom particle velocity while keeping direction not changed
    void heat(unsigned int, unsigned int, unsigned int, double, PicParams &);

    //! Maximum charge at initialization
    double max_charge;

    //! Method used to save all Particles properties for the considered Species
    void dump(std::ofstream &);

    //! Method used to sort particles
    void sort_part();

    //! Vector containing all Particles of the considered Species
    Particles particles;
    //std::vector<int> index_of_particles_to_exchange;

    // particles reaching boundaries, stored in psi_particles to performPSI
    // after performPSI, the particles should be erased
    Particles psi_particles;

    //! to keep rack of ionized electrons
    Species *electron_species;

    //! Cluster width in number of cells
    unsigned int clrw; //Should divide the number of cells in X of a single MPI domain. Should default to 1.
    //! first and last index of each particle bin
    std::vector<int> bmin, bmax;

    //! Oversize (copy from picparams)
    std::vector<unsigned int> oversize;

    //! Cell_length (copy from picparams)
    std::vector<double> cell_length;

    inline void clearExchList()
    {
        indexes_of_particles_to_exchange.clear();
    }
    inline void addPartInExchList(int iPart)
    {
        indexes_of_particles_to_exchange.push_back(iPart);
    }
    std::vector<int> indexes_of_particles_to_exchange;

    //Copy of the species parameters from picparams
    SpeciesStructure species_param;

    //! Method to know if we have to project this species or not.
    //bool  isProj(double time_dual, SimWindow* simWindow);

    double getLostNrjBC() const { return species_param.mass * nrj_bc_lost; }
    double getLostNrjMW() const { return species_param.mass * nrj_mw_lost; }

    double getNewParticlesNRJ() const { return species_param.mass * nrj_new_particles; }
    void reinitDiags()
    {
        nrj_bc_lost = 0;
        nrj_mw_lost = 0;
        nrj_new_particles = 0;
    }

    inline int getMemFootPrint()
    {
        int speciesSize = (2 * ndim + 3 + 1) * sizeof(double) + sizeof(short);
        if (particles.isTestParticles)
            speciesSize += sizeof(unsigned int);
        //speciesSize *= getNbrOfParticles();
        speciesSize *= getParticlesCapacity();
        return speciesSize;
    }

    void printAvgVelocity();

    // record the particles reaching boundary to perform psi
    std::vector<std::vector<int>> indexes_of_particles_to_perform_psi;
    //> iDirection = 0,1,2,3,4,5 ==> west, east, north, south, front, back
    inline void clearPsiList(int iDirection)
    {
        indexes_of_particles_to_perform_psi[iDirection].clear();
    }
    inline void clearPsiList()
    {
        for (int iDirection = 0; iDirection < indexes_of_particles_to_perform_psi.size(); iDirection++)
        {
            indexes_of_particles_to_perform_psi[iDirection].clear();
        }
    }
    inline void addPartInPsiList(int iDirection, int iPart)
    {
        indexes_of_particles_to_perform_psi[iDirection].push_back(iPart);
    }

    // after particle pushing, calculate the depositted charge on the boundary
    void calDepCharge(ElectroMagn *EMfields, PicParams &params);


    void interpolator1D1Order(double *Ex, int i_space, int i_particle, LocalFields &E_local);
    void projector1D1Order(double *rho, int i_space, int i_particle);
    void pusherBoris(int i_space, int i_particle, LocalFields &E_local);
    void pusherBoris(int i_space, int i_particle, LocalFields &E_local, LocalFields &B_local);
    void sort_particles();



    // insert and erase particles for bins: mainly used in Collision and PSI
    void insert_particles_to_bins(Particles &insert_Particles, std::vector<int> &count_in_bins);
    void insert_particles(Particles &insert_particles);
    void erase_particles_from_bins(std::vector<int> &indexs_to_erase);

    std::vector<int> indexes_of_particles_to_absorb;

    vector<vector<Particle>> particles;

    vector<vector<int>> list_removed_particle;

    vector<int> n_particle_in_cell;

    int n_particle_total;
    int n_cell;

  private:
    double q;
    double m;
    double q_over_m;
    double weight;

    //! Method used to apply boundary-condition for the Particles of the considered Species
    PartBoundCond *partBoundCond;

    //! Method used to Push the particles (change momentum & change position)
    Pusher *Push;

    //! Method to create new particles.
    int createParticles(std::vector<unsigned int> n_space_to_create, std::vector<double> cell_index, int new_bin_idx, PicParams &param);
};

#endif
