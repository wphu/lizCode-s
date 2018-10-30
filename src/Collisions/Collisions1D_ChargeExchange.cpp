#include "Collisions1D_ChargeExchange.h"
#include "Field2D.h"



#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <sstream>

using namespace std;


// Constructor
Collisions1D_ChargeExchange::Collisions1D_ChargeExchange(PicParams& params, vector<Species*>& vecSpecies,
                       unsigned int n_col,
                       vector<unsigned int> sg1,
                       vector<unsigned int> sg2,
                       string CS_fileName)
: Collisions1D(params)
{

    n_collisions    = n_col;
    species_group1  = sg1;
    species_group2  = sg2;
    crossSection_fileName = CS_fileName;




    // Calculate total number of bins
    nbins = vecSpecies[0]->bmin.size();
    totbins = nbins;

    readCrossSection();



}

Collisions1D_ChargeExchange::~Collisions1D_ChargeExchange()
{

}


// Calculates the collisions for a given Collisions1D object
void Collisions1D_ChargeExchange::collide(PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime)
{
    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2;

    vector<int> index1, index2;
    vector<int> n1, n2;
    vector<double> density1, density2;
    double n1_max, n2_max;
    vector<double> temp(3, 0.0);
    int idNew;
    int totNCollision = 0;
    vector<int> bmin1, bmax1, bmin2, bmax2, bmin3, bmax3;
    unsigned int npairs; // number of pairs of macro-particles
    double npairs_double;
    unsigned int i1, i2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, m12, W1, W2;

    double  sigma_cr, sigma_cr_max, v_square, v_magnitude, ke1, ke_primary, ke_secondary,
            ran, P_collision;
    double atomic_mass;

    int iBin_global;
    double ke_radiative;

    Diagnostic1D *diag1D = static_cast<Diagnostic1D*>(diag);

    sg1 = &species_group1;
    sg2 = &species_group2;

    // ions                          atoms
    s1 = vecSpecies[(*sg1)[0]];      s2 = vecSpecies[(*sg2)[0]];
    p1 = &(s1->particles);           p2 = &(s2->particles);
    m1 = s1->species_param.mass;     m2 = s2->species_param.mass;
    W1 = s1->species_param.weight;   W2 = s2->species_param.weight;
    bmin1 = s1->bmin;                bmin2 = s2->bmin;
    bmax1 = s1->bmax;                bmax2 = s2->bmax;

    atomic_mass = s1->species_param.atomic_mass;

    // Loop on bins
    n1.resize(bmin1.size());
    density1.resize(bmin1.size());
    n1_max = 0.0;
    for(unsigned int ibin=0 ; ibin<nbins ; ibin++)
    {
        n1[ibin] = bmax1[ibin] - bmin1[ibin];
    }

    n2.resize(bmin2.size());
    density2.resize(bmin2.size());
    n2_max = 0.0;
    for(unsigned int ibin=0 ; ibin<nbins ; ibin++)
    {
        n2[ibin] = bmax2[ibin] - bmin2[ibin];
        density2[ibin] = n2[ibin] * W2;
        if( density2[ibin] > n2_max ) { n2_max = density2[ibin]; };
    }

    sigma_cr_max = maxCV(s1, s2);
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {

        if(  (ibin+1) * params.cell_length[0] < params.region_collision_zoom[0]
          || ibin * params.cell_length[0] > params.region_collision_zoom[1] )
        {
            collision_zoom_factor = 1.0;
        }
        else
        {
            collision_zoom_factor = params.collision_zoom_factor;
        }

        //>calculate the particle number of species1 in each cell, and the indexs of particles in the cell
        index1.resize( n1[ibin] );

        for(int iPart = 0; iPart < n1[ibin]; iPart++)
        {
            index1[iPart] = bmin1[ibin] + iPart;
        }
        random_shuffle(index1.begin(), index1.end());

        //>calculate the particle number of species2 in each cell, and the indexs of particles in the cell
        index2.resize( n2[ibin] );

        for(int iPart = 0; iPart < n2[ibin]; iPart++)
        {
            index2[iPart] = bmin2[ibin] + iPart;
        }
        random_shuffle(index2.begin(), index2.end());

        // Now start the real loop
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------

        npairs_double = n1[ibin] * (1 - exp(-density2[ibin] * sigma_cr_max * timesteps_collision * timestep * collision_zoom_factor) );
        npairs = npairs_double;
        npairsRem[ibin] += ( npairs_double - npairs );
        if(npairsRem[ibin] >= 1.0)
        {
            npairsRem[ibin] = npairsRem[ibin] - 1.0;
            npairs++;
        }

        if(npairs > n1[ibin])
        {
            cout<<"npairs is larger than the particle number in a cell!!!"<<endl;
            cout<<"npairs, n1 are: "<<npairs<<" "<<n1[ibin]<<endl;
            npairs = n1[ibin];
        }
        if(npairs > n2[ibin])
        {
            cout<<"npairs is larger than the particle number in a cell!!!"<<endl;
            cout<<"npairs, n2 are: "<<npairs<<" "<<n2[ibin]<<endl;
            npairs = n2[ibin];
        }

        for(int i = 0; i < npairs; i++)
        {
            i1 = index1[i];
            i2 = index2[i];

            v_square = (p1->momentum(0,i1) - p2->momentum(0,i2) ) * ( p1->momentum(0,i1) - p2->momentum(0,i2) ) +
                       (p1->momentum(1,i1) - p2->momentum(1,i2) ) * ( p1->momentum(1,i1) - p2->momentum(1,i2) ) +
                       (p1->momentum(2,i1) - p2->momentum(2,i2) ) * ( p1->momentum(2,i1) - p2->momentum(2,i2) );
            v_magnitude = sqrt(v_square);
            // kinetic energy of species1 (ion)
            ke_primary = 0.5 * m1 * v_square / const_e;
            ke1 =  ke_primary;

            sigma_cr = v_magnitude * interpCrossSection(ke1);
            P_collision = sigma_cr / sigma_cr_max;

            v_square = p1->momentum(0,i1) * p1->momentum(0,i1) +
                       p1->momentum(1,i1) * p1->momentum(1,i1) +
                       p1->momentum(2,i1) * p1->momentum(2,i1);
            v_magnitude = sqrt(v_square);
            // kinetic energy of species1 (ion)
            ke_primary = 0.5 * m1 * v_square / const_e;

            // Generate a random number between 0 and 1
            double ran_p = (double)rand() / RAND_MAX;

            if(ran_p < P_collision)
            {
                calculate_scatter_velocity(s1, s2, i1, i2);
                temp[0] = p2->momentum(0, i2);
                temp[1] = p2->momentum(1, i2);
                temp[2] = p2->momentum(2, i2);
                p2->momentum(0, i2) = p1->momentum(0, i1);
                p2->momentum(1, i2) = p1->momentum(1, i1);
                p2->momentum(2, i2) = p1->momentum(2, i1);
                p1->momentum(0,i1) = temp[0];
                p1->momentum(1,i1) = temp[1];
                p1->momentum(2,i1) = temp[2];

                temp[0] = p2->position(0, i2);
                p2->position(0, i2) = p1->position(0, i1);
                p1->position(0,i1) = temp[0];

                temp[0] = p2->position_old(0, i2);
                p2->position_old(0, i2) = p1->position_old(0, i1);
                p1->position_old(0,i1) = temp[0];

                totNCollision++;

                v_square = p1->momentum(0,i1) * p1->momentum(0,i1) +
                           p1->momentum(1,i1) * p1->momentum(1,i1) +
                           p1->momentum(2,i1) * p1->momentum(2,i1);
                v_magnitude = sqrt(v_square);
                // kinetic energy of new species1 (ion)
                ke_secondary = 0.5 * m1 * v_square / const_e;
                ke_radiative = ke_primary - ke_secondary;
                iBin_global = ibin;
                diag1D->radiative_energy_collision[n_collisions][iBin_global] += ke_radiative;
            } // end if
        }


    } // end loop on bins
}


double  Collisions1D_ChargeExchange::cross_section(double ke)
{

}

double Collisions1D_ChargeExchange::maxCV(Species* s1, Species* s2)
{
    int nPart1, nPart2, npairs;
    Particles *p1, *p2;

    p1 = &(s1->particles);
    p2 = &(s2->particles);
    nPart1 = p1->size();
    nPart2 = p2->size();

    double v_square;
    double v_magnitude;
    double mass;
    double atomic_mass;
    double ke;
    double maxCrossSectionV;
    double crossSectionV;

    maxCrossSectionV = 0.0;
    mass = s1->species_param.mass;
    atomic_mass = s1->species_param.atomic_mass;
    npairs = (nPart1 < nPart2) ? nPart1 : nPart2;
    for(unsigned int iPart = 0; iPart < npairs; iPart++)
    {
        v_square = ( p1->momentum(0,iPart) - p2->momentum(0,iPart) ) * ( p1->momentum(0,iPart) - p2->momentum(0,iPart) ) +
                   ( p1->momentum(1,iPart) - p2->momentum(1,iPart) ) * ( p1->momentum(1,iPart) - p2->momentum(1,iPart) ) +
                   ( p1->momentum(2,iPart) - p2->momentum(2,iPart) ) * ( p1->momentum(2,iPart) - p2->momentum(2,iPart) );
        v_magnitude = sqrt(v_square);
        ke = 0.5 * mass * v_square  / const_e;
        crossSectionV = v_magnitude * interpCrossSection(ke);
        if(crossSectionV > maxCrossSectionV) {maxCrossSectionV = crossSectionV;};
    }
    return maxCrossSectionV;
}


// elastic collision is assumed, same as Collisions1D_DSMC::ELASTIC
void Collisions1D_ChargeExchange::calculate_scatter_velocity(Species* s1, Species* s2, int i1, int i2)
{
    double VRCP[3];		//VRCP(3) are the post-collision components of the relative velocity
    double VCCM[3];		//VCCM(3) are the components of the centre of mass velocity
    double RML = 0.5;
    double RMM = 0.5;
    double A,B,C,D;
    double OC,SC;
    double VR;

    Particles *p1, *p2;
    p1 = &(s1->particles);
    p2 = &(s2->particles);

    for (int i=0;i<3;i++)
	{
		VCCM[i] = RML * p1->momentum(i,i1) + RMM * p2->momentum(i,i2);
	}

    VR = sqrt( pow(p1->momentum(0, i1) - p2->momentum(0, i2), 2)
              +pow(p1->momentum(1, i1) - p2->momentum(1, i2), 2)
              +pow(p1->momentum(2, i1) - p2->momentum(2, i2), 2) );

    //use the VHS logic
    double RF = (double)random() / RAND_MAX;
	B = 2.0 * RF - 1.0;	//B is the cosine of a random elevation angle
	A = sqrt(1.0 - B * B);
	VRCP[0] = B * VR;
    
    RF = (double)random() / RAND_MAX;
	C = 2.0 * const_pi * RF; //C is a random azimuth angle
	VRCP[1] = A * cos(C) * VR;
	VRCP[2] = A * sin(C) * VR;

    for (int i = 0; i < 3; i++)
	{
        p1->momentum(i,i1) = VCCM[i] + VRCP[i] * RMM;
        p2->momentum(i,i2) = VCCM[i] - VRCP[i] * RML;
	}
}