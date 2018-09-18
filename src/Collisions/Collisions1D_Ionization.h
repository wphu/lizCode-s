/*
Collisions1D_Ionization class
*/

#ifndef COLLISIONS1D_IONIZATION_H
#define COLLISIONS1D_IONIZATION_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"

using namespace std;

class Particles;

class Collisions1D_Ionization : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Ionization(PicParams&,std::vector<Species*>&,unsigned int,
                            std::vector<unsigned int>,std::vector<unsigned int>,std::vector<unsigned int>,
                            string );
    ~Collisions1D_Ionization();


    double cross_section(double ke);


    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    // get the maximum value of crossSection*velocity
    double maxCV(Particles* particles, double eMass);

private:
    //>the ionization threshold energy
    double energy_ionization_threshold;

};


#endif
