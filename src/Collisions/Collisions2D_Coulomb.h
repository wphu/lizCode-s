/*
Collisions2D_Coulomb class
*/

#ifndef COLLISIONS2D_COULOMB_H
#define COLLISIONS2D_COULOMB_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions2D.h"

class Collisions2D_Coulomb : public Collisions2D
{

public:
    //! Constructor for Collisions2D between two species
    Collisions2D_Coulomb(PicParams&,std::vector<Species*>&,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions2D_Coulomb();


    //! Coulomb logarithm (zero or negative means automatic)
    double coulomb_log;
    //! Contains the debye length in each cluster, computed each timestep
    std::vector<double> debye_length_squared;

    //! Method to calculate the Debye length in each cluster
    void calculate_debye_length(PicParams&,std::vector<Species*>&);

    //! Method called in the main smilei loop to apply collisions at each timestep
    void collide(PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    virtual double cos_chi(double);
private:


};


#endif
