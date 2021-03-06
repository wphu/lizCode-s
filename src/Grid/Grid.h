#ifndef GRID_H
#define GRID_H

#include <cmath>

#include <vector>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"
#include "PicParams.h"
#include "Field2D.h"



//! Class Grid: generic class allowing to define complex geometry and boundary, also used
//! to decide if the particle hit the wall
class Grid
{

public:



    //! Constructor for Grid: with no input argument
    Grid(){};

    Grid(PicParams &params){};

    //! Destructor for Grid
    virtual ~Grid() {
        ;
    } ;

    //! Virtual method used to allocate Grid
    virtual void allocateDims(){};
    virtual void geometry(){};
    virtual void computeNcp(){};
    virtual void compute(){};

    //! vector containing the dimensions of the Grid
    //! \todo private/friend/modify SmileiMPI* (JD)
    std::vector<unsigned int> dims_;

    //! returns the dimension of the Grid
    inline std::vector<unsigned int> dims () {return dims_;}
    //! All arrays may be viewed as a 1D array
    //! Linearized diags


    //! pointer to the linearized array
    int* iswall_;
    int* bndr_;
    double* bndrVal_;

    //! The number of the current point in the discrete Poisson Eqution left coefficient matrix
    int* numcp_;
    //>>>total number of numcp_ points
    int ncp;
    std::vector<int> dims_source;
    int nx,ny,nz;
    std::string gridType;


protected:

private:

};

#endif
