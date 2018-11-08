#include "Grid3D.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <fstream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Grid3D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Grid3D::Grid3D() : Grid()
{


}

// with the dimensions as input argument
Grid3D::Grid3D(
    PicParams &params,
    string grid_type,
    string gap_kind,
    int ny_source_temp,
    int ny_gapHeight_temp,
    int nx_gapWeight_temp,
    int ny_bevel_depth_temp,
    double potential_wall_temp):
    Grid(params)
{
    gridType = grid_type;
    potential_wall = potential_wall_temp;

    // number of nodes of the grid in the x-direction
    dims_.resize(3);

    dims_[0] = params.n_space[0]+1;
    dims_[1] = params.n_space[1]+1;
    dims_[2] = params.n_space[2]+1;


    nx=dims_[0];
    ny=dims_[1];
    nz=dims_[2];

    iswall_3D.allocate_dims(dims_);
    iswall_3D.allocate_dims(dims_);
    bndr_3D.allocate_dims(dims_);
    bndrVal_3D.allocate_dims(dims_);
    numcp_3D.allocate_dims(dims_);

    if(gridType == "rectangle")
    {
            geometry();
    }
    else if(gridType == "gap")
    {
        geometry_gap();
    }

    DEBUG(0, "Grid3D construction function end ==============");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Grid3D
// ---------------------------------------------------------------------------------------------------------------------
Grid3D::~Grid3D()
{

}

void Grid3D::compute()
{
    computeNcp();
}

void Grid3D::allocateDims( )
{
}


//>>>no gap geometry, with source in x direction
void Grid3D::geometry( )
{

}


//>>>classical gap geometry, with source in x direction
void Grid3D::geometry_gap( )
{
}



void Grid3D::computeNcp()
{
    ncp = 0;
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                if( bndr_3D(i,j,k) == 0 || bndr_3D(i,j,k) == 1
                 || bndr_3D(i,j,k) == 2 || bndr_3D(i,j,k) == 8)
                {
                    ncp++;
                    numcp_3D(i,j,k) = ncp-1;
                } 
            }
            
        }
    }

}
