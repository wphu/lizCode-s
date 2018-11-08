#include "Grid2D.h"
#include "Field2D.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <fstream>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creators for Grid2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Grid2D::Grid2D() : Grid()
{


}

// with the dimensions as input argument
Grid2D::Grid2D(
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
    gapKind = gap_kind;
    ny_source = ny_source_temp;
    ny_gapHeight = ny_gapHeight_temp;
    nx_gapWeight = nx_gapWeight_temp;
    ny_bevel_depth = ny_bevel_depth_temp;
    potential_wall = potential_wall_temp;

    dx = params.cell_length[0];
    dy = params.cell_length[1];

    bevel_depth = ny_bevel_depth * dy;

    // number of nodes of the grid in the x-direction
    dims_.resize(2);
    dims_.resize(2);

    dims_[0] = params.n_space[0]+1;
    dims_[1] = params.n_space[1]+1;


    nx = dims_[0];
    ny = dims_[1];

    allocateDims();
    if(gridType == "rectangle")
    {
            geometry();
    }
    else if(gridType == "gap")
    {
        geometry_gap();
    }
    else if(gridType == "iter_gap")
    {
        geometry_iter_gap();
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Grid2D
// ---------------------------------------------------------------------------------------------------------------------
Grid2D::~Grid2D()
{

}

void Grid2D::compute()
{
    computeNcp();

    n_line = lines.size();
    n_segments.resize(n_line);
    n_segment_total = 0;
    for(int iLine = 0; iLine < n_line; iLine++)
    {
        n_segments[iLine] = lines[iLine].size();
        n_segment_total += n_segments[iLine];
    }
}

void Grid2D::allocateDims( )
{
    iswall_      = new int[dims_[0]*dims_[1]];
    bndr_        = new int[dims_[0]*dims_[1]];
    bndrVal_     = new double[dims_[0]*dims_[1]];
    numcp_       = new int[dims_[0]*dims_[1]];

    iswall_2D    =new int*[dims_[0]];
    for (int i=0; i<dims_[0]; i++)  {
        iswall_2D[i] = iswall_ + i*dims_[1];
        for (int j=0; j<dims_[1]; j++) iswall_2D[i][j] = 0;
    }

    bndr_2D      =new int*[dims_[0]];
    for (int i=0; i<dims_[0]; i++)  {
        bndr_2D[i] = bndr_ + i*dims_[1];
        for (int j=0; j<dims_[1]; j++) bndr_2D[i][j] = 0;
    }

    bndrVal_2D   =new double*[dims_[0]];
    for (int i=0; i<dims_[0]; i++)  {
        bndrVal_2D[i] = bndrVal_ + i*dims_[1];
        for (int j=0; j<dims_[1]; j++) bndrVal_2D[i][j] = 0.0;
    }

    numcp_2D    =new int*[dims_[0]];
    for (int i=0; i<dims_[0]; i++)  {
        numcp_2D[i] = numcp_ + i*dims_[1];
        for (int j=0; j<dims_[1]; j++) numcp_2D[i][j] = 0;
    }

}


// rectangle region
void Grid2D::geometry( )
{

}


// classical gap geometry, with source in y direction, the source region is a rectangle below the upper boundary
void Grid2D::geometry_gap( )
{

}

// iter divetor gap geometry, with bevel top surface in toroidal direction
void Grid2D::geometry_iter_gap( )
{

}


void Grid2D::computeNcp()
{
    ncp=0;
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            if( bndr_2D[i][j]==0 || bndr_2D[i][j]==1
             || bndr_2D[i][j]==2 || bndr_2D[i][j]==8)
            {
                ncp++;
                numcp_2D[i][j]=ncp-1;
            }
        }
    }
}
