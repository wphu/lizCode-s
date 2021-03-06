#ifndef GRIDFACTORY_H
#define GRIDFACTORY_H

#include "Grid.h"
#include "Grid2D.h"
#include "Grid3D.h"
#include "PicParams.h"
#include "InputData.h"
#include "SmileiIO.h"

#include "Tools.h"

#include <string>

using namespace std;

//  --------------------------------------------------------------------------------------------------------------------

//  --------------------------------------------------------------------------------------------------------------------
class GridFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------

    //  --------------------------------------------------------------------------------------------------------------------
    static Grid* create(PicParams& params, InputData &ifile, SmileiIO* sio) {
        Grid* grid = NULL;
        MESSAGE("Geometry:" << params.geometry);
        if ( params.geometry == "1d3v" ) 
        {
        }
        else if ( params.geometry == "2d3v" ) 
        {
            string gridType;
    		string gapKind;
            int ny_source;
            int ny_gapHeight;
            int nx_gapWeight;
            int ny_bevel_depth;
            double potential_wall;

    	    // Loop over each binary Grid group and parse info
    	    unsigned int numGrid=ifile.nComponents("Grid");
    	    for (unsigned int n_Grid = 0; n_Grid < numGrid; n_Grid++) {

                if(n_Grid > 0) {
                    ERROR("Now the number of Grid can only be 1 !!!");
                    return grid;
                }
    			ifile.extract("gridType",gridType,"Grid",n_Grid);

    	        gapKind = "divertor"; // default
    	        ifile.extract("gapKind",gapKind,"Grid",n_Grid);

    	        ny_source = 20; // default
    	        ifile.extract("ny_source",ny_source,"Grid",n_Grid);

    	        ny_gapHeight = 100; // default
    	        ifile.extract("ny_gapHeight",ny_gapHeight,"Grid",n_Grid);

                nx_gapWeight = 100;
    	        ifile.extract("nx_gapWeight",nx_gapWeight,"Grid",n_Grid);

                ny_bevel_depth = 20;
    	        ifile.extract("ny_bevel_depth",ny_bevel_depth,"Grid",n_Grid);

                potential_wall = 0.0;
    	        ifile.extract("potential_wall",potential_wall,"Grid",n_Grid);

                grid = new Grid2D(params, gridType, gapKind, ny_source, ny_gapHeight, nx_gapWeight, ny_bevel_depth, potential_wall);

                if(grid->gridType == "from_file")
                {
                    sio->readGrid(grid);
                }
                grid->compute();
            }
        }
        else if ( params.geometry == "3d3v" ) 
        {
            string gridType;
    		string gapKind;
            int ny_source;
            int ny_gapHeight;
            int nx_gapWeight;
            int ny_bevel_depth;
            double potential_wall;

    	    // Loop over each binary Grid group and parse info
    	    unsigned int numGrid=ifile.nComponents("Grid");
    	    for (unsigned int n_Grid = 0; n_Grid < numGrid; n_Grid++) {

                if(n_Grid > 0) {
                    ERROR("Now the number of Grid can only be 1 !!!");
                    return grid;
                }
    			ifile.extract("gridType",gridType,"Grid",n_Grid);

    	        gapKind = "divertor"; // default
    	        ifile.extract("gapKind",gapKind,"Grid",n_Grid);

    	        ny_source = 20; // default
    	        ifile.extract("ny_source",ny_source,"Grid",n_Grid);

    	        ny_gapHeight = 100; // default
    	        ifile.extract("ny_gapHeight",ny_gapHeight,"Grid",n_Grid);

                nx_gapWeight = 100;
    	        ifile.extract("nx_gapWeight",nx_gapWeight,"Grid",n_Grid);

                ny_bevel_depth = 20;
    	        ifile.extract("ny_bevel_depth",ny_bevel_depth,"Grid",n_Grid);

                potential_wall = 0.0;
    	        ifile.extract("potential_wall",potential_wall,"Grid",n_Grid);

                grid = new Grid3D(params, gridType, gapKind, ny_source, ny_gapHeight, nx_gapWeight, ny_bevel_depth, potential_wall);
    	        if(grid->gridType == "from_file")
                {
                    sio->readGrid(grid);
                    DEBUG(0, "Reading grid ends ================");
                }
                grid->compute();
                DEBUG(1, "Grid compute ends ================");
            }
            //DEBUG(2, "Scatter grid ends ================");
        }
        else 
        {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }
        return grid;
    }

};

#endif
