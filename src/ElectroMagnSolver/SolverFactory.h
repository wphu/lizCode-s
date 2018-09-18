#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include "EF_Solver1D_TDMA.h"
#include "EF_Solver1D_TDMA_imp.h"
#include "EF_Solver1D_GeneralThomas.h"
#include "MF_Solver1D_Yee.h"
#include "MF_Solver2D_Yee.h"
#include "MF_Solver2D_Cowan.h"
#include "EF_Solver2D_SLU.h"
#include "EF_Solver2D_SLU_DIST.h"
#include "EF_Solver3D_SLU.h"
#include "EF_Solver3D_SLU_DIST.h"
#include "EF_Solver3D_PETSc_KSP.h"
#include "EF_Solver3D_UMFPACK.h"
#include "InputData.h"
#include "PicParams.h"


#include "Tools.h"

class SolverFactory {
public:
    static Solver* create(PicParams& params, InputData &ifile, Grid* grid) {
        Solver* solver = NULL;
        int nx_source_left;
        if ( params.geometry == "1d3v" ) 
        {
            nx_source_left = 0.0; // default
            ifile.extract("nx_source_left", nx_source_left);
            if(params.method == "explicit")
            {
                if(params.solver_type == "GeneralThomas")
                {
                  solver = new EF_Solver1D_GeneralThomas(params, nx_source_left);
                }
                else
                {
                  solver = new EF_Solver1D_TDMA(params, nx_source_left);
                }

            }
            else if(params.method == "implicit")
            {
                solver = new EF_Solver1D_TDMA_imp(params, nx_source_left);
            }

        }
        else if ( params.geometry == "2d3v" ) 
        {
            if(params.poisson_solver == "SuperLU_serial")
            {
                #ifdef SuperLU_serial
                solver = new EF_Solver2D_SLU(params, grid);
                #endif
            }   
        }
        else if ( params.geometry == "3d3v" ) 
        {
            if(params.poisson_solver == "SuperLU_serial")
            {
                #ifdef SuperLU_serial
                solver = new EF_Solver3D_SLU(params, grid);
                #endif
            }
        
            if(params.poisson_solver == "petsc")
            {
                #ifdef petsc
                solver = new EF_Solver3D_PETSc_KSP(params, grid);
                #endif
            }
            if(params.poisson_solver == "umfpack")
            {
                #ifdef umfpack
                solver = new EF_Solver3D_UMFPACK(params, grid);
                #endif
            }
        }
        else 
        {
        }

        if(solver == NULL)
        {
            ERROR("Creating solver failed, poisson_solver is: "<<params.poisson_solver);
        }

        return solver;
    }

};


#endif
