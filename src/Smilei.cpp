////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                                                ////
////                                                                                                                ////
////                                   PARTICLE-IN-CELL CODE SMILEI                                                 ////
////                    Simulation of Matter Irradiated by Laser at Extreme Intensity                               ////
////                                                                                                                ////
////                          Cooperative OpenSource Object-Oriented Project                                        ////
////                                      from the Plateau de Saclay                                                ////
////                                          started January 2013                                                  ////
////                                                                                                                ////
////                                                                                                                ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Smilei.h"

#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <iomanip>

#include "InputData.h"
#include "PicParams.h"
#include "GridFactory.h"
#include "SpeciesFactory.h"
#include "PartSourceFactory.h"
#include "CollisionsFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "SolverFactory.h"
#include "SmileiIOFactory.h"
#include "PSIFactory.h"
#include "ElectroMagnBC_Factory.h"
#include "DiagnosticFactory.h"
#include "Timer.h"


#ifdef  __DEBUG
#include <gperftools/profiler.h>
#endif

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    //cout.setf( ios::fixed,  ios::floatfield ); // floatfield set to fixed

    // -------------------------
    // Simulation Initialization
    // -------------------------

    // Check for namelists (input files)
    vector<string> namelists(argv + 1, argv + argc);
    if (namelists.size()==0) ERROR("No namelists given!");



    TITLE("Input data info");
    // Read the namelists file (no check!)
    InputData input_data(namelists);
    
    TITLE("PicParams info");
    // Read simulation & diagnostics parameters
    PicParams params(input_data);
    params.print();

    // Count timer
    vector<Timer> timer(13);
    timer[0].init("Total time");
    timer[1].init("EmitLoad");
    timer[2].init("Collide");
    timer[3].init("Interpolate and Move");
    
    timer[5].init("Absorb paritcle (2D)");
    timer[6].init("Project Particle");
    timer[7].init("PSI");
    timer[8].init("Diagnostic");
    timer[9].init("Fields Solve");
    timer[10].init("Write IO");
    timer[11].init("poisson solver init, like: SuperLU factorize, petsc_ksp init");
    timer[12].init("particle sort");

    timer[0].restart();


    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------

    // ---------------------------
    // Initialize Species & Fields
    // ---------------------------

    TITLE("Initializing particles");
    // Initialize the vecSpecies object containing all information of the different Species
    // ------------------------------------------------------------------------------------
    // vector of Species (virtual)
    vector<Species*> vecSpecies = SpeciesFactory::createVector(params);


    // Initialize the electromagnetic fields and interpolation-projection operators
    // according to the simulation geometry
    // ----------------------------------------------------------------------------

    TITLE("Initializing ElectroMagn Fields");
    // object containing the electromagnetic fields (virtual)
    ElectroMagn* EMfields = ElectroMagnFactory::create(params, input_data);

    TITLE("Initializing Fields Bounary Condition");
    vector<ElectroMagnBC*> vecEmBC = ElectroMagnBCFactory::create(params);

    TITLE("Creating PSI");
    vector<PSI*> vecPSI = PSIFactory::create(params, input_data, vecSpecies);

    TITLE("Creating PartSource");
    vector<PartSource*> vecPartSource = PartSourceFactory::create(params, input_data, vecSpecies);


    // Initialize the collisions (vector of collisions)
    // ------------------------------------------------------------------------------------
    TITLE("Creating Collisions");
    vector<Collisions*> vecCollisions = CollisionsFactory::create(params, input_data, vecSpecies);

    TITLE("Creating Interp/Proj");
    // interpolation operator (virtual)
    Interpolator* Interp = InterpolatorFactory::create(params);

    // projection operator (virtual)
    Projector* Proj = ProjectorFactory::create(params);

    //Create i/o environment
    TITLE("Creating IO output environment");
    SmileiIO*  sio  = SmileiIOFactory::create(params, EMfields, vecSpecies);

    //>Initialize Grid
    TITLE("Creating grid");
    Grid* grid = NULL;
    grid = GridFactory::create(params, input_data, sio);

    if(grid != NULL && grid->gridType != "from_file")
    {
        sio->writeGrid(grid);
    }

    TITLE("Creating Diagnostic");
    Diagnostic*  diag  = DiagnosticFactory::create(params, grid, EMfields, vecPSI, vecCollisions);
    if(params.geometry == "3d3v")
    {
        sio->createDiagsPattern(params, diag);
    }
    


    TITLE("Creating Solver");
    timer[11].restart();
    Solver* solver = NULL;
    solver = SolverFactory::create(params, input_data, grid);
    timer[11].update();

    // ------------------------------------------------------------------------
    // Initialize the simulation times time_prim at n=0 and time_dual at n=+1/2
    // ------------------------------------------------------------------------
    unsigned int stepStart = sio->stepStart, stepStop = params.n_time;
    // time at integer time-steps (primal grid)
    double time_prim = stepStart * params.timestep;
    // time at half-integer time-steps (dual grid)
    double time_dual = (stepStart +0.5) * params.timestep;

    int itime = stepStart;

    // control timesteps, use larger timesteps for some heavy particles
    vector<int> timestep_control;
    timestep_control.resize( params.species_param.size() );
    for(int i = 0; i < params.species_param.size(); i++)
    {
        timestep_control[i] = 0;
    }

    TITLE("Solve the field first time before PIC loop");
    for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
    {
        EMfields->restartRhoJs(ispec, 0);
        vecSpecies[ispec]->Project(time_dual, ispec, EMfields, Proj, params);
    }
    EMfields->computeTotalRhoJ();
    (*solver)(EMfields);


    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    TITLE("Time-Loop is started: number of time-steps n_time = " << params.n_time);

    #ifdef  __DEBUG
    
    string prof_name = "./liz.prof";
    ProfilerStart(prof_name.c_str());
    
    /*
    string prof_name = "./liz.prof";
    ProfilerStart(prof_name.c_str());
    */
    #endif

    if(params.method == "explicit")
    {
        while(itime <= stepStop)
        {
            itime++;
            time_prim += params.timestep;
            time_dual += params.timestep;
            //MESSAGE("timestep = "<<itime);

            // ================== EmitLoad =========================================
            //> add Particle Source: emit from boundary or load in some region
            timer[1].restart();
            for (unsigned int iPS=0 ; iPS<vecPartSource.size(); iPS++)
            {
                vecPartSource[iPS]->emitLoad(params,vecSpecies,itime, EMfields);
            }
            timer[1].update();


            // ================== Collide =========================================
            timer[2].restart();
            if(itime % params.timesteps_collision == 0)
            {
                for (unsigned int icoll=0 ; icoll<vecCollisions.size(); icoll++)
                {
                    vecCollisions[icoll]->collide(params, EMfields, vecSpecies, diag, itime);
                }
            }
            timer[2].update();

            MESSAGE("Interpolate and Move ");
            // ================== Interpolate and Move ===============================
            int tid(0);
            timer[3].restart();
            for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
            {
                timestep_control[ispec]++;
                if(timestep_control[ispec] == params.species_param[ispec].timestep_zoom)
                {
                    vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, params);
                }
            }
            timer[3].update();

            MESSAGE("Sort Particle ");
            // ================== Sort Particle ============================================
            timer[4].restart();
            for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
            {
                if(timestep_control[ispec] == params.species_param[ispec].timestep_zoom)
                {
                    timer[12].restart();
                    vecSpecies[ispec]->sort_part(); // Should we sort test particles ?? (JD)
                    timer[12].update();
                    timestep_control[ispec] = 0;
                }
                vecSpecies[ispec]->clearExchList();
            }
            timer[4].update();

            MESSAGE("Run Diagnostic ");
            // ================== Run Diagnostic =============================================
            // absorb particles and calculate particle flux, heat flux, and average angle for 2D and 3D
            timer[8].restart();
            //diag->run(grid, vecSpecies, EMfields, vecPSI, itime);
            timer[8].update();

            MESSAGE("Project Particle ");
            // ================== Project Particle =========================================
            timer[6].restart();
            for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
            {
                EMfields->restartRhoJs(ispec, 0);
                vecSpecies[ispec]->Project(time_dual, ispec, EMfields, Proj, params);
            }
            timer[6].update();


            MESSAGE("PSI ");
            // ================== Plasma Surface Interacton ==================================
            timer[7].restart();
            for (unsigned int ipsi=0 ; ipsi<vecPSI.size(); ipsi++)
            {
                vecPSI[ipsi]->performPSI(params, grid, vecSpecies, EMfields, diag, itime);
                vecSpecies[vecPSI[ipsi]->species2]->sort_part(); // Should we sort test particles ?? (JD)
            }
            for (unsigned int ispec=0 ; ispec<params.species_param.size(); ispec++)
            {
                // clear psi_particles to avoid unnecessary repeated PSI performs for multiple ion timesteps
                vecSpecies[ispec]->clearPsiList();
                vecSpecies[ispec]->clearExchList(); 
            }
            timer[7].update();

            MESSAGE("Solve Eields");
            // ================== Solve Electromagnetic Fields ===============================
            timer[9].restart();
            EMfields->restartRhoJ();
            EMfields->computeTotalRhoJ();
            (*solver)(EMfields);
            timer[9].update();

            MESSAGE("Write IO");
            // ================== Write IO ====================================================
            timer[10].restart();
            if(params.ntime_step_avg)
            {
                EMfields->incrementAvgFields(itime);
            }
            if(itime % params.dump_step == 0)
            {
                MESSAGE("time step = "<<itime);
            }
            sio->write(params, EMfields, vecSpecies, diag, itime);
            if(itime % params.timesteps_restore == 0)
            {
                sio->storeP(params, vecSpecies, itime);
            }
            timer[10].update();

        }//END of the time loop
    }
    else if(params.method == "implicit")
    {
        while(itime <= stepStop)
        {
            itime++;
            time_prim += params.timestep;
            time_dual += params.timestep;
            //MESSAGE("timestep = "<<itime);

        }//END of the time loop
    }

    #ifdef  __DEBUG
    ProfilerStop();
    #endif


    sio->endStoreP(params, vecSpecies, itime);




    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    timer[0].update();
    for( int i = 0; i < timer.size(); i++)
    {
        if(i == 0)
        {
            TITLE(timer[i].name() << " = " << timer[i].getDateTime());
        }
        else 
        {
            TITLE(timer[i].name() << " = " << timer[i].getTime() <<" s");
        }
    }


    // ------------------------------------------------------------------------
    // check here if we can close the python interpreter
    // ------------------------------------------------------------------------
    TITLE("Cleaning up python runtime environement");
    input_data.cleanup();

    TITLE("Time profiling :");
    //timer[0].update();
    //timer[0].print();


    // ------------------------------------------------------------------
    //                      Temporary validation diagnostics
    // ------------------------------------------------------------------

    // temporary EM fields dump in Fields.h5


    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    delete Proj;
    delete Interp;
    delete EMfields;

    for (unsigned int iPS=0 ; iPS<vecPartSource.size(); iPS++) delete vecPartSource[iPS];
    vecPartSource.clear();

    for(unsigned int i=0; i<vecCollisions.size(); i++) delete vecCollisions[i];
    vecCollisions.clear();

    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    TITLE("END");
    delete sio;
    return 0;

}//END MAIN
