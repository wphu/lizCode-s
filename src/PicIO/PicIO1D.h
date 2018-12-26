#ifndef PICIO1D_H
#define PICIO1D_H

#include <string>
#include <vector>

#include "PicIO.h"
#include "Array4D.h"
#include "Diagnostic1D.h"

using namespace std;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PicIO1D
//  --------------------------------------------------------------------------------------------------------------------
class PicIO1D : public PicIO {
public:
    //! Create // HDF5 environment
    PicIO1D( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies);
    //! Destructor for PicIO
    ~PicIO1D();

    virtual void write( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    //! Build memory and file space for // HDF5 write/read
    void createFieldsPattern( PicParams& params, ElectroMagn* fields );


    // Create particles h5 file pattern
    void createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies );

    // Create particles h5 file pattern
    void createDiagsPattern( PicParams& params, Diagnostic* diag);

    void initVDF( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies );
    // calculate velocity distribution function
    void calVDF( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime);

    // write grid to grid.h5 file
    virtual void writeGrid(Grid* grid){};

    // read grid from grid.h5 file
    virtual void readGrid(Grid* grid){};

private:


};

#endif /* SMILEIO_CART1D_H_ */
