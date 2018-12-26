#ifndef PICIO2D_H
#define PICIO2D_H

#include <string>
#include <vector>

#include "PicIO.h"
#include "Diagnostic2D.h"
#include "Grid2D.h"
#include "Segment.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PicIO2D
//  --------------------------------------------------------------------------------------------------------------------
class PicIO2D : public PicIO {
public:
    //! Create // HDF5 environment
    PicIO2D( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies);
    //! Destructor for PicIO
    ~PicIO2D();

    virtual void write( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime);

    //! Build memory and file space for // HDF5 write/read
    void createFieldsPattern( PicParams& params, ElectroMagn* fields );
    // Create particles h5 file pattern
    void createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies );

    // Create particles h5 file pattern
    void createDiagsPattern( PicParams& params, Diagnostic* diag);

    void initVDF( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies );
    // calculate velocity distribution function
    void calVDF( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies);

    // write grid to grid.h5 file
    virtual void writeGrid(Grid* grid);

    // read grid from grid.h5 file
    virtual void readGrid(Grid* grid);
private:



};

#endif /* SMILEIO_CART2D_H_ */
