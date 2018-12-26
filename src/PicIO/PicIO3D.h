#ifndef PICIO3D_H
#define PICIO3D_H

#include <string>
#include <vector>

#include "PicIO.h"
#include "Diagnostic3D.h"
#include "Grid3D.h"
#include "Segment.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PicIO3D
//  --------------------------------------------------------------------------------------------------------------------
class PicIO3D : public PicIO {
public:
    //! Create // HDF5 environment
    PicIO3D( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies);
    //! Destructor for PicIO
    ~PicIO3D();

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

#endif /* PICIO3D_H_ */
