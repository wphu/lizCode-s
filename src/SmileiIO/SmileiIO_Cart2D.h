/*
 * SmileIO_Cart2D.h
 *
 *  Created on: 3 juil. 2013
 */
#ifndef SMILEIO_CART2D_H
#define SMILEIO_CART2D_H

#include <string>
#include <vector>

#include "SmileiIO.h"
#include "Diagnostic2D.h"
#include "Grid2D.h"
#include "Segment.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart2D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart2D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart2D( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies);
    //! Destructor for SmileiIO
    ~SmileiIO_Cart2D();

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
