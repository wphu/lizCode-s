#ifndef SMILEIIOFACTORY_H
#define SMILEIIOFACTORY_H

#include "SmileiIO.h"
#include "SmileiIO_Cart1D.h"

#include "PicParams.h"
#include "Diagnostic.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIOFactory
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIOFactory {
public:
    static SmileiIO* create(PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies) 
    {
        SmileiIO* sio = NULL;
        if ( params.geometry == "1d3v" ) 
        {
            sio = new  SmileiIO_Cart1D(params, fields, vecSpecies);
        }
        else if ( params.geometry == "2d3v" ) 
        {
            sio = new  SmileiIO_Cart2D(params, fields, vecSpecies);
        }
        else if ( params.geometry == "3d3v" ) 
        {
            sio = new  SmileiIO_Cart3D(params, fields, vecSpecies);
        }
        else 
        {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }

        return sio;
    }

};

#endif
