#ifndef PROJECTORFACTORY_H
#define PROJECTORFACTORY_H

#include "Projector.h"
#include "Projector1D1Order.h"
#include "Projector1D2Order.h"
#include "Projector1D4Order.h"
#include "Projector2D1Order.h"
#include "Projector2D2Order.h"
#include "Projector2D4Order.h"
#include "Projector3D1Order.h"

#include "PicParams.h"
#include "Tools.h"

class ProjectorFactory {
public:
    static Projector* create(PicParams& param) 
    {
        Projector* Proj = NULL;
        // ---------------
        // 1d3v simulation
        // ---------------
        if ( ( params.geometry == "1d3v" ) && ( params.projection_order == (unsigned int)1 ) ) 
        {
            Proj = new Projector1D1Order(params);
        }
        else if ( ( params.geometry == "1d3v" ) && ( params.projection_order == (unsigned int)2 ) ) 
        {
            Proj = new Projector1D2Order(params);
        }

        // ---------------
        // 2d3v simulation
        // ---------------
        else if ( ( params.geometry == "2d3v" ) && ( params.interpolation_order == (unsigned int)1 ) ) 
        {
            Proj = new Projector2D1Order(params);
        }

        // ---------------
        // 3d3v simulation
        // ---------------
        else if ( ( params.geometry == "3d3v" ) && ( params.interpolation_order == (unsigned int)1 ) ) 
        {
            Proj = new Projector3D1Order(params);
        }
        else 
        {
            ERROR( "Unknwon parameters : " << params.geometry << ", Order : " << params.interpolation_order );
        }

        return Proj;
    }

};

#endif
