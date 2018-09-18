/*
 * SmileiIO_Cart3D.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO_Cart3D.h"

#include <sstream>

#include "PicParams.h"
#include "Field3D.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

SmileiIO_Cart3D::SmileiIO_Cart3D(PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies)
: SmileiIO(params)
{
    reloadP(params, vecSpecies);
    createFieldsPattern(params, fields);
}

SmileiIO_Cart3D::~SmileiIO_Cart3D()
{
}

//> create hdf5 file, datespace, dateset and so on
void SmileiIO_Cart3D::createFieldsPattern( PicParams& params, ElectroMagn* fields )
{
    fieldsGroup.dims_global[2] = params.n_space_global[2] + 1;
    fieldsGroup.dims_global[1] = params.n_space_global[1] + 1;
    fieldsGroup.dims_global[0] = params.n_space_global[0] + 1;
 

    fieldsGroup.ndims_[0] = fieldsGroup.dims_global[0];
    fieldsGroup.ndims_[1] = fieldsGroup.dims_global[1];
    fieldsGroup.ndims_[2] = fieldsGroup.dims_global[2];

    fieldsGroup.offset[0] = 0;
    fieldsGroup.offset[1] = 0;
    fieldsGroup.offset[2] = 0;

    fieldsGroup.stride[0] = 1;
    fieldsGroup.stride[1] = 1;
    fieldsGroup.stride[2] = 1;

    fieldsGroup.block[0] = 1;
    fieldsGroup.block[1] = 1;
    fieldsGroup.block[2] = 1;

    // For attribute
    fieldsGroup.aDims = 3;

    createFieldsGroup(fields);


} // END createPattern




void SmileiIO_Cart3D::createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies )
{

}


void SmileiIO_Cart3D::createDiagsPattern(PicParams& params, Diagnostic* diag)
{
    Diagnostic3D* diag3D = static_cast<Diagnostic3D*>(diag);

    string diag_name;
    const char* h5_name;
    hid_t dataset_id;
    int n_dim_data = 3;

    diagsGroup.dataset_stringName.resize(0);

    // ======= set dimension, offset, stride, and block =================
    diagsGroup.dims_global[2] = params.n_space_global[2] + 1;
    diagsGroup.dims_global[1] = params.n_space_global[1] + 1;
    diagsGroup.dims_global[0] = params.n_space_global[0] + 1;
 

    diagsGroup.ndims_[0] = diagsGroup.dims_global[0];
    diagsGroup.ndims_[1] = diagsGroup.dims_global[1];
    diagsGroup.ndims_[2] = diagsGroup.dims_global[2];

    diagsGroup.offset[0] = 0;
    diagsGroup.offset[1] = 0;
    diagsGroup.offset[2] = 0;

    diagsGroup.stride[0] = 1;
    diagsGroup.stride[1] = 1;
    diagsGroup.stride[2] = 1;

    diagsGroup.block[0] = 1;
    diagsGroup.block[1] = 1;
    diagsGroup.block[2] = 1;

    // ======= create diagsGroup ================================
    for(int i_species = 0; i_species < diag3D->n_species; i_species++)
    {
        diagsGroup.dataset_stringName.push_back(diag3D->particleFlux_global[i_species]->name);
        diagsGroup.dataset_data.push_back(diag3D->particleFlux_global[i_species]->data_);

        diagsGroup.dataset_stringName.push_back(diag3D->heatFlux_global[i_species]->name);
        diagsGroup.dataset_data.push_back(diag3D->heatFlux_global[i_species]->data_);
    }
    diagsGroup.dataset_id.resize( diagsGroup.dataset_stringName.size() );

}


void SmileiIO_Cart3D::initVDF( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies )
{

}



void SmileiIO_Cart3D::calVDF( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies)
{
    Species *s;
    Particles *p;
}



//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO_Cart3D::write( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime)
{
    const char* h5_name;
    int iDiag;
    int n_dim_data = 3;
    Diagnostic3D* diag3D = static_cast<Diagnostic3D*>(diag);
    if(params.is_calVDF)
    {
        //calVDF( params, smpi, fields, vecSpecies, itime);
    }

    if( itime % params.dump_step == 0)
    {
        ndims_t = itime / params.dump_step - 1;
        long long ndims_t_temp = ndims_t;

        // create file at current output step
        data_file_name = "data/data" + to_string(ndims_t_temp) + ".h5";
        data_file_id = H5Fcreate( data_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // ============= write attributes, n_dim: number of dimension ======================
        hid_t attrs_dataspace_id, attrs_id;
        int n_dim = 3;
        hsize_t attrs_dims[1];
        attrs_dims[0] = 1;
        attrs_dataspace_id = H5Screate_simple(1, attrs_dims, NULL);
        attrs_id           = H5Acreate2(data_file_id, "n_dim", H5T_STD_I32BE, attrs_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attrs_id, H5T_NATIVE_INT, &n_dim);
        H5Sclose(attrs_dataspace_id);
        H5Aclose(attrs_id);

        // =============write fields============================================
        fieldsGroup.group_id = H5Gcreate(data_file_id, "/Fields", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        for(int i = 0; i < fieldsGroup.dataset_stringName.size(); i++)
        {
            fieldsGroup.dataspace_id = H5Screate_simple(n_dim_data, fieldsGroup.dims_global, NULL);
            h5_name = fieldsGroup.dataset_stringName[i].c_str();
            fieldsGroup.dataset_id[i] = H5Dcreate2(fieldsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, fieldsGroup.dataspace_id,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            fieldsGroup.status = H5Dwrite(fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fieldsGroup.dataset_data[i]);
            fieldsGroup.status = H5Sclose(fieldsGroup.dataspace_id);
            fieldsGroup.status = H5Dclose(fieldsGroup.dataset_id[i]);
        }
        fieldsGroup.status = H5Gclose(fieldsGroup.group_id);


        // =============write Diagnostics ============================================
        diagsGroup.group_id = H5Gcreate(data_file_id, "/Diagnostic", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        for(int i = 0; i < diagsGroup.dataset_stringName.size(); i++)
        {
            diagsGroup.dataspace_id = H5Screate_simple(n_dim_data, diagsGroup.dims_global, NULL);
            h5_name = diagsGroup.dataset_stringName[i].c_str();
            diagsGroup.dataset_id[i] = H5Dcreate2(fieldsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            diagsGroup.status = H5Dwrite(diagsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diagsGroup.dataset_data[i]);
            diagsGroup.status = H5Sclose(diagsGroup.dataspace_id);
            diagsGroup.status = H5Dclose(diagsGroup.dataset_id[i]);
        }
        diagsGroup.status = H5Gclose(diagsGroup.group_id);

        status = H5Fclose(data_file_id);
    }



} // END write


void SmileiIO_Cart3D::writeGrid(Grid* grid)
{

}


void SmileiIO_Cart3D::readGrid(Grid* grid)
{
    hid_t       grid_dataspace_id;
    hid_t       grid_dataset_id;
    herr_t      grid_status;
    string      grid_dataset_name;

    int ii;
    int grid_ndim = 3;
    hsize_t grid_dims_global[3];

    Grid3D* grid3D = static_cast<Grid3D*>(grid);
    grid_file_name  = "data/grid.h5";
    grid_file_id    = H5Fopen( grid_file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    grid_dims_global[0] = grid3D->globalDims_[0];
    grid_dims_global[1] = grid3D->globalDims_[1];
    grid_dims_global[2] = grid3D->globalDims_[2];


    // =============read grid============================================
    grid_dataset_name = "is_wall";
    grid_dataset_id = H5Dopen2(grid_file_id, grid_dataset_name.c_str(), H5P_DEFAULT);
    grid_status = H5Dread(grid_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(grid3D->iswall_global_3D(0,0,0)));
    grid_status = H5Dclose(grid_dataset_id);

    grid_dataset_name = "bndr_type";
    grid_dataset_id = H5Dopen2(grid_file_id, grid_dataset_name.c_str(), H5P_DEFAULT);    
    grid_status = H5Dread(grid_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(grid3D->bndr_global_3D(0,0,0)));
    grid_status = H5Dclose(grid_dataset_id);

    grid_dataset_name = "bndr_val";
    grid_dataset_id = H5Dopen2(grid_file_id, grid_dataset_name.c_str(), H5P_DEFAULT);        
    grid_status = H5Dread(grid_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(grid3D->bndrVal_global_3D(0,0,0)));
    grid_status = H5Dclose(grid_dataset_id);

    grid_status = H5Fclose(grid_file_id);

}