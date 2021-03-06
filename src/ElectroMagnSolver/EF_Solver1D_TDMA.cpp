
#include "EF_Solver1D_TDMA.h"

#include "ElectroMagn.h"
#include "Field1D.h"

EF_Solver1D_TDMA::EF_Solver1D_TDMA(PicParams &params, int nx_sou_left)
    : Solver1D(params)
{
    dx = params.cell_length[0];
    dx_inv_ = 1.0 / dx;
    dx_sq = dx * dx;
    nx = params.n_space_global[0]+1;
    nx_source_left = nx_sou_left;

    if(params.bc_em_type_x[0] == "Dirichlet"){
        bc_x_left = 1;
        bc_e_value[0][0] = params.bc_em_value_x[0];
    }
    else if(params.bc_em_type_x[0] == "Neumann"){
        bc_x_left = 2;
        bc_e_derivative[0][0] = params.bc_em_value_x[0];
    }

    if(params.bc_em_type_x[1] == "Dirichlet"){
        bc_x_right = 1;
        bc_e_value[0][1] = params.bc_em_value_x[1];
    }
    else if(params.bc_em_type_x[1] == "Neumann"){
        bc_x_right = 2;
        bc_e_derivative[0][1] = params.bc_em_value_x[1];
    }

    initTDMA();

}

EF_Solver1D_TDMA::~EF_Solver1D_TDMA()
{
}

void EF_Solver1D_TDMA::operator()( ElectroMagn* fields)
{
    Field1D* Ex1D           = static_cast<Field1D*>(fields->Ex_);
    Field1D* rho1D          = static_cast<Field1D*>(fields->rho_);
    Field1D* phi1D          = static_cast<Field1D*>(fields->phi_);

    solve_TDMA(rho1D, phi1D);
    solve_Ex(phi1D, Ex1D);
}


void EF_Solver1D_TDMA::initTDMA()
{
    a = new double[nx - nx_source_left];
    b = new double[nx - nx_source_left];
    c = new double[nx - nx_source_left];
    f = new double[nx - nx_source_left];
    e = new double[nx - nx_source_left];
    d = new double[nx - nx_source_left];
    for(int i =1; i < nx-1-nx_source_left; i++)
    {
        a[i] = 1.0;
        b[i] = -2.0;
        c[i] = 1.0;
    }

    if(bc_x_left == 1){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = 0.0;
    }
    else if(bc_x_left == 2){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = -1.0;
    }

    if(bc_x_right == 1){
        a[nx-1-nx_source_left] = 0.0;
        b[nx-1-nx_source_left] = 1.0;
        c[nx-1-nx_source_left] = 0.0;
    }
    else if(bc_x_right == 2){
        a[nx-1-nx_source_left] = -1.0;
        b[nx-1-nx_source_left] = 1.0;
        c[nx-1-nx_source_left] = 0.0;
    }



}


void EF_Solver1D_TDMA::solve_TDMA(Field* rho, Field* phi)
{
    Field1D* rho1D = static_cast<Field1D*>(rho);
    Field1D* phi1D = static_cast<Field1D*>(phi);

    //> The boundary value can be changed with time
    for(int i = 1; i < nx-1-nx_source_left; i++)
    {
        f[i] = -dx_sq * const_ephi0_inv * (*rho1D)(i+nx_source_left);
    }

    if(bc_x_left == 1){
        f[0] = bc_e_value[0][0];
    }
    else if(bc_x_left == 2){
        f[0] = bc_e_derivative[0][0];
    }

    if(bc_x_right == 1){
        f[nx-1-nx_source_left] = bc_e_value[0][1];
    }
    else if(bc_x_right == 2){
        f[nx-1-nx_source_left] = -bc_e_derivative[0][1];
    }

    e[0] = c[0] / b[0];
    d[0] = f[0] / b[0];
    for(int i =1; i < nx-1-nx_source_left; i++)
    {
        e[i] = c[i] / ( b[i] - a[i] * e[i-1] );
        d[i] = ( f[i] -a[i] * d[i-1] ) / ( b[i] - a[i] * e[i-1] );
    }

    (*phi1D)(nx-1) = ( f[nx-1-nx_source_left] - a[nx-1-nx_source_left] * d[nx-2-nx_source_left] ) / ( b[nx-1-nx_source_left] - a[nx-1-nx_source_left] * e[nx-2-nx_source_left] );
    for(int i = nx-2-nx_source_left; i >= 0; i--)
    {
        (*phi1D)(i+nx_source_left) = d[i] - e[i] * (*phi1D)(i+1+nx_source_left);
    }

    for(int i=0; i<nx_source_left; i++)
    {
        (*phi1D)(i) = (*phi1D)(nx_source_left);
    }


}


void EF_Solver1D_TDMA::solve_Ex(Field* phi, Field* Ex)
{
    Field1D* phi1D = static_cast<Field1D*>(phi);
    Field1D* Ex1D = static_cast<Field1D*>(Ex);


    for(int i = 1; i < nx-1; i++)
    {
        (*Ex1D)(i) = - ((*phi1D)(i+1) - (*phi1D)(i-1)) *0.5 * dx_inv_;
    }

    (*Ex1D)(0) = -(-3.0 * (*phi1D)(0) + 4.0 * (*phi1D)(1) - (*phi1D)(2)) *0.5 * dx_inv_;
    (*Ex1D)(nx-1) = -((*phi1D)(nx-3) - 4.0 * (*phi1D)(nx-2) + 3.0 * (*phi1D)(nx-1)) *0.5 * dx_inv_;

}


// no source region for electric field
void EF_Solver1D_TDMA::initTDMA_org()
{
    a = new double[nx];
    b = new double[nx];
    c = new double[nx];
    f = new double[nx];
    e = new double[nx];
    d = new double[nx];
    for(int i =1; i < nx-1; i++)
    {
        a[i] = 1.0;
        b[i] = -2.0;
        c[i] = 1.0;
    }
    if(bc_x_left == 1){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = 0.0;
    }
    else if(bc_x_left == 2){
        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = -1.0;
    }
    if(bc_x_right == 1){
        a[nx-1] = 0.0;
        b[nx-1] = 1.0;
        c[nx-1] = 0.0;
    }
    else if(bc_x_right == 2){
        a[nx-1] = -1.0;
        b[nx-1] = 1.0;
        c[nx-1] = 0.0;
    }
}

// no source region for electric field
void EF_Solver1D_TDMA::solve_TDMA_org(Field* rho, Field* phi)
{
    Field1D* rho1D = static_cast<Field1D*>(rho);
    Field1D* phi1D = static_cast<Field1D*>(phi);
    //> The boundary value can be changed with time
    for(int i = 1; i < nx-1; i++)
    {
        f[i] = -dx_sq * const_ephi0_inv * (*rho1D)(i);
    }
    if(bc_x_left == 1){
        f[0] = bc_e_value[0][0];
    }
    else if(bc_x_left == 2){
        f[0] = bc_e_derivative[0][0];
    }
    if(bc_x_right == 1){
        f[nx-1] = bc_e_value[0][1];
    }
    else if(bc_x_right == 2){
        f[nx-1] = -bc_e_derivative[0][1];
    }
    e[0] = c[0] / b[0];
    d[0] = f[0] / b[0];
    for(int i =1; i < nx-1; i++)
    {
        e[i] = c[i] / ( b[i] - a[i] * e[i-1] );
        d[i] = ( f[i] -a[i] * d[i-1] ) / ( b[i] - a[i] * e[i-1] );
    }
    (*phi1D)(nx-1) = ( f[nx-1] - a[nx-1] * d[nx-2] ) / ( b[nx-1] - a[nx-1] * e[nx-2] );
    for(int i = nx-2; i >= 0; i--)
    {
        (*phi1D)(i) = d[i] - e[i] * (*phi1D)(i+1);
    }
}
