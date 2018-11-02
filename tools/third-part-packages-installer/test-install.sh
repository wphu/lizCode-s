install_path_header=/home/huwanpeng/opt
compiler_c=gcc
compiler_cxx=g++
compiler_fortran=gfortran
compiler_mpicc=mpicc
compiler_mpicxx=mpicxx
source_codes_root_path=$(pwd)

# install mpich3
CC=compiler_c
package=mpich-3.2.1
install_path=mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
    #rm -rf ${packages[i]}
fi
export PATH=${install_path_header}/${install_path}/bin:$PATH
which mpicc


# install lapack
FORTRAN=${compiler_fortran}
export FFLAGS="-fPIC"
OPTS="-O2 -frecursive"
DRVOPTS=${OPTS}
NOOPT="-O0 -frecursive"
package=lapack-3.8.0
install_path=lapack
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    mkdir build
    cd build
    mkdir ${install_path_header}/${install_path}
    cmake ../ -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path}
    make
    make install
    cd ../..
fi
export FFLAGS=""

# install umfpack included in SuiteSparse
CC=${compiler_c}
FC=${compiler_fortran}
export CFLAGS="-fPIC"
export CPPFLAGS="-fPIC"
package=SuiteSparse-5.3.0
install_path=SuiteSparse
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mv SuiteSparse ${package}
    cd ${package}
    make BLAS="${install_path_header}/lapack/lib/libblas.a -lgfortran" LAPACK=${install_path_header}/lapack/lib/liblapack.a
    #make BLAS="/home/huwanpeng/source-codes/lapack-3.8.0/librefblas.a -lgfortran" LAPACK=/home/huwanpeng/source-codes/lapack-3.8.0/liblapack.a
    rm -rf ${install_path_header}/${install_path} 
    mkdir ${install_path_header}/${install_path}
    cp -r bin ${install_path_header}/${install_path}/bin
    cp -r lib ${install_path_header}/${install_path}/lib
    cp -r include ${install_path_header}/${install_path}/include
    cd ..
fi