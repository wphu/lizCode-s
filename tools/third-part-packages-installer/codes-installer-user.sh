install_path_header=/home/huwanpeng/opt
compiler_c=gcc
compiler_cxx=g++
compiler_fortran=gfortran

# install mpich3
CC=compiler_c
package=mpich-3.2.1
install_path=mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
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


# install hdf5
CC=compiler_c
CFLAGS=-fPIC
package=hdf5-1.8.20
install_path=hdf5
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
fi


# install hdf5-mpich
CC=mpicc
package=hdf5-1.8.20
install_path=hdf5-mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --enable-parallel --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
fi



# install fftw
CC=gcc
package=fftw-3.3.4
install_path=fftw
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
fi


# install netcdf(netcdf-c)
CC=gcc
export CPPFLAGS="-I${install_path_header}/hdf5/include"
export LDFLAGS="-L${install_path_header}/hdf5/lib"
package=netcdf-4.6.0
install_path=netcdf
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

# install netcdf-cxx
CC=gcc
CXX=g++
export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=netcdf-cxx4-4.3.0
install_path=netcdf
if [ -d ${install_path_header}/${install_path} ];then
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

# install SuperLU-4.3
CC=compiler_c
FC=compiler_fortran
package=superlu_4.3
install_path=SuperLU
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mv SuperLU_4.3 ${package}
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
    chmod a+r -R ${install_path_header}/${install_path}
fi


# install SuperLU-DIST
CC=compiler_c
CXX=compiler_cxx
FC=compiler_fortran
package=superlu_dist_5.3.0
install_path=superlu_dist
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
fi


# install PETSc
CC=compiler_c
FC=compiler_fortran
package=petsc-3.8.3
install_path=petsc
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    python2 ./configure --with-mpi-dir=${install_path_header}/mpich --download-fblaslapack --prefix=${install_path_header}/${install_path}
    make PETSC_DIR=${install_path_header}/${package} PETSC_ARCH=arch-linux2-c-debug all
    make PETSC_DIR=${install_path_header}/${package} PETSC_ARCH=arch-linux2-c-debug install
    cd ..
fi


# sundials
CC=compiler_c
FC=compiler_fortran
package=sundials-3.1.0
install_path=sundials
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mkdir examples
    mkdir lib
    mkdir sundials-build

    cmake \
    -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path} \
    -DEXAMPLES_INSTALL_PATH=${install_path_header}/${install_path}/examples \
    -DCMAKE_LINKER=${install_path_header}/${install_path}/lib \
    -DLAPACK_ENABLE=ON \
    -DOPENMP_ENABLE=ON \
    -DMPI_ENABLE=ON \
    ./${package}

    make
    make install
    rm -rf ${install_path_header}examples ${install_path_header}lib ${install_path_header}sundials-build
fi



# install umfpack included in SuiteSparse
CC=compiler_c
FC=compiler_fortran
package=superlu_4.3
install_path=SuperLU_4.3
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make
    make install
    cd ..
fi

