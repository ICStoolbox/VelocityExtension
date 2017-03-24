# velext [![Build Status](https://travis-ci.org/ISCDtoolbox/VelocityExtension.svg?branch=master)](https://travis-ci.org/ISCDtoolbox/VelocityExtension)
velext is a program to solve a simplified version of Helmholtz equation ```-h(x) ∆u(x) + u(x) = f(x)```. The physical interpretation of the unknown ```u(x)``` and the function ```f(x)```depends usually on what the equation models (wave propagation, quantum mechanics, etc.). Here, the program is mainly intended to extend a velocity field defined at domain boundaries to interior points, and can be considered as an extrapolation technique. Velext has many applications, one of the main is related to the level set method that requires a velocity field for which the values must be known for at least at mesh vertices that are updated in time.

#### Installation

1. you will need to install the [ISCD Commons Library](https://github.com/ISCDtoolbox/Commons) on your system. 
Please refer to the instructions provided on the ICS Commons Library page in order to install this library.

2. download the zip archive of VelocityExtension or clone this repository:

   ` git clone https://github.com/ISCDtoolbox/VelocityExtension.git `

   navigate to the downloaded directory: 

   ` cd VelocityExtension `

   then create build directory and compile the project with cmake
   ```
   mkdir build
   cd build
   cmake ..
   make
   make install
   ```

#### Usage
After compiling ```velext``` as described above, you should have an executable file in your $HOME/bin directory. If your PATH variable is correctly set to this directory, velext can be called with the following syntax:

    usage: velext [+/-v | -h] [-a val] [-n nit] [-r res] source_file[.mesh] [-c chi_file.[sol]] [-p param_file[.elas]] [-s data_file[.sol]] [-o output_file[.sol]]

The square braces indicate optional arguments. Some commands have flags, some others do not.

The options and flags are:

    --help       show the syntax and exit.
    --version    show the version and date of release and exit.

    -a val       value of diffusion coefficient h
    -n nit       number of iterations max for convergence 
    -r res       value of the residual (Krylov space) for convergence
    -v           suppress any message (for use with function call).
    +v           increase the verbosity level for output.

    source_file.mesh      name of the mesh file
    param_file            name of file containing elasticity parameters
    chi_file              name of file containing characteristic function (LS)
    data_file.sol         name of file containing the initial solution or boundary conditions
    output_file.sol       name of the output file (displacement field)
    
#### Quickstart
You can test the installation and look at examples by entering the [demos](demos) directory and running the program:

    cd demos
    velext disk.mesh -o disk.new.sol        # or equivalently:  velext disk.mesh -p disk.velext -o disk.new.sol

that will produce an output that will look like:

     - VELEXT, Release 2.0a, Jan 19, 2016
       (C) Copyright 2014- , ICS-SU

     - LOADING DATA
        disk.mesh: 1758 vertices, 253 edges, 3301 triangles
        disk.velext: 2 parameters
        Compressing mesh: 908 vertices, 1664 edges, 1664 triangles
     - COMPLETED: 0.006s

     ** MODULE VELEXT: 2.0a
        Matrix and right-hand side assembly
        Solving linear system: 6.553947E-07 in 28 iterations
     ** COMPLETED: 0.018s

     - WRITING DATA
        Uncompressing data: 1758 data vectors
        disk.new.sol: 1758 data vectors
     - COMPLETED: 0.002s

     ** Cumulative time: 0.026s.

#### Authors & contributors
* velext has been initiated by Charles Dapogny (Université Joseph Fourier) and Pascal Frey (Université Pierre et Marie Curie). Current team includes Dena Kazerani (Université Pierre et Marie Curie) and Loic Norgeot (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
velext is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).
