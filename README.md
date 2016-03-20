# TKDGenerator

## What is this?
A plugin for UG4 to create domains consisting of [truncated octahedrons](https://en.wikipedia.org/wiki/Truncated_octahedron) or Tetra-kai-do-deca-hedron (greek 14-areas, abbr. TKD). 
One can achieve tesselation via stacking TKDs. This is of interest for instance in modeling "stratum corneum", the most outer skin layer,
in which the cells are very flat and have a large overlap like in a brick and mortar model.

## Installation
# Requirements
0. C++ compiler
1. CMake >= 2.8

First install ughub, which is a package manager for the UG4 environment.
Instructions can be found [here](https://github.com/UG4/ughub/wiki).

Now add the source of TKDGenerator to your ug4 installation:
    ughub addsource https://github.com/marscher/tkd_generator
    

Now you are now ready to install/compile ugcore and TKDGenerator plugin from your ug4 home:

    cd $UG4_ROOT
    mkdir build
    cmake -DTKDGenerator=ON ..
    make
