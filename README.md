QGDsolver is OpenFOAM framework for simulation of fluid flows using regularized equations approach. It contains library for approximation of partial derivatives at face centers of unstructured grids and a set of OpenFOAM solvers:

1. **QGDFoam** - solver for compressible viscous perfect gas flows in a wide Mach number range - from 0 to infinity
2. **QHDFoam** - solver for incompressible viscous fluid flows with buoyancy force
3. **particlesQGDFOam** - solver for compressible viscous perfect gas flows in a wide Mach number range with particles - from 0 to infinity
4. **particlesQHDFoam** - solver for incompressible viscous fluid flows with buoyancy force with particles
5. **SRFQHDFoam** - solver for incompressible viscous fluid flows in rotating frame of reference  with buoyancy force
6. **QHDDyMFoam** -  solver for incompressible viscous fluid flows in domains with deforming boundary and with buoyancy force
7. **interQHDFoam** - solver for incompressible 2-phase viscous fluid flows with buoyancy force and surface tension
8. **reactingLagrangianQGDFoam** - solver for reacting multicomponent compressible viscous perfect gas flows in a wide Mach number range with particles - from 0 to infinity
9. **scalarTransportQHDFoam** - solver for scalar transport equation to demonstrate **the very basics** of QGD/QHD equations principles

Brief description of the framework is presented HERE: https://github.com/unicfdlab/QGDsolver/blob/master/qgd-framework-2020-final.pdf

The repository is organized as follows:
* *master* branch is used for the Doxygen-generated documents (not finished yet) and last test report
* *digitef-dev-ABCD* - latest working version of the framework for OpenFOAM+ version *ABCD*
* other branches are for internal use and are not intended for compilation

Source code of releases for OpenFOAM+ are stored in *releases* section, the naming conventions are the same as for repository's branches

To compile sources, run *./Allwmake*

To clean sources, run *./Allwclean*

To change libraries and binaries destination, run: a) *./SwitchDest USER* - this will set installation paths to $FOAM_USER_LIBBIN and $FOAM_USER_APPBIN; b) run *./SwitchDest* - this will set installation paths to $FOAM_LIBBIN and $FOAM_APPBIN

Unfinished, but refining Doxygen documentation: https://unicfdlab.github.io/QGDsolver/html/index.html

In case of questions, please, write to:

* the corresponding www.cfd-online.com Forum threads: https://www.cfd-online.com/Forums/openfoam-news-announcements-other/227336-qgdsolver-openfoam-computational-framework-fluid-flows-based-regularized-equ.html
* Issues of this repository: https://github.com/unicfdlab/QGDsolver/issues


If you found these solver/solvers to be useful or you want to read about implemented numerical algorithms, please cite or refer to

* M. V. Kraposhin, E. V. Smirnova, T. G. Elizarova, and M. A. Istomina, 
Computers & Fluids 166, 163–175 (2018). https://doi.org/10.1016/j.compfluid.2018.02.010

* M.  Kraposhin, D. Ryazanov, T. Elizarova, I. Sibgatullin, M. Kalugin, V. Velikhov, Ey. Ryabinkin
OpenFOAM High Performance Computing Solver for Simulation of Internal Wave Attractors in Stratified
Flows Using Regularized Hydrodynamic Equations // in 2018 Ivannikov Ispras Open Conference (ISPRAS) 
Proceedings, 2018, https://doi.org/10.1109/ISPRAS.2018.00027 

