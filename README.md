QGDsolver is OpenFOAM framework for simulation of fluid flows using regularized equations approach. It contains library for approximation of partial derivatives at face centers of unstructured grids and a set of OpenFOAM solvers for:

1. QGDFoam - solver for compressible viscous perfect gas flows in a wide Mach number range - from 0 to infinity
2. QHDFoam - solver for incompressible viscous fluid flows with buoyancy force
3. particlesQGDFOam - solver for compressible viscous perfect gas flows in a wide Mach number range with particles - from 0 to infinity
4. particlesQHDFoam - solver for incompressible viscous fluid flows with buoyancy force with particles
5. SRFQHDFoam - solver for incompressible viscous fluid flows in rotating frame of reference  with buoyancy force
6. QHDDyMFoam -  solver for incompressible viscous fluid flows in domains with deforming boundary and with buoyancy force
7. interQHDFoam - solver for incompressible 2-phase viscous fluid flows with buoyancy force and surface tension
8. reactingLagrangianQGDFoam - solver for reacting multicomponent compressible viscous perfect gas flows in a wide Mach number range with particles - from 0 to infinity
9. scalarTransportQHDFoam - solver for scalar transport equation to demonstrate **the very basics** of QGD/QHD equations principles



If you found these solver/solvers to be useful or you want to read about implemented numerical algorithms, please cite or refer to

* M. V. Kraposhin, E. V. Smirnova, T. G. Elizarova, and M. A. Istomina, 
Computers & Fluids 166, 163–175 (2018). https://doi.org/10.1016/j.compfluid.2018.02.010

* M.  Kraposhin, D. Ryazanov, T. Elizarova, I. Sibgatullin, M. Kalugin, V. Velikhov, Ey. Ryabinkin
OpenFOAM High Performance Computing Solver for Simulation of Internal Wave Attractors in Stratified
Flows Using Regularized Hydrodynamic Equations // in 2018 Ivannikov Ispras Open Conference (ISPRAS) 
Proceedings, 2018, https://doi.org/10.1109/ISPRAS.2018.00027 

