# Contents

1. [QGDsolver library brief](#QGDsolver-library-brief)
2. [QGDsolver installation](#QGDsolver-installation)
3. [Meeting points for users and developers](#Meeting-points-for-users-and-developers)
4. [Published papers related to QGDsolver technology](#Published-papers-related-to-QGDsolver-technology)
5. [For citation](#For-citation)

# QGDsolver library brief
[To the contents](#Contents)

QGDsolver is OpenFOAM framework for simulation of fluid flows using regularized (QGD/QHD) equations approach. It contains library for approximation of partial derivatives at face centers of unstructured grids and a set of OpenFOAM solvers:

1. **QGDFoam** - solver for compressible viscous perfect gas flows in a wide Mach number range - from 0 to infinity
2. **QHDFoam** - solver for incompressible viscous fluid flows with buoyancy force
3. **particlesQGDFOam** - solver for compressible viscous perfect gas flows in a wide Mach number range with particles - from 0 to infinity
4. **particlesQHDFoam** - solver for incompressible viscous fluid flows with buoyancy force with particles
5. **SRFQHDFoam** - solver for incompressible viscous fluid flows in rotating frame of reference  with buoyancy force
6. **QHDDyMFoam** -  solver for incompressible viscous fluid flows in domains with deforming boundary and with buoyancy force
7. **interQHDFoam** - solver for incompressible 2-phase viscous fluid flows with buoyancy force and surface tension
8. **reactingLagrangianQGDFoam** - solver for reacting multicomponent compressible viscous perfect gas flows in a wide Mach number range with particles - from 0 to infinity
9. **scalarTransportQHDFoam** - solver for scalar transport equation to demonstrate **the very basics** of QGD/QHD equations principles
10. **rhoQGDFoam** - solver for compressible viscous flow with arbitrary equation of state (EoS) and in a wide Mach number range - from 0 to infinity

Brief description of the framework is presented here: https://github.com/unicfdlab/QGDsolver/blob/master/qgd-framework-2020-final.pdf

# QGDsolver installation
[To the contents](#Contents)

The repository is organized as follows:
* *master* branch is used for the Doxygen-generated documents (not finished yet) and last test report
* Next branches of the library correspond to OpenFOAM+ versions as follows:
    - [digitef-dev-1912](https://github.com/unicfdlab/QGDsolver/tree/digitef-dev-v1912) for OpenFOAM+ v1912
    - [digitef-dev-2012](https://github.com/unicfdlab/QGDsolver/tree/digitef-dev-v2012) for OpenFOAM+ v2012
    - [digitef-dev-2112](https://github.com/unicfdlab/QGDsolver/tree/digitef-dev-2112) for OpenFOAM+ v2112
    - [v2212](https://github.com/unicfdlab/QGDsolver/tree/v2212) for OpenFOAM+ v2212
    - [v2312](https://github.com/unicfdlab/QGDsolver/tree/v2312) for OpenFOAM+ v2312
* other branches are for internal use and are not intended for compilation

Source code of releases for OpenFOAM+ are stored in *releases* section, the naming conventions are the same as for repository's branches

To compile sources, run *./Allwmake*

To clean sources, run *./Allwclean*

To change libraries and binaries destination, run script *SwitchDest*: a) *./SwitchDest USER* - this will set installation paths to $FOAM_USER_LIBBIN and $FOAM_USER_APPBIN; b) run *./SwitchDest* - this will set installation paths to $FOAM_LIBBIN and $FOAM_APPBIN


# Meeting points for users and developers
[To the contents](#Contents)

Unfinished, but refining Doxygen documentation: https://unicfdlab.github.io/QGDsolver/html/index.html

In case of questions, please, write to:

* the corresponding www.cfd-online.com Forum threads: https://www.cfd-online.com/Forums/openfoam-news-announcements-other/227336-qgdsolver-openfoam-computational-framework-fluid-flows-based-regularized-equ.html
* Issues of this repository: https://github.com/unicfdlab/QGDsolver/issues
* [The telegram channel](https://t.me/hybridCentralSolvers)
* The [ResearchGate](https://www.researchgate.net/) project dedicated to the development of [QGDsolver library](https://www.researchgate.net/project/QGDsolver-OpenFOAM-framework-for-simulation-of-fluid-flows-using-regularized-equations-approach)


# Published papers related to QGDsolver technology
[To the contents](#Contents)

| Title | Description |
|------|-------------|
|[Knudsen pump model](https://github.com/tarminik/KnudsenPump): **Tutorial**| ![The Knudsen pump flow field](https://github.com/mkraposhin/QGDsolver/blob/master/photo_5228881387778596247_y.jpg) |
|[Regularized equations of shallow water for inhomogenous and free surface flows in geophysical problems (in Russian: Регуляризованные уравнения мелкой воды для моделирования неоднородных течений и течений со свободной поверхностью в задачах геофизики](https://keldysh.ru/council/3/D00202403/ivanov_av_diss.pdf): **PhD Thesis**| ![Shallow water flow over a cone](https://github.com/mkraposhin/QGDsolver/blob/master/shallow-water-vs-cone.png) |
|[On a new method for regularizing equations two phase incompressible fluid (in Russian)](https://keldysh.ru/papers/2021/prep2021_61.pdf): **Article**|![Filament visualization](https://github.com/unicfdlab/QGDsolver/blob/master/filament-qhd.jpg)|
|[The Eulerian–Lagrangian Approach for the Numerical Investigation of an Acoustic Field Generated by a High-Speed Gas-Droplet Flow](https://www.mdpi.com/2311-5521/6/8/274):  **Article** | ![Jet with particles Logo](https://www.mdpi.com/fluids/fluids-06-00274/article_deploy/html/images/fluids-06-00274-ag-550.jpg)|
|[Simulation of transonic low-Reynolds jets using quasi-gas dynamics equations](https://iopscience.iop.org/article/10.1088/1742-6596/1382/1/012019): **Article**|![QGDFoam vs experiment](https://www.researchgate.net/publication/337709457/figure/fig1/AS:832046201073664@1575386681051/Time-averaged-jet-centreline-Mach-number-distribution-1-QGDFoam-with-i-i14-03-2_W640.jpg)|
|[Prediction of the Free Jet Noise Using Quasi-gas Dynamic Equations and Acoustic Analogy](https://link.springer.com/chapter/10.1007/978-3-030-50436-6_16): **Article**|![QGDFoam instant jet velocities](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-50436-6_16/MediaObjects/500810_1_En_16_Fig5_HTML.png)|
|[Numerical simulation of disk pump problems using OpenFOAM implementation of regularized equations (in Russian)](https://keldysh.ru/papers/2020/prep2020_66.pdf): **Article**|![QHDFoam fluid flow in disk pump](https://github.com/unicfdlab/QGDsolver/blob/master/QHDFoam-diskpump.png)|
|[Numerical modelling of hydrodynamical structures using quasi-gasdynamics algrotithms and its implementation in OpenFOAM (in Russian)](https://keldysh.ru/council/3/D00202403/istomina_diss.pdf):  **PhD Thesis** |![Accreation disk](https://github.com/unicfdlab/PhDTheses/blob/main/Istomina_diss.png)|
|[Biharmonic attractors of internal waves (in Russian)](https://github.com/unicfdlab/PhDTheses/blob/main/ryazanov_da_diss.pdf): **PhD Thesis**|![A scheme of biharmonic attractor computed with QHDFoam](https://github.com/unicfdlab/PhDTheses/blob/main/ryazanov_da.png)|
|[Development of the New OpenFOAM Solver for Shallow Water Simulation Using QGD/QHD Library](https://www.researchgate.net/publication/352222832_Development_of_the_New_OpenFOAM_Solver_for_Shallow_Water_Simulation_Using_QGDQHD_Library): **Presentation**|![Streamlines](https://github.com/unicfdlab/QGDsolver/blob/master/RSWEFoam_zip.png)|


# For citation
[To the contents](#Contents)

If you have found these library/solvers useful, please cite or refer to

* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3878453.svg)](https://doi.org/10.5281/zenodo.3878453)

* M. V. Kraposhin, E. V. Smirnova, T. G. Elizarova, and M. A. Istomina Development of a new OpenFOAM solver using regularized gas dynamic equations //
Computers & Fluids 166, 163–175 (2018). https://doi.org/10.1016/j.compfluid.2018.02.010

* M.  Kraposhin, D. Ryazanov, T. Elizarova Numerical algorithm based on regularized equations for 
incompressible flow modeling and its implementation in OpenFOAM // Computer Physics Communications (2021). https://doi.org/10.1016/j.cpc.2021.108216

