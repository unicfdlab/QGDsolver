## Building

This branch is laid out as an OpenFOAM (openfoam.com) plugin.

- Standalone build into user locations (`$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN`):
  ```bash
  ./Allwmake
  ```
- As a plugin: place (or add as a git submodule) this repository in
  `$WM_PROJECT_DIR/plugins/QGDsolver` and run `./Allwmake-plugins` from
  `$WM_PROJECT_DIR`. It is then installed into `$FOAM_APPBIN`, `$FOAM_LIBBIN`.
- Custom location:
  `./Allwmake -prefix=user|group|openfoam|/some/path`

If you found these solver/solvers useful or you want to read about implemented numerical algorithms, please cite or refer to

* M. V. Kraposhin, E. V. Smirnova, T. G. Elizarova, and M. A. Istomina, 
Computers & Fluids 166, 163–175 (2018). https://doi.org/10.1016/j.compfluid.2018.02.010

* M.  Kraposhin, D. Ryazanov, T. Elizarova Numerical algorithm based on regularized equations for 
incompressible flow modeling and its implementation in OpenFOAM // Computer Physics Communications (2021). https://doi.org/10.1016/j.cpc.2021.108216
