# RadM1: General Relativistic Radiation Transport (M1 Closure)

This package implements the **M1 Closure Scheme** for radiation transport within the KHARMA/Parthenon framework. The implementation is based on the formulation described by **Sadowski et al. (2013/2014)**.

## References

The primary reference for the physics and numerical implementation of this module is:

* **Sadowski et al. (2013)**: *"Semi-implicit scheme for treating radiation under M1 closure in general relativistic conservative fluid dynamics codes"* [Sadowski et al. (2013)](https://academic.oup.com/mnras/article/429/4/3533/1020840)

Geometric corrections for coordinate singularities follow the methodology established in:

* **McKinney et al. (2012)**: *"General relativistic magnetohydrodynamic simulations of magnetically choked accretion flows around black holes"* [MNRAS 423, 3083](https://ui.adsabs.harvard.edu/abs/2012MNRAS.423.3083M/abstract)


## Implementation Details

All the modifications in the source code due to radiation has been tagged with the comment *"Out of the package modification RADM1."* such that it can be easily found when navigating through the code. All the radiation specific calculations were put in the radM1 package found in **`radM1.cpp`** and **`radM1.hpp`** files.

### 1. Primitive and Conserved Variables
To support the evolution of the radiation field, we have introduced a new set of variables to the state vectors. These variables are initialized in the RadM1::Initialize function in the radM1.cpp file. These are registered in the `VarMap` structure located in **`types.hpp`**.

* **Primitive Variables:**
    * `UU_RAD`: Radiation energy density in the fluid frame ($\hat{E}$).
    * `U1_RAD`, `U2_RAD`, `U3_RAD`: Radiation fluid-frame flux components ($\hat{F}^i$).
* **Conserved Variables:**
    * `UU_RAD` (Conserved): Time-component of the radiation stress-energy tensor ($\sqrt{-g} R^t_0$).
    * `U1_RAD`, `U2_RAD`, `U3_RAD`: Spatial components of the radiation stress-energy tensor ($\sqrt{-g} R^t_i$).

### 2. Correction of the Connection Coeficients
When calculating the evolution of the electromagnetic and radiation tensor, the terms involving the Christoffel symbols calculated at cell centres will not balance out with the corresponding spatial derivatives and it can lead to catastrophic errors. To deal with that, we follow the approach in Mckinney et al. (2012), Appendix A. We recalculated the christoffel symbols to match the derivatives. This is done in gr_coordinates.cpp and therefore, the parameter `correct_connections = true` is necessary when running RadM1.

### 3. Calculating the flux of radiation variables

The evolution of the radiation variables is performed as:

$$\partial_t (\sqrt{-g} R^t_\nu) + \partial_i (\sqrt{-g} R^i_\nu) = \sqrt{-g} R^\kappa_\lambda \Gamma^\lambda_{\nu \kappa} - \sqrt{-g} G_\nu$$

The flux evolution of all the variables is calculated in **`get_flux.hpp`**. The auxiliary functions used are defined either in **`flux.cpp`** or **`flux_functions.cpp`**. For the evolution of the radiative variables, we added the calculation of the geometric source term $\sqrt{-g} R^\kappa_\lambda \Gamma^\lambda_{\nu \kappa}$ into the function `Flux::AddGeoSource`.

After calculating the geometric source term, we add the radiation four-force ($G_\nu$) as a source term as well into the function `RadM1::AddSource`.


We can calculate the characteristic wavespeed for the GRMHD and radiation variables separately and use it for evaluating its corresponding advective fluxes. To account for fluxes of nearby cells, we calculate

$$ \frac{dU}{dt} = - \frac{\mathcal{F}_R^1 - \mathcal{F}_L^1}{dx^1} - \frac{\mathcal{F}_R^2 - \mathcal{F}_L^2}{dx^2} -\frac{\mathcal{F}_R^3 - \mathcal{F}_L^3}{dx^3}$$

This is done inside **`get_flux.hpp`** in `GetFlux` function. The characteristic wavespeed of the radiation is calculated in `vchar_rad` in **`flux_functions.hpp`**.

We then set the timestep of the simulation to be defined by the radiation characteristic wavespeed into `EstimateTimestep` found in **`grmhd.cpp`**.


### 4. Setting Initial Conditions for Fishbone-Montecrief Torus