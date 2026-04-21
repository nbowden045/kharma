/* 
 *  File: radM1.hpp
 *  
 *  BSD 3-Clause License
 *  
 *  Copyright (c) 2020, AFD Group at UIUC
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, this
 *     list of conditions and the following disclaimer.
 *  
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once

#include "decs.hpp"

#include "gr_coordinates.hpp"
#include "types.hpp"
#include "kharma_utils.hpp"
#include "grmhd_functions.hpp"

#include <parthenon/parthenon.hpp>


namespace RadM1 {

/**
 * Initialize the radM1 package with several options from the input deck
 */
std::shared_ptr<KHARMAPackage> Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages);

/**
 * Add the radM1 source term.  Applied in Flux::AddSource, just after the FluxDivergence calculation
 */
//TaskStatus AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain);
TaskStatus AddImplicitRadiationSourceTerms(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain);

/**
 * Convert from conserved to primitive variables for the radiation field.
 */
TaskStatus BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse=false);

/**
 * Apply floors to the radiation energy variables. 
 */
void ApplyRadM1Floors(MeshBlockData<Real> *rc, IndexDomain domain);
/*
* These are just place holders to calculate G^\nu following Eq.16 Mckinney et al 2014.
* Should check if it should be G^\nu or G_\nu (ASK BEN).
*/

KOKKOS_INLINE_FUNCTION Real calc_kabs(Real rho, Real T) {
    return 1.0; 
}

// Scattering Opacity (kappa_s)
KOKKOS_INLINE_FUNCTION Real calc_kscattering(Real rho, Real T) {
    return 0.4;
}



// Planck function for blackbody radiation
KOKKOS_INLINE_FUNCTION Real calc_lambda(Real T) {
    const Real sigma_SB = 5.67e-5; // Stefan-Boltzmann constant probably defined somewhere! (ASK BEN) (CGS)
    return (sigma_SB * T * T * T * T) / M_PI;
}


// Calculate radiation four velocity in the lab frame
KOKKOS_INLINE_FUNCTION void calc_4vecs(const GRCoordinates& G, const VariablePack<Real>& P, const VarMap& m,
                                     const int& k, const int& j, const int& i, const Loci loc, FourVectors& D_rad)
{
    const Real gamma = GRMHD::lorentz_calc(G, P, m, k, j, i, loc);
    const Real alpha = 1. / m::sqrt(-G.gcon(loc, j, i, 0, 0));

    D_rad.ucon[0] = gamma / alpha;
    
    VLOOP D_rad.ucon[v+1] = P(m.U1_RAD + v, k, j, i) - gamma * alpha * G.gcon(loc, j, i, 0, v+1);

    G.lower(D_rad.ucon, D_rad.ucov, k, j, i, loc);
}


// Calculate radiation four velocity in the lab frame, with Local variables
template <typename Local>
KOKKOS_INLINE_FUNCTION void calc_4vecs(const GRCoordinates& G, const Local& P, const VarMap& m, const int& j, const int& i, const Loci loc, FourVectors& D_rad)
{
    const Real gamma = GRMHD::lorentz_calc(G, P, m, j, i, loc);
    const Real alpha = 1. / m::sqrt(-G.gcon(loc, j, i, 0, 0));

    D_rad.ucon[0] = gamma / alpha;
    
    VLOOP D_rad.ucon[v+1] = P(m.U1_RAD + v) - gamma * alpha * G.gcon(loc, j, i, 0, v+1);

    G.lower(D_rad.ucon, D_rad.ucov, 0, j, i, loc);
}

// It will calculate the lab frame radiation tensor following
// Equation $$R^{\mu\nu} = \frac{4}{3}E_{rf}u^\mu_{rf}u^\nu_{rf} + \frac{1}{3}E_{rf}g^{\mu\nu}$$
// First index up, second index down. 
// Note that ucon and ucov here are ucon_rad and ucov_rad
KOKKOS_INLINE_FUNCTION void calc_tensor(const Real& UU_rad, const FourVectors& D, const int dir, Real mhd_rad[GR_DIM])
{
   DLOOP1 {
        // mhd_rad[mu] = (4.0/3.0) * UU_rad * D.ucon[dir] * D.ucov[mu] + (1.0/3.0) * UU_rad * (dir == mu ? 1.0 : 0.0);
        mhd_rad[mu] = (4.0/3.0) * UU_rad * D.ucon[dir] * D.ucov[mu] + (1.0/3.0) * UU_rad * (dir == mu ? 1.0 : 0.0);
    }
}

KOKKOS_INLINE_FUNCTION void initialize_radiation_pressure(Real UU, Real * UU_rad) {
    //Here we assume that Pgas + Prad = Ptot
    //This translates to rho * T + 1/3 a_rad * T^4 - Ptot = 0
    //The derivative gives us rho + 4/3 a_rad * T^3 = 0, which we can use to find the root of the equation and solve for T given rho and Ptot. 
    // This should be done if we're simulating high accretion rates, bnecause then we should not start with a low radiation pressure, but for all purposes
    // we are gonna assume here that the radiation pressure is negligible at the start of the simulation, so we can just set it to a small value.

    // radiation pressure is 0.1% of the gas pressure at the start of the simulation.
    *UU_rad = UU * 0.001;

    return;
}
}