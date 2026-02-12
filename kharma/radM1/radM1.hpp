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
TaskStatus AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain);


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
KOKKOS_INLINE_FUNCTION Real calc_B(Real T) {
    const Real sigma_SB = 5.67e-5; // Stefan-Boltzmann constant probably defined somewhere! (ASK BEN) (CGS)
    return (sigma_SB * T * T * T * T) / M_PI;
}
}
