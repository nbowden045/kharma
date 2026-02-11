/* 
 *  File: radM1.cpp
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
#include "radM1.hpp"


std::shared_ptr<KHARMAPackage> RadM1::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    auto pkg = std::make_shared<KHARMAPackage>("RadM1");
    Params &params = pkg->AllParams();

    // Wind term in funnel
    Real n = pin->GetOrAddReal("RadM1", "ne", 2.e-4);
    params.Add("ne", n);
    Real Tp = pin->GetOrAddReal("RadM1", "Tp", 10.0);
    params.Add("Tp", Tp);
    Real u1 = pin->GetOrAddReal("RadM1", "u1", 0.0);
    params.Add("u1", u1);
    int power = pin->GetOrAddInteger("RadM1", "power", 4);
    params.Add("power", power);
    Real ramp_start = pin->GetOrAddReal("RadM1", "ramp_start", 0.0);
    params.Add("ramp_start", ramp_start);
    Real ramp_end = pin->GetOrAddReal("RadM1", "ramp_end", 0.0);
    params.Add("ramp_end", ramp_end);
    Real const_G0 = pin->GetOrAddReal("RadM1", "const_G0", 1.05e-5);
    params.Add("const_G0", const_G0);
    Real const_G1 = pin->GetOrAddReal("RadM1", "const_G1", 1.05e-5);
    params.Add("const_G1", const_G1);
    Real const_G2 = pin->GetOrAddReal("RadM1", "const_G2", 1.05e-5);
    params.Add("const_G2", const_G2);
    Real const_G3 = pin->GetOrAddReal("RadM1", "const_G3", 1.05e-5);
    params.Add("const_G3", const_G3);

    pkg->AddSource = RadM1::AddSource;

    return pkg;
}


// TaskStatus RadM1::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
// {
//     // Pointers
//     auto pmesh = mdudt->GetMeshPointer();
//     auto pmb0 = mdudt->GetBlockData(0)->GetBlockPointer();
//     // Options
//     const auto& gpars = pmb0->packages.Get("GRMHD")->AllParams();
//     const auto& pars = pmb0->packages.Get("RadM1")->AllParams();
//     const auto& globals = pmb0->packages.Get("Globals")->AllParams();
//     const Real gam = gpars.Get<Real>("gamma");
//     const Real n = pars.Get<Real>("ne");
//     const Real Tp = pars.Get<Real>("Tp");
//     const Real u1 = pars.Get<Real>("u1");
//     const int power = pars.Get<int>("power");
//     const Real ramp_start = pars.Get<Real>("ramp_start");
//     const Real ramp_end = pars.Get<Real>("ramp_end");
//     const Real time = globals.Get<Real>("time");

//     // Pack variables
//     PackIndexMap cons_map;
//     auto dUdt = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
//     const VarMap m_u(cons_map, true);
//     // Get sizes
//     const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
//     const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
//     const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
//     const IndexRange block = IndexRange{0, dUdt.GetDim(5) - 1};

//     // Set the wind via linear ramp-up with time, if enabled
//     const Real current_n = (ramp_end > 0.0) ? m::min(m::max(time - ramp_start, 0.0) / (ramp_end - ramp_start), 1.0) * n : n;

//     pmb0->par_for("add_radM1", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
//         KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
//             const auto& G = dUdt.GetCoords(b);
//             // Need coordinates to evaluate particle addtn rate
//             // Note that makes the wind spherical-only, TODO ensure this
//             GReal Xembed[GR_DIM];
//             G.coord_embed(k, j, i, Loci::center, Xembed);
//             GReal r = Xembed[1], th = Xembed[2];

//             // Particle addition rate: concentrate at poles
//             Real drhopdt = current_n * m::pow(m::cos(th), power) / SQR(1. + r * r);

//             // Insert fluid moving in positive U1, without B field
//             // Ramp up like density, since we're not at a set proportion
//             const Real uvec[NVEC] = {current_n / n * u1, 0, 0};
//             const Real B_P[NVEC] = {0};

//             // Add plasma to the T^t_a component of the stress-energy tensor
//             // Notice that U already contains a factor of sqrt{-g}
//             Real rho_ut, T[GR_DIM];
//             GRMHD::p_to_u_mhd(G, drhopdt, drhopdt * Tp * 3., uvec, B_P, gam, k, j, i, rho_ut, T);

//             dUdt(b, m_u.RHO, k, j, i) += rho_ut;
//             dUdt(b, m_u.UU, k, j, i) += T[0];
//             dUdt(b, m_u.U1, k, j, i) += T[1];
//             dUdt(b, m_u.U2, k, j, i) += T[2];
//             dUdt(b, m_u.U3, k, j, i) += T[3];
//         }
//     );

//     return TaskStatus::complete;
// }

TaskStatus RadM1::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{
    // 1. Pointers and Setup
    auto pmesh = mdudt->GetMeshPointer();
    auto pmb0 = mdudt->GetBlockData(0)->GetBlockPointer();

    // Get Parameters (for your constant G test)
    const auto& pars = pmb0->packages.Get("RadM1")->AllParams();
    const Real G0 = pars.Get<Real>("const_G0");
    const Real G1 = pars.Get<Real>("const_G1");
    const Real G2 = pars.Get<Real>("const_G2");
    const Real G3 = pars.Get<Real>("const_G3");

    // 2. Pack Variables
    // We only need the 'dUdt' (Conserved) array to add our source term to.
    PackIndexMap cons_map;
    auto dUdt = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    const VarMap m_u(cons_map, true);

    // 3. Get Loop Bounds
    const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, dUdt.GetDim(5) - 1};

    // 4. Parallel Loop
    pmb0->par_for("add_radM1_source", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            
            // A. Get Geometry
            const auto& G = dUdt.GetCoords(b);
            Real gdet = G.gdet(Loci::center, j, i); // sqrt(-g)

            // B. Calculate G_nu (Covariant Components)
            // In the future, you will pass (rho, u, B, R_mn) here.
            Real Gnu_lower[GR_DIM];
            calc_Gnu(G0, G1, G2, G3, Gnu_lower);

            // C. Apply Source to Gas
            // The gas evolution equation is: partial_t( U ) + ... = sqrt(-g) * G_nu
            // Note: U1, U2, U3 in KHARMA correspond to T^t_1, T^t_2, T^t_3 (covariant spatial index)
            // UU corresponds to the energy equation (T^t_t or T^t_0)
            
            // Energy (Component 0)
            dUdt(b, m_u.UU, k, j, i) += gdet * Gnu_lower[0];

            // Momentum (Components 1, 2, 3)
            dUdt(b, m_u.U1, k, j, i) += gdet * Gnu_lower[1];
            dUdt(b, m_u.U2, k, j, i) += gdet * Gnu_lower[2];
            dUdt(b, m_u.U3, k, j, i) += gdet * Gnu_lower[3];
        }
    );

    return TaskStatus::complete;
}