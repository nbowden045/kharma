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
#include "kharma_driver.hpp"
#include "units.hpp"
#include "kharma.hpp"
#include "inverter.hpp"
#include <limits> 
std::shared_ptr<KHARMAPackage> RadM1::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    bool units_enabled = pin->GetOrAddBoolean("units", "on", false);
    bool correct_connections = pin->GetOrAddBoolean("coordinates", "correct_connections", false);
    // Check if the Units package is initialized, since we need it for the radiation four-force calculations.
    // if (!units_enabled) {
    //     printf("\033[1;31mError: Units package not enabled! It must be enabled with/BEFORE RadM1.\033[0m\n");
    //     exit(1);
    // }
    if (!correct_connections) {
        printf("\033[1;33mError: Connection coefficient corrections for GRMHD are disabled. RadM1 requires these connections to evolve the fields properly.\033[0m\n");
        exit(1);
    }

    auto pkg = std::make_shared<KHARMAPackage>("RadM1");
    Params &params = pkg->AllParams();

    //InitializeRadPrims();
    // Adding the conserved and primitive variables for the radiation field.
    auto& driver = packages->Get("Driver")->AllParams();
    auto driver_type = driver.Get<DriverType>("type");

    // I believe this is the boolean to determine whether we are doing an implicit evolution of the radiation field, but I haven't tested it yet.
    // It checks if the driver is imex and if the user has set RadM1/implicit = true in the input file.
    // I think in general we want implicit evolution of the radiation field? Can we have a smarter method than just setting them all to implicit? Maybe check Mckinney et al 2014
    bool implicit_radm1 = (driver_type == DriverType::imex && pin->GetOrAddBoolean("RadM1", "implicit", false));

    // Based on the previous boolean, we set the appropriate evolution flag for the radiation variables.
    MetadataFlag areWeImplicit = (implicit_radm1) ? Metadata::GetUserFlag("Implicit")
                                              : Metadata::GetUserFlag("Explicit");

    //Registers a new label called "RADM1" for the variables to be grouped together later?
    Metadata::AddUserFlag("RADM1");
    // Properties of these new flags
    // I believe Metadata::Cell means these variables are defined at cell centers, but I should ask Cora.
    // We also add the "areWeImplicit" flag, which is either "Implicit" or "Explicit" based on the user's choice in the input file.
    std::vector<MetadataFlag> flags_radm1 = {Metadata::Real, Metadata::Cell, areWeImplicit, Metadata::GetUserFlag("RADM1")};


    //Retrieves the existing flags for the primitive and conserved variables, and adds the new radM1 flags to them.
    //Then adds the new variables for the radiation primitives and conserved variables with these flags.
    // auto flags_prim = driver.Get<std::vector<MetadataFlag>>("prim_flags");
    // flags_prim.insert(flags_prim.end(), flags_radm1.begin(), flags_radm1.end());

    auto flags_prim = driver.Get<std::vector<MetadataFlag>>("prim_flags");
    flags_prim.insert(flags_prim.end(), flags_radm1.begin(), flags_radm1.end());

    // Explicitly tag these as Primitive so your PackVariables call finds them!
    flags_prim.push_back(Metadata::GetUserFlag("Primitive")); 

    // I think this will push this metadata to restart files
    flags_prim.push_back(Metadata::Restart);

    // sync variables across boundaries (ASK CORA)
    if (pin->GetOrAddBoolean("RadM1", "sync_utop_seed", true)) { 
        flags_prim.push_back(Metadata::FillGhost);
    }
    //
    auto flags_cons = driver.Get<std::vector<MetadataFlag>>("cons_flags");
    flags_cons.insert(flags_cons.end(), flags_radm1.begin(), flags_radm1.end());

    auto m_prim_scalar = Metadata(flags_prim);
    pkg->AddField("prims.u_rad", m_prim_scalar);

    auto m_cons_scalar = Metadata(flags_cons);
    pkg->AddField("cons.u_rad", m_cons_scalar);

    auto flags_prim_vec(flags_prim);
    flags_prim_vec.push_back(Metadata::Vector);

    auto flags_cons_vec(flags_cons);
    flags_cons_vec.push_back(Metadata::Vector);

    std::vector<int> s_vector({NVEC});
    auto m_prim_vector = Metadata(flags_prim_vec, s_vector);
    pkg->AddField("prims.uvec_rad", m_prim_vector);

    auto m_cons_vector = Metadata(flags_cons_vec, s_vector);
    pkg->AddField("cons.uvec_rad", m_cons_vector);

    // Get the total number of variables you need to store.
    // This is the sum of the fluid and radiation variables (ask cora how to do it properly without hardcoding)
    int nvar = 12; 

    parthenon::Metadata::AddUserFlag("RadGuessU");
    parthenon::Metadata::AddUserFlag("RadGuessP");
    //Define the custom User Flags
    MetadataFlag flag_guess_u = Metadata::GetUserFlag("RadGuessU");
    MetadataFlag flag_guess_p = Metadata::GetUserFlag("RadGuessP");

    //Create the flags
    std::vector<MetadataFlag> flags_u_guess = {Metadata::Cell, Metadata::Derived, flag_guess_u};
    std::vector<MetadataFlag> flags_p_guess = {Metadata::Cell, Metadata::Derived, flag_guess_p};

    // Add the fields to the package. 
    pkg->AddField("U_guess", Metadata(flags_u_guess, std::vector<int>({nvar})));
    pkg->AddField("P_guess", Metadata(flags_p_guess, std::vector<int>({nvar})));

    Real u_rad_floor = pin->GetOrAddReal("radM1", "u_rad_floor", 1.e-8);
    pkg->AllParams().Add("u_rad_floor", u_rad_floor);


    // Real kappa_ep = pin->GetOrAddReal("radM1", "kappa_ep", 0.0);
    // Real sigma = pin->GetOrAddReal("radM1", "sigma", 0.0);
    // pkg->AllParams().Add("kappa_ep", kappa_ep);
    // pkg->AllParams().Add("sigma", sigma);

    //Right now, to execute the torus problem with radM1, we need to initialize the radiation primitives in fm_torus.cpp (this is stupid) (ASK Cora)
    //I think this should be moved to RadM1 method, maybe call an initialization method like a task straight after initializing the torus?
    //Especially because we'll want to initialize the radiation field for other problems. 

    //This method should allow you to add source terms to both plasma and radiation variables separately.
    pkg->AddSource = RadM1::AddImplicitRadiationSourceTerms;

    //Add inversion to the tasks
    pkg->BlockUtoP = RadM1::BlockUtoP;

    pkg->BlockApplyFloors = RadM1::ApplyRadM1Floors;


    return pkg;
}

void RadM1::ApplyRadM1Floors(MeshBlockData<Real> *rc, IndexDomain domain)
{
    auto pmb = rc->GetBlockPointer();
    auto& params = pmb->packages.Get("RadM1")->AllParams();

    const Real erad_floor   = params.Get<Real>("u_rad_floor");
    PackIndexMap prims_map;
    auto P = rc->PackVariables({Metadata::GetUserFlag("Primitive")}, prims_map);
    const VarMap m_p(prims_map, false);

    // We need to check if we actually have B fields enabled to avoid segfaults
    const bool has_b_field = pmb->packages.AllPackages().count("B_FluxCT") || 
                             pmb->packages.AllPackages().count("B_CD");

    auto bounds = pmb->cellbounds;
    const IndexRange ib = bounds.GetBoundsI(domain);
    const IndexRange jb = bounds.GetBoundsJ(domain);
    const IndexRange kb = bounds.GetBoundsK(domain);

    pmb->par_for("ApplyRadM1Floors", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            
            // Absolute Floor
            Real ehat = P(m_p.UU_RAD, k, j, i);
            if (ehat < erad_floor) {
                ehat = erad_floor;
                P(m_p.UU_RAD, k, j, i) = erad_floor;
            }
        }
    );
}


TaskStatus RadM1::BlockPtoU(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;

    // Pack Conserved Variables (Destination)
    auto& U = rc->PackVariables(std::vector<std::string>{"cons.u_rad", "cons.uvec_rad"});

    // Pack Primitive Variables (Source)
    PackIndexMap prim_map;
    auto P = rc->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    // Get Loop Bounds
    IndexRange3 b = KDomain::GetRange(rc, domain, coarse);

    // Parallel Loop
    pmb->par_for("RadM1_PtoU", b.ks, b.ke, b.js, b.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {

        Real Erf = P(m_p.UU_RAD, k, j, i);
        Real uvec_radframe[4] = {0, P(m_p.U1_RAD, k, j, i), P(m_p.U2_RAD, k, j, i), P(m_p.U3_RAD, k, j, i)};
        const Real gamma = GRMHD::lorentz_calc(G, uvec_radframe, k, j, i, Loci::center);
        Real ucon_rad[GR_DIM];
        //GRMHD::calc_ucon(G, uvec_radframe, k, j, i, Loci::center, ucon_rad);
        calc_ucon_rad(G, P, m_p, k, j, i, Loci::center, ucon_rad);
        // R^t^mu
        Real R_t_con[GR_DIM];

        for (int mu=0; mu<GR_DIM; ++mu) {
            R_t_con[mu] = 4./3. * Erf * ucon_rad[0] * ucon_rad[mu] + 1./3. * Erf * G.gcon(Loci::center, j, i, 0, mu);
        }

        // Calculat R^t_mu
        Real R_t_cov[GR_DIM];
        G.lower(R_t_con, R_t_cov, k, j, i, Loci::center);

        U(0, k, j, i) = R_t_cov[0] * G.gdet(Loci::center, j, i); // cons.u_rad
        U(1, k, j, i) = R_t_cov[1] * G.gdet(Loci::center, j, i); // cons.uvec_rad 1
        U(2, k, j, i) = R_t_cov[2] * G.gdet(Loci::center, j, i); // cons.uvec_rad 2
        U(3, k, j, i) = R_t_cov[3] * G.gdet(Loci::center, j, i); // cons.uvec_rad 3
    });
    return TaskStatus::complete;
}


KOKKOS_INLINE_FUNCTION void ApplyColdClosureFix(const GRCoordinates& G, 
                                                const Real R_t_cov_orig[GR_DIM], 
                                                const double gammarel2_fixed,
                                                const int& j, const int& i, 
                                                Real& new_R_t_t, 
                                                Real& Erf) {
    
    Real gcon_tt = G.gcon(Loci::center, j, i, 0, 0);
    
    // Time-Space cross term: g^{ti} R_i
    Real dot_t_i = G.gcon(Loci::center, j, i, 0, 1) * R_t_cov_orig[1] +
                   G.gcon(Loci::center, j, i, 0, 2) * R_t_cov_orig[2] +
                   G.gcon(Loci::center, j, i, 0, 3) * R_t_cov_orig[3];

    // Purely Spatial term: g^{ij} R_i R_j
    Real dot_spatial = G.gcon(Loci::center, j, i, 1, 1) * (R_t_cov_orig[1] * R_t_cov_orig[1]) +
                 2.0 * G.gcon(Loci::center, j, i, 1, 2) * (R_t_cov_orig[1] * R_t_cov_orig[2]) +
                       G.gcon(Loci::center, j, i, 2, 2) * (R_t_cov_orig[2] * R_t_cov_orig[2]) +
                 2.0 * G.gcon(Loci::center, j, i, 1, 3) * (R_t_cov_orig[1] * R_t_cov_orig[3]) +
                 2.0 * G.gcon(Loci::center, j, i, 2, 3) * (R_t_cov_orig[2] * R_t_cov_orig[3]) +
                       G.gcon(Loci::center, j, i, 3, 3) * (R_t_cov_orig[3] * R_t_cov_orig[3]);

    // In HARM/Kharma, alpha = 1.0 / sqrt(-gcon_tt). 
    // Therefore utsq = gammarel2 / alpha^2 simplifies to:
    Real utsq = -gammarel2_fixed * gcon_tt; 

    // The massive C-string expands exactly to (dot_t_i^2 - gcon_tt * dot_spatial).
    Real radical_inside = (dot_t_i * dot_t_i - gcon_tt * dot_spatial) 
                        * utsq * (gcon_tt + utsq) 
                        * m::pow(gcon_tt + 4.0 * utsq, 2);
                        
    Real radical = m::sqrt(m::max(0.0, radical_inside));

    // Calculate new covariant R_t (Avcov[0] in C)
    new_R_t_t = 0.25 * (-4.0 * dot_t_i * utsq * (gcon_tt + utsq) + radical) 
              / (gcon_tt * utsq * (gcon_tt + utsq));

    // Calculate fluid-frame energy Erf
    Erf = 0.75 * radical / (utsq * (gcon_tt + utsq) * (gcon_tt + 4.0 * utsq));
}

KOKKOS_INLINE_FUNCTION double CalculateGammaRel2(const GRCoordinates& G, const Real R_t_cov[GR_DIM], const Real invariant_scalar, const int& j, const int& i) {
    Real gcon_tt = G.gcon(Loci::center, j, i, 0, 0);
    Real R_t_t = R_t_cov[0];
    
    // Time-Space cross term: g^{ti} R^t_i
    Real dot_t_i = G.gcon(Loci::center, j, i, 0, 1) * R_t_cov[1] +
                   G.gcon(Loci::center, j, i, 0, 2) * R_t_cov[2] +
                   G.gcon(Loci::center, j, i, 0, 3) * R_t_cov[3];

    // Purely Spatial term: g^{ij} R^t_i R^t_j
    Real dot_spatial = G.gcon(Loci::center, j, i, 1, 1) * (R_t_cov[1] * R_t_cov[1]) +
                   2.0 * G.gcon(Loci::center, j, i, 1, 2) * (R_t_cov[1] * R_t_cov[2]) +
                         G.gcon(Loci::center, j, i, 2, 2) * (R_t_cov[2] * R_t_cov[2]) +
                   2.0 * G.gcon(Loci::center, j, i, 1, 3) * (R_t_cov[1] * R_t_cov[3]) +
                   2.0 * G.gcon(Loci::center, j, i, 2, 3) * (R_t_cov[2] * R_t_cov[3]) +
                         G.gcon(Loci::center, j, i, 3, 3) * (R_t_cov[3] * R_t_cov[3]);

    // Calculate Roots for (u^t_R)^2
    Real radical_inside = 4.0 * (gcon_tt * gcon_tt) * (R_t_t * R_t_t)
                        + (dot_t_i * dot_t_i)
                        + gcon_tt * (8.0 * R_t_t * dot_t_i + 3.0 * dot_spatial);
    Real radical = m::sqrt(m::max(0.0, radical_inside)); 

    Real num_a = 2.0 * (gcon_tt * gcon_tt) * (R_t_t * R_t_t)
                + dot_t_i * (dot_t_i + radical)
                + gcon_tt * (4.0 * R_t_t * dot_t_i + dot_spatial + R_t_t * radical);

    Real num_b = -2.0 * (gcon_tt * gcon_tt) * (R_t_t * R_t_t)
                - gcon_tt * (4.0 * R_t_t * dot_t_i + dot_spatial)
                + gcon_tt * R_t_t * radical
                + dot_t_i * (-dot_t_i + radical);

    Real gamma2a = -0.25 * num_a / invariant_scalar;
    Real gamma2b =  0.25 * num_b / invariant_scalar;

    // Sadowski 2013: Select the largest real root for (u^t)^2
    Real gamma2 = gamma2a;
    if (m::isnan(gamma2a) || m::isinf(gamma2a)) {
        gamma2 = gamma2b;
    } else if (!m::isnan(gamma2b) && !m::isinf(gamma2b)) {
        gamma2 = m::max(gamma2a, gamma2b);
    }

    Real alpha_sq = -1.0 / gcon_tt;
    Real gammarel2 = gamma2 * alpha_sq; 

    // Hard floor for physical bounds (Lorentz factor squared MUST be >= 1.0)
    // TODO: FLOOR! CHANGE THIS
    Real GAMMA_SMALL_LIMIT = (1.0-1e-10);
    if (gammarel2 < 1.0 && gammarel2 > GAMMA_SMALL_LIMIT) {
        gammarel2 = 1.0;
    }
    
    return gammarel2;
}


TaskStatus RadM1::BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;

    //Pack Variables
    auto U = rc->PackVariables(std::vector<std::string>{"cons.u_rad", "cons.uvec_rad"});

    // Pack Primitive Variables
    PackIndexMap prim_map;
    auto P = rc->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    //Get Loop Bounds
    IndexRange3 b = KDomain::GetRange(rc, domain, coarse);

    auto& params = pmb->packages.Get("RadM1")->AllParams();
    //TODO: FLOOR! CHANGE THIS
    const Real min_erad   = 10.*1.e-80;

    //Parallel Loop
    pmb->par_for("RadM1_UtoP", b.ks, b.ke, b.js, b.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {


            
            Real gdet = G.gdet(Loci::center, j, i);

            // Reconstruct Radiation Stress-Energy Tensor R^{mu, nu} 
            // For the M1 scheme, we assume the radiation is isotropic and satisfies the Eddington approximation (P^{ij} = (1/3) E delta^{ij} in the fluid frame)
            // but in the radiation frame. Therefore, \bar{R}^{tt} = E, \bar{R}^{ii} = \bar{E}/3, and every other component is zero. 
            // So we have Equation 27 in Sadowski et al. 2013 in the radiation rest frame, but since it's covariant, it's valid for every other frame (including lab frame).
            // The equation goes as follows:
            // R^{mu nu} = (4/3) \bar{E} u_R^mu u_R^nu + (1/3) \bar{E} g^{mu nu}, \bar{E} is always in the radiation rest frame.


            //Recover R^t_mu
            // U(0) is R^t_t, U(1..3) are R^t_i
            Real R_t_cov[GR_DIM] = {
                U(0, k, j, i) / gdet,
                U(1, k, j, i) / gdet,
                U(2, k, j, i) / gdet,
                U(3, k, j, i) / gdet
            };

            // Get R^{t\mu} from R^t_\mu by doing R^{t\mu} = g^{t\nu} R_{\nu\mu}
            // This is R^t_\mu
            Real R_t_con[GR_DIM];
            G.raise(R_t_cov, R_t_con, k, j, i, Loci::center);

            //Calculate Invariant Scalar S = R^t_mu * R^{t\mu}
            // Real invariant_scalar = 0.0;
            // for(int mu=0; mu<4; ++mu) {
            //      // Note: R_t_cov is R^t_mu (mixed), R_t_con is R^{t\mu} (upper). 
            //      // S = g_{mu nu} R^{t mu} R^{t nu} 
            //      for(int nu=0; nu<4; ++nu) {
            //         invariant_scalar += G.gcov(Loci::center, j, i, mu, nu) * R_t_con[mu] * R_t_con[nu];
            //      }
            // }
            Real invariant_scalar = 0.0;
            for(int mu=0; mu<4; ++mu) {
                invariant_scalar += R_t_cov[mu] * R_t_con[mu];
            }

            // From Sadowski et al. Eq 32 and 33
            // Solving by isolating Isolaring u^t_R^2 results in a much better equation, but I think (and I could be wrong)
            // that solving by isolating \bar{E} results in a much more numerical stable solution.

            // Calculating u^t_R^2
            Real gammarel2 = CalculateGammaRel2(G, R_t_cov, invariant_scalar, j, i);

            // Pre-calculate alpha_sq so E_rf can use it
            Real alpha_sq = -1.0 / G.gcon(Loci::center, j, i, 0, 0); 
            Real alpha = m::sqrt(alpha_sq);
            Real E_rf = (3.0 * R_t_con[0] * alpha_sq) / (4.0 * gammarel2 - 1.0);
            
            // TODO: FLOOR! CHANGE THIS
            Real GAMMAMAX = 20.;
            int nonfailure = gammarel2 >= 1.0 && E_rf > min_erad && gammarel2 <= (GAMMAMAX * GAMMAMAX)/(min_erad * min_erad);
            Real uvec_radframe_con[4] = {0};

            if(nonfailure){
                for(int mu=0; mu<4; ++mu) {
                    uvec_radframe_con[mu] = alpha * (R_t_con[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2 - 1.0)) / (4./3. * E_rf * m::sqrt(gammarel2));
                }
            } else {
                // Attempt Cold Closure
                Real gammarel2_slow = m::pow(1.0 +10.0 * std::numeric_limits<double>::epsilon(), 2.0);
                Real gammarel2_fast = GAMMAMAX * GAMMAMAX;

                Real R_t_t_slow, Erf_slow;
                ApplyColdClosureFix(G, R_t_cov, gammarel2_slow, j, i, R_t_t_slow, Erf_slow);

                Real R_t_t_fast, Erf_fast;
                ApplyColdClosureFix(G, R_t_cov, gammarel2_fast, j, i, R_t_t_fast, Erf_fast);

                Real R_t_t_new, gammarel2_new;
                
                if (m::abs(R_t_t_slow - R_t_cov[0]) > m::abs(R_t_t_fast - R_t_cov[0])) {
                    R_t_t_new = R_t_t_fast;
                    E_rf = Erf_fast;
                    gammarel2_new = gammarel2_fast;
                } else {
                    R_t_t_new = R_t_t_slow;
                    E_rf = Erf_slow;
                    gammarel2_new = gammarel2_slow;
                }

                Real R_t_cov_new[GR_DIM] = {R_t_t_new, R_t_cov[1], R_t_cov[2], R_t_cov[3]};

                Real R_t_con_new[GR_DIM];
                G.raise(R_t_cov_new, R_t_con_new, k, j, i, Loci::center);

                // Calculate primitive velocities with the new state
                if (E_rf > 0.0) {
                    for(int mu=0; mu<4; ++mu) {
                        uvec_radframe_con[mu] = alpha * (R_t_con_new[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2_new - 1.0)) / (4./3. * E_rf * m::sqrt(gammarel2_new));
                    }
                } else {
                    for(int mu=0; mu<4; ++mu) uvec_radframe_con[mu] = 0.0;
                }
            }

            // Pack Results
            P(m_p.UU_RAD, k, j, i) = E_rf;
            P(m_p.U1_RAD, k, j, i) = uvec_radframe_con[1];
            P(m_p.U2_RAD, k, j, i) = uvec_radframe_con[2];
            P(m_p.U3_RAD, k, j, i) = uvec_radframe_con[3];

    });
    return TaskStatus::complete;
}


KOKKOS_INLINE_FUNCTION void calc_Gnu(
    const GRCoordinates& G, const VariablePack<Real>& P_guess, const VarMap& m_p, 
    const int k, const int j, const int i, 
    const Real kappa_ep, const Real sigma, const Real gam, 
    GReal* Gnu_gdet)
{
    Gnu_gdet[0] = 0.0; // Energy exchange term
    Gnu_gdet[1] = 0.0; // Momentum exchange term in x
    Gnu_gdet[2] = 0.0; // Momentum exchange term in
    Gnu_gdet[3] = 0.0; // Momentum exchange term in z
    
}

KOKKOS_INLINE_FUNCTION void calc_Rtt_ff(const GRCoordinates& G, const VariablePack<Real>& P, const VarMap& m_p, const int& k, const int& j, const int& i, const Loci loc, Real * Rtt, Real ucon_gas[GR_DIM])
{
    Real Rij[GR_DIM][GR_DIM];

    Real ucov_gas[GR_DIM];
    G.lower(ucon_gas, ucov_gas, k, j, i, loc);

    for (int mu= 0; mu < GR_DIM; mu++){
        calc_tensor_m1(G, P, m_p, k, j, i, loc, Rij[mu]);
    }
    Real Rtt = 0.0;

    for (int mu=0; mu<GR_DIM; ++mu) {
        for (int nu=0; nu<GR_DIM; ++nu) {
            Rtt += - Rij[mu][nu] * ucon_gas[nu] * ucov_gas[mu]; 
        }
    }
}

KOKKOS_INLINE_FUNCTION int SolveImplicitLab4DPrimitives(
    const GRCoordinates& G, const const VariablePack<Real>& P_guess, const VarMap& m_p, 
    const int k, const int j, const int i, 
    const Real kappa_ep, const Real sigma, const Real gam, 
    Real Gnu_gdet[GR_DIM], Real Rtt_ff)
{
    
}

// This function will be a wrapper for the different implicit methods used to calculate the radiation four-force.
// For now I'm following Cora's advice to solely implement the secant method. It has no safe guards.
TaskStatus RadM1::RadiationImplicitRoutine(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
    const auto& G = pmb0->coords;

    const Real gam = pmb0->packages.Get("GRMHD")->Param<Real>("gamma");
    // Get from inverter package. Errors and tolerances
    auto &pars_inverter = pmb0->packages.Get("Inverter")->AllParams();

    //This is the threshold for determining if MHD or RAD dominate.
    const Real RADIMPLICITTHREASHOLD = 1e-2;

    // Pack Conserved Variables
    PackIndexMap cons_map;
    auto U = md->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    const VarMap m_u(cons_map, true); // Now map is safely populated

    PackIndexMap prim_map;
    auto P = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    const Real dt =  pmb0->packages.Get("Globals")->Param<Real>("dt_last"); // I have no idea if this is the way to get the last dt, but I assume it is. Copied straight out of hubble.cpp
    auto dU = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved});

    auto U_guess = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("RadGuessU")});
    auto P_guess = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("RadGuessP")});
    auto U_guess_init = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("RadGuessU")});
    auto P_guess_init = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("RadGuessP")});


    // Get Loop Bounds
    const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, dU.GetDim(5) - 1};

    pmb0->par_for("Implicit_Gnu_calculation", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            
            //get state after advection and geometry source terms for conserved variables.
            //The P variables will be a reasonable estimate of primitives.
            //Now I need to perform both PtoU convertion in gas and primitives
            BlockPtoU(pmb0->GetBlockData(b), IndexDomain::interior, false); //Radiation
            GRMHD::p_to_u(G, P, m_p, gam, k, j, i, U, m_u, Loci::center); //GRMHD

            P_guess = P;
            U_guess = U;
            U_guess_init = U;
            P_guess_init = P;


            //Now we calculate Ehat, which is E in the fluid frame.
            //Depending on it's value, we'll define if the cell is RAD or MHD dominated.
            Real ugas_con[GR_DIM];
            GRMHD::calc_ucon(G, ugas_con, k, j, i, Loci::center);
            Ehat = - calc_Rtt_ff(G, P, m_p, k, j, i, Loci::center, ugas_con); // This is R^{tt} in the fluid frame.

            // 0 means start by assuming MHD dominates and 1 means that assume RAD dominates.
            int start_with = 0;
            if(Ehat < RADIMPLICITTHREASHOLD * P_guess(m_p.UU, k, j, i)) start_with = 1;
            int nimplicit_iter = 0;
            int iter_max = 1;
            int success;
            for(nimplicit_iter = 0; nimplicit_iter < iter_max; ++nimplicit_iter){
                success = 0;

                //First we are gonna try the first case;
                if (success != 1){
                    P_guess = P_guess_init;
                    U_guess = U_guess_init;

                    success = SolveImplicitLab4DPrimitives(G, U_guess, P_guess, m_p, k, j, i, kappa_ep, sigma, gam, Gnu_gdet, Ehat, P);
                }

                //switch params MHD <--> RAD
                if (success != 1){
                    P_guess = P_guess_init;
                    U_guess = U_guess_init;

                    success = SolveImplicitLab4DPrimitives(G, U_guess, P_guess, m_p, k, j, i, kappa_ep, sigma, gam, Gnu_gdet, Ehat, P);
                }

                // use entropy fluid frame equations and switch MHD <->RAD
                if (sucess != 1){
                    P_guess = P_guess_init;
                    U_guess = U_guess_init;

                    success = SolveImplicitLab4DPrimitives(G, U_guess, P_guess, m_p, k, j, i, kappa_ep, sigma, gam, Gnu_gdet, Ehat, P);
                }
            }
        }
    );

    return TaskStatus::complete;
}

