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
#include "inverter.hpp"  // Add this include

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
        GRMHD::calc_ucon(G, uvec_radframe, k, j, i, Loci::center, ucon_rad);

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
    if (gammarel2 < 1.0 || m::isnan(gammarel2)) {
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
            Real E_rf = (3.0 * R_t_con[0] * alpha_sq) / (4.0 * gammarel2 - 1.0);
            //There were a lot of checks in this part, so I should test this very carefully.
            if(E_rf < 0.0 || m::isnan(E_rf) || m::isinf(E_rf)) {
                printf("Something Happened!, E_rf = %e, R_t_con[0]: %e, alpha_sq %e, gammarel2 %e\n", E_rf, R_t_con[0], alpha_sq, gammarel2);
            }
            Real uvec_radframe_con[4];
            Real alpha = m::sqrt(alpha_sq);
            for(int mu=0; mu<4; ++mu) {
                uvec_radframe_con[mu] = alpha * (R_t_con[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2 - 1.0))/(4./3. * E_rf * sqrt(gammarel2));
            }

            P(m_p.UU_RAD, k, j, i) = E_rf;
            P(m_p.U1_RAD, k, j, i) = uvec_radframe_con[1];
            P(m_p.U2_RAD, k, j, i) = uvec_radframe_con[2];
            P(m_p.U3_RAD, k, j, i) = uvec_radframe_con[3];

    });
    return TaskStatus::complete;
}


// KOKKOS_INLINE_FUNCTION void calc_Gnu(
//     const GRCoordinates& G, const VariablePack<Real>& P_guess, const VarMap& m_p, 
//     const int k, const int j, const int i, 
//     const Real kappa_ep, const Real sigma, const Real gam, 
//     GReal* Gnu_gdet) {

//     Gnu_gdet[0] = 0.0; // Energy exchange term
//     Gnu_gdet[1] = 0.0; // Momentum exchange term in x
//     Gnu_gdet[2] = 0.0; // Momentum exchange term in
//     Gnu_gdet[3] = 0.0; // Momentum exchange term in z
    
// }

KOKKOS_INLINE_FUNCTION void calc_Gnu(
    const GRCoordinates& G, const VariablePack<Real>& P_guess, const VarMap& m_p, 
    const int k, const int j, const int i, 
    const Real kappa_ep, const Real sigma, const Real gam, 
    GReal* Gnu_gdet) 
{
    Real gdet = G.gdet(Loci::center, j, i);
    
    // Get Gas primitives & Compute Temperature (T = P_gas / rho)
    Real rho = P_guess(m_p.RHO, k, j, i);
    Real uu_gas = P_guess(m_p.UU, k, j, i);
    Real P_gas = (gam - 1.0) * uu_gas; 
    Real T = P_gas / rho;
    
    //Get 4-velocity and FourVectors for the M1 Tensor
    FourVectors D_gas;
    GRMHD::calc_4vecs(G, P_guess, m_p, k, j, i, Loci::center, D_gas);
    
    // Construct the full lab-frame Radiation Tensor R^{\mu\nu}
    // calc_tensor_m1 computes R^{\mu\nu} for a specific direction \mu. 
    // We call it 4 times to fill the whole 4x4 matrix.
    Real R_tensor[4][4];
    RadM1::calc_tensor_m1(G, P_guess, m_p, 0, k, j, i, Loci::center, R_tensor[0]);
    RadM1::calc_tensor_m1(G, P_guess, m_p, 1, k, j, i, Loci::center, R_tensor[1]);
    RadM1::calc_tensor_m1(G, P_guess, m_p, 2, k, j, i, Loci::center, R_tensor[2]);
    RadM1::calc_tensor_m1(G, P_guess, m_p, 3, k, j, i, Loci::center, R_tensor[3]);
    
    // Compute G^\mu = - \rho \kappa ( R^{\mu\alpha} u_\alpha + \sigma T^4 u^\mu )
    Real kappa = kappa_ep * rho;
    Real emission_term = sigma * pow(T, 4);
    
    for (int mu = 0; mu < 4; ++mu) {
        Real R_u = 0.0;
        for (int alpha = 0; alpha < 4; ++alpha) {
            R_u += R_tensor[mu][alpha] * D_gas.ucov[alpha];
        }
        
        Real Gnu = -rho * kappa * (R_u + emission_term * D_gas.ucon[mu]);
        Gnu_gdet[mu] = Gnu * gdet;
    }
}


KOKKOS_INLINE_FUNCTION int SetImplicitInitialGuess(
    const GRCoordinates& G, const VarMap& m_p, const VarMap& m_u, 
    VariablePack<Real>& U, VariablePack<Real>& U_guess, VariablePack<Real>& P_guess, VariablePack<Real>& dU, 
    const Real dt, const int k, const int j, const int i, 
    GReal* Gnu_gdet, const int iter_max, const Real err_tol, const Real gam, const Real kappa_ep, const Real sigma)
{
// The function will follow:
// 1. Calculate U->P for gas and radiation variables
// 2. Then Recalculate U from P for gas variables.
// 3. Calculate the source term Gnu from this P.
// 4. Calculate \Delta U = (sqrt(-g)) * Gnu * dt
// 5. Calculate relative error for the gas variables:
    // Energy: 
        // Norm_0 = |U_n_0| + |U_{n+1}_0| + |dt * dU_0|
        // Error_0 = 0.25 * |U_{n+1}_0 - U_n_0 - dt * dU_0| / Norm_0

    // Momentum:
        // Combine the scale of all 3 momentum directions into one norm
        // Norm_mom = (|sqrt(g^11) * U_n_1| + |sqrt(g^11) * U_{n+1}_1| + |sqrt(g^11) * dt * dU_1|) +
        //            (|sqrt(g^22) * U_n_2| + |sqrt(g^22) * U_{n+1}_2| + |sqrt(g^22) * dt * dU_2|) +
        //            (|sqrt(g^33) * U_n_3| + |sqrt(g^33) * U_{n+1}_3| + |sqrt(g^33) * dt * dU_3|)
    
        // Fallback if the gas is perfectly stationary:
        // if (Norm_mom == 0) Norm_mom = Norm_0; 
    
        // Error_1 = 0.25 * sqrt(g^11) * |U_{n+1}_1 - U_n_1 - dt * dU_1| / Norm_mom
        // Error_2 = 0.25 * sqrt(g^22) * |U_{n+1}_2 - U_n_2 - dt * dU_2| / Norm_mom
        // Error_3 = 0.25 * sqrt(g^33) * |U_{n+1}_3 - U_n_3 - dt * dU_3| / Norm_mom
    
    // Total_Gas_Error = Error_0 + Error_1 + Error_2 + Error_3
// 6. Calculate the relative error for the radiation variables
    // Energy:
        // Norm_rad_0 = |U_n_RAD_0| + |U_{n+1}_RAD_0| + |dt * dU_RAD_0|
        // Error_rad_0 = 0.25 * |U_{n+1}_RAD_0 - U_n_RAD_0 + dt * dU_RAD_0| / Norm_rad_0
    // Momentum:
        // Combine the scale of all 3 momentum directions into one norm
        // Norm_rad_mom = (|sqrt(g^11) * U_n_RAD_1| + |sqrt(g^11) * U_{n+1}_RAD_1| + |sqrt(g^11) * dt * dU_RAD_1|) +
        //                (|sqrt(g^22) * U_n_RAD_2| + |sqrt(g^22) * U_{n+1}_RAD_2| + |sqrt(g^22) * dt * dU_RAD_2|) +
        //                (|sqrt(g^33) * U_n_RAD_3| + |sqrt(g^33) * U_{n+1}_RAD_3| + |sqrt(g^33) * dt * dU_RAD_3|)
        // Fallback if the radiation is perfectly stationary:
        // if (Norm_rad_mom == 0) Norm_rad_mom = Norm_rad_0;
        // Error_rad_1 = 0.25 * sqrt(g^11) * |U_{n+1}_RAD_1 - U_n_RAD_1 + dt * dU_RAD_1| / Norm_rad_mom
        // Error_rad_2 = 0.25 * sqrt(g^22) * |U_{n+1}_RAD_2 - U_n_RAD_2 + dt * dU_RAD_2| / Norm_rad_mom
        // Error_rad_3 = 0.25 * sqrt(g^33) * |U_{n+1}_RAD_3 - U_n_RAD_3 + dt * dU_RAD_3| / Norm_rad_mom
        // Total_Rad_Error = Error_rad_0 + Error_rad_1 + Error_rad_2 + Error_rad_3
    
    // Total error = Total_Gas_Error + Total_Rad_Error. If err < 1.e-12. We don't even need to do the implicit solve. This might be for cases where Gnu is
    // not stiff.

    // Note for future readers: If we evaluated the $\theta$ or $\phi$ direction individually when computing the error, 
    // the local momentum in that specific direction might be zero, causing a division by zero. By summing all three directions into Norm_mom,
    // we are essentially dividing by the "total momentum scale" of the fluid cell, which is much safer and more physically meaningful.


    // Here we are getting the first u_to_p conversion for the gas variables.
    int pflag = Inverter::u_to_p<Inverter::Type::kastaun>(G, U, m_u, gam, k, j, i, P_guess, m_p, Loci::center, iter_max, err_tol, false);

    //Now we have to recompute U_guess from P_guess for the gas variables.
    GRMHD::p_to_u(G, P_guess, m_p, gam, k, j, i, U_guess, m_u, Loci::center);
    
    // Calculate source terms from the initial guess of the primitives.
    calc_Gnu(G, P_guess, m_p, k, j, i, kappa_ep, sigma, gam, Gnu_gdet);
    
    // Compute the gas errors
    Real sqrt_g11 = sqrt(fabs(G.gcon(Loci::center, j, i, 1, 1)));
    Real sqrt_g22 = sqrt(fabs(G.gcon(Loci::center, j, i, 2, 2)));
    Real sqrt_g33 = sqrt(fabs(G.gcon(Loci::center, j, i, 3, 3)));

    Real norm_0 = fabs(U(m_u.UU, k, j, i)) + fabs(U_guess(m_u.UU, k, j, i)) + fabs(dt * dU(m_u.UU, k, j, i));
    Real err_0 = 0.0;
    if (norm_0 > 0.0) {
        err_0 = 0.25 * fabs(U_guess(m_u.UU, k, j, i) - U(m_u.UU, k, j, i) - dt * dU(m_u.UU, k, j, i)) / norm_0;
    }

    Real norm_mom = 
        (sqrt_g11 * fabs(U(m_u.U1, k, j, i)) + sqrt_g11 * fabs(U_guess(m_u.U1, k, j, i)) + sqrt_g11 * fabs(dt * dU(m_u.U1, k, j, i))) +
        (sqrt_g22 * fabs(U(m_u.U2, k, j, i)) + sqrt_g22 * fabs(U_guess(m_u.U2, k, j, i)) + sqrt_g22 * fabs(dt * dU(m_u.U2, k, j, i))) +
        (sqrt_g33 * fabs(U(m_u.U3, k, j, i)) + sqrt_g33 * fabs(U_guess(m_u.U3, k, j, i)) + sqrt_g33 * fabs(dt * dU(m_u.U3, k, j, i)));

    if (norm_mom == 0.0) norm_mom = norm_0;

    // Momentum Errors
    Real err_1 = 0.25 * sqrt_g11 * fabs(U_guess(m_u.U1, k, j, i) - U(m_u.U1, k, j, i) - dt * dU(m_u.U1, k, j, i)) / norm_mom;
    Real err_2 = 0.25 * sqrt_g22 * fabs(U_guess(m_u.U2, k, j, i) - U(m_u.U2, k, j, i) - dt * dU(m_u.U2, k, j, i)) / norm_mom;
    Real err_3 = 0.25 * sqrt_g33 * fabs(U_guess(m_u.U3, k, j, i) - U(m_u.U3, k, j, i) - dt * dU(m_u.U3, k, j, i)) / norm_mom;

    Real total_gas_error = err_0 + err_1 + err_2 + err_3;

    //Radiation errors
    Real total_rad_error = 0.0;

    //Final error
    Real total_error = total_gas_error + total_rad_error;
    

    if(total_error < 1e-12){
        // Skip implicit solve for this cell, for this we use flag = 0;
        return 0;
    }else{
        // Proceed with implicit solve, for this we use flag = 1;
        return 1;
    }
}

KOKKOS_INLINE_FUNCTION void ImplicitSolver(
    const GRCoordinates& G, const VarMap& m_p, const VarMap& m_u, 
    VariablePack<Real>& P_guess, VariablePack<Real>& U_guess, const VariablePack<Real>& U, const VariablePack<Real>& dU, 
    const Real dt, const int k, const int j, const int i, 
    GReal* Gnu_gdet, const int iter_max, const Real err_tol, const Real gam, const Real kappa_ep, const Real sigma)
{
    // This function will use the secant method to update the radiation variables.
    // There might be a better way to do this, I'm sure...
    
    Real R_old[4], R_curr[4], R_new[4];
    Real Residual_old[4], Residual_curr[4];
    
    int rad_indices[4] = {m_u.UU_RAD, m_u.U1_RAD, m_u.U2_RAD, m_u.U3_RAD};
    int gas_indices[4] = {m_u.UU, m_u.U1, m_u.U2, m_u.U3};
    
    for(int d = 0; d < 4; ++d) {
        // Set the initial guess for R to the current inverted value from the initial guess.
        R_curr[d] = U_guess(rad_indices[d], k, j, i);

        // For the secant method's first step, we need a slight perturbation for R_old
        R_old[d] = R_curr[d] * 0.999 + 1e-15; 
        
        // Residual = R_{guess} - R_n + dt * Gnu_guess
        // dR/dt = -Gnu, so R_np1 = R_n - dt * Gnu
        Residual_curr[d] = R_curr[d] - U(rad_indices[d], k, j, i) + dt * Gnu_gdet[d];
        Residual_old[d] = Residual_curr[d]; // Will be updated in loop
    }

    int iter = 0;
    int iter_max_secant = 10;
    int err_tol_secant = 1.e-12;
    Real max_res = 1.0;

    while (iter < iter_max_secant && max_res > err_tol_secant) {
        max_res = 0.0;
        
        for(int d = 0; d < 4; ++d) {
            Real dRes = Residual_curr[d] - Residual_old[d];
            
            // This could be zero, therefore I will redefine it as 1e-16, just cause it's the double precision limit, but we can change it if we want.
            if (fabs(dRes) < 1e-16) {
                R_new[d] = R_curr[d]; 
            } else {
                // Here we update the new try with the secant method formula
                // x_{new} = x - Res(x) * (x - x_{old}) / (Res(x) - Res(x_{old}))
                R_new[d] = R_curr[d] - Residual_curr[d] * (R_curr[d] - R_old[d]) / dRes;
            }
            
            // Update history
            R_old[d] = R_curr[d];
            Residual_old[d] = Residual_curr[d];
            R_curr[d] = R_new[d];
            
            // Update local state for next evaluation
            U_guess(rad_indices[d], k, j, i) = R_curr[d];
        }

        // Re-evaluate Physics with R_curr (New Guess)
        // Sadowski Eq 56 principle: U_gas_{new} is linked to R_{new} via conservation
        for(int d = 0; d < 4; ++d) {
            U_guess(gas_indices[d], k, j, i) = U(gas_indices[d], k, j, i) + dt * dU(gas_indices[d], k, j, i) - (R_curr[d] - U(rad_indices[d], k, j, i));
        }

        // Invert to Primitives and Re-calculate Gnu
    
        Inverter::u_to_p<Inverter::Type::kastaun>(G, U_guess, m_u, gam, k, j, i, P_guess, m_p, Loci::center, iter_max, err_tol, false);
        GRMHD::p_to_u(G, P_guess, m_p, gam, k, j, i, U_guess, m_u, Loci::center);
        
        calc_Gnu(G, P_guess, m_p, k, j, i, kappa_ep, sigma, gam, Gnu_gdet); // Recalculate Gnu with the new guess. Here we are just using a placeholder since we don't have the actual physics implemented yet.

        //Compute New Residuals
        for(int d = 0; d < 4; ++d) {
            Residual_curr[d] = R_curr[d] - U(rad_indices[d], k, j, i) + dt * Gnu_gdet[d];
            if (fabs(Residual_curr[d]) > max_res) max_res = fabs(Residual_curr[d]);
        }
        
        iter++;
    }
}

// This function will be a wrapper for the different implicit methods used to calculate the radiation four-force.
// For now I'm following Cora's advice to solely implement the secant method. It has no safe guards.
TaskStatus RadM1::AddImplicitRadiationSourceTerms(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
    const auto& G = pmb0->coords;

    const Real gam = pmb0->packages.Get("GRMHD")->Param<Real>("gamma");
    // Get from inverter package. Errors and tolerances
    auto &pars_inverter = pmb0->packages.Get("Inverter")->AllParams();
    const int iter_max = pmb0->packages.Get("Inverter")->Param<int>("iter_max");
    const Real err_tol = pmb0->packages.Get("Inverter")->Param<Real>("err_tol");

    // For the shock tube
    // const Real kappa_ep = pmb0->packages.Get("radM1")->Param<Real>("kappa_ep");
    // const Real sigma = pmb0->packages.Get("radM1")->Param<Real>("sigma");

    // get units from units package:
    const Real length_unit = pmb0->packages.Get("Units")->Param<Real>("length_unit_cgs");
    const Real rho_unit = pmb0->packages.Get("Units")->Param<Real>("rho_scale");
    const Real kappa_ep = 0.4 * rho_unit * length_unit; //kappa_ep has units of cm^2/g, 
    const Real sigma = 3.085e9 * length_unit; // sigma has units of cm^(-1)

    // Pack Conserved Variables
    PackIndexMap cons_map;
    auto U = md->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    const VarMap m_u(cons_map, true); // Now map is safely populated

    // Pack Primitive Variables
    PackIndexMap prim_map;
    auto P = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    const Real dt =  pmb0->packages.Get("Globals")->Param<Real>("dt_last"); // I have no idea if this is the way to get the last dt, but I assume it is. Copied straight out of hubble.cpp
    auto dU = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved});

    
    auto U_guess = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("RadGuessU")});
    auto P_guess = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("RadGuessP")});
    // Get Loop Bounds
    const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, dU.GetDim(5) - 1};

    pmb0->par_for("Implicit_Gnu_calculation", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            GReal local_Gnu_gdet[GR_DIM];

            int proceed_with_implicit = SetImplicitInitialGuess( 
                G, m_p, m_u, U(b), U_guess(b), P_guess(b), dU(b), dt, k, j, i, local_Gnu_gdet, iter_max, err_tol, gam, kappa_ep, sigma
            );

            if(proceed_with_implicit) {
                // Refine Gnu using the implicit solver.
                ImplicitSolver(G, m_p, m_u, P_guess(b), U_guess(b), U(b), dU(b), dt, k, j, i, local_Gnu_gdet, iter_max, err_tol, gam, kappa_ep, sigma);
            }

            dU(b, m_u.UU, k, j, i) += local_Gnu_gdet[0];
            dU(b, m_u.U1, k, j, i) += local_Gnu_gdet[1];
            dU(b, m_u.U2, k, j, i) += local_Gnu_gdet[2];
            dU(b, m_u.U3, k, j, i) += local_Gnu_gdet[3];

            dU(b, m_u.UU_RAD, k, j, i) -= local_Gnu_gdet[0];
            dU(b, m_u.U1_RAD, k, j, i) -= local_Gnu_gdet[1];
            dU(b, m_u.U2_RAD, k, j, i) -= local_Gnu_gdet[2];
            dU(b, m_u.U3_RAD, k, j, i) -= local_Gnu_gdet[3];
        }
    );

    return TaskStatus::complete;
}

