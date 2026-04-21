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
    if (!units_enabled) {
        printf("\033[1;31mError: Units package not enabled! It must be enabled with/BEFORE RadM1.\033[0m\n");
        exit(1);
    }
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
    // I believe Metadata::Cell means these variables are defined at cell centers, but I should ask Ben.
    // We also add the "areWeImplicit" flag, which is either "Implicit" or "Explicit" based on the user's choice in the input file.
    std::vector<MetadataFlag> flags_radm1 = {Metadata::Cell, areWeImplicit, Metadata::GetUserFlag("RADM1")};


    //Retrieves the existing flags for the primitive and conserved variables, and adds the new radM1 flags to them.
    //Then adds the new variables for the radiation primitives and conserved variables with these flags.
    auto flags_prim = driver.Get<std::vector<MetadataFlag>>("prim_flags");
    flags_prim.insert(flags_prim.end(), flags_radm1.begin(), flags_radm1.end());
    auto flags_cons = driver.Get<std::vector<MetadataFlag>>("cons_flags");
    flags_cons.insert(flags_cons.end(), flags_radm1.begin(), flags_radm1.end());

    // I think this will push this metadata to restart files
    flags_prim.push_back(Metadata::Restart);

    // sync variables across boundaries (ASK BEN)
    if (pin->GetOrAddBoolean("RadM1", "sync_utop_seed", true)) { 
        flags_prim.push_back(Metadata::FillGhost);
    }
    //

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

    Real u_rad_floor = pin->GetOrAddReal("radM1", "u_rad_floor", 1.e-50);
    Real u_rad_max_floor = pin->GetOrAddReal("radM1", "u_rad_max_floor", 1.e20);

    pkg->AllParams().Add("u_rad_floor", u_rad_floor);
    pkg->AllParams().Add("u_rad_max_floor", u_rad_max_floor);
   

    // New Ratio Parameters (matching my legacy code)
    // Default values are placeholders; you should set reasonable defaults or require them in input.
    pkg->AllParams().Add("rad_rho_min", pin->GetOrAddReal("radM1", "rad_rho_min", 1.e-20));
    pkg->AllParams().Add("rad_rho_max", pin->GetOrAddReal("radM1", "rad_rho_max", 1.e20));
    pkg->AllParams().Add("rad_u_min", pin->GetOrAddReal("radM1", "rad_u_min", 1.e-20));
    pkg->AllParams().Add("rad_u_max", pin->GetOrAddReal("radM1", "rad_u_max", 1.e20));

    // Print all the parameters
    if (MPIRank0()){
        printf("RadM1 floor Parameters:\n");
        printf("u_rad_floor: %e\n", u_rad_floor);
        printf("u_rad_max_floor: %e\n", u_rad_max_floor);
        printf("rad_rho_min: %e\n", params.Get<Real>("rad_rho_min"));
        printf("rad_rho_max: %e\n", params.Get<Real>("rad_rho_max"));
        printf("rad_u_min: %e\n", params.Get<Real>("rad_u_min"));
        printf("rad_u_max: %e\n", params.Get<Real>("rad_u_max"));
    }

    // Magnetic ratio (if needed)
    pkg->AllParams().Add("rad_b_max", pin->GetOrAddReal("radM1", "rad_b_max", 100.0));
    
    //Right now, to execute the torus problem with radM1, we need to initialize the radiation primitives in fm_torus.cpp (this is stupid) (ASK BEN)
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
    const Real erad_max_floor = params.Get<Real>("u_rad_max_floor");
    const Real erad_rho_min = params.Get<Real>("rad_rho_min");
    const Real erad_rho_max = params.Get<Real>("rad_rho_max");
    const Real erad_u_min   = params.Get<Real>("rad_u_min");
    const Real erad_u_max   = params.Get<Real>("rad_u_max");
    const Real erad_b_max   = params.Get<Real>("rad_b_max");

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

            // //Radiation vs Density
            // Real rho = P(m_p.RHO, k, j, i);
            
            // // Radiation too small compared to mass
            // if (ehat < erad_rho_min * rho) {
            //     // Boost Radiation
            //     ehat = erad_rho_min * rho;
            //     P(m_p.UU_RAD, k, j, i) = ehat;
            // }

            // // Radiation dominates mass too much
            // if (ehat > erad_rho_max * rho) {
            //     // Boost Density -> modifying fluid var from Rad package? I'm blindly following Korals floors checks, should talk to ben
            //     // P(m_p.RHO, k, j, i) = ehat / erad_rho_max;
            //     P(m_p.UU_RAD, k, j, i) = erad_rho_max * rho;
            // }

            // //Radiation vs Internal Energy
            // Real u_gas = P(m_p.UU, k, j, i);

            // if (ehat < erad_u_min * u_gas) {
            //     // Boost Radiation
            //     ehat = erad_u_min * u_gas;
            //     P(m_p.UU_RAD, k, j, i) = ehat;
            // }

            // // if (ehat > erad_u_max * u_gas) {
            // //     // Boost internal energy -> modifying fluid var from Rad package? I'm blindly following Korals floors checks, should talk to ben
            // //     // P(m_p.UU, k, j, i) = ehat / erad_u_max;
            // //     P(m_p.UU_RAD, k, j, i) = erad_u_max * u_gas;
            // // }

            // //Radiation and Magnetic Pressure
            // if (has_b_field) {
            //     FourVectors Dtmp;
            //     GRMHD::calc_4vecs(G, P, m_p, k, j, i, Loci::center, Dtmp);

            //     GReal bsq = 0;
            //     DLOOP2 bsq += G.gcov(Loci::center, j, i, mu, nu) * Dtmp.bcon[mu] * Dtmp.bcon[nu];
            //     GReal mag_pressure = 0.5 * bsq;

            //     //Apply Magnetic Floor
            //     // Simplified: Boost radiation to match magnetic floor
            //     if (mag_pressure > erad_b_max * ehat) {
            //         P(m_p.UU_RAD, k, j, i) = mag_pressure / erad_b_max;
            //     }
            // }
        }
    );
}

TaskStatus RadM1::BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;

    //Pack Variables
    auto& U = rc->PackVariables(std::vector<std::string>{"cons.u_rad", "cons.uvec_rad"});

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

            //print if any nans in R_t_con, in any of the arrays
            // if (1) {
            //     for(int mu=0; mu<4; ++mu) {
            //         if (m::isnan(R_t_con[mu]) || m::isinf(R_t_con[mu])) {
            //             printf("NaN or Inf detected in R_t_con at block %d, k=%d, j=%d, i=%d, component %d\n", pmb->gid, k, j, i, mu);
            //                 return TaskStatus::complete;
            //         }
            //     }
            // }

            //Calculate Invariant Scalar S = R^t_mu * R^{t\mu}
            Real invariant_scalar = 0.0;
            for(int mu=0; mu<4; ++mu) {
                 // Note: R_t_cov is R^t_mu (mixed), R_t_con is R^{t\mu} (upper). 
                 // S = g_{mu nu} R^{t mu} R^{t nu} 
                 for(int nu=0; nu<4; ++nu) {
                    invariant_scalar += G.gcov(Loci::center, j, i, mu, nu) * R_t_con[mu] * R_t_con[nu];
                 }
            }

            
            // Isolaring u^t_R^2 in Equation 33 to find u^t_R in Equation 32 from Sadowski et al. 2013.
            // It gives g^{tt}\bar{E}^2 - 2 R^{tt} \bar{E} - 3 invariant_scalar = 0
            // It yields the solution \bar{E} = (R^{tt} +- sqrt((R^{tt})^2 + 3 g^tt invariant_scalar) )/ g^{tt}
            // We are gonna take the negative root since g^{tt} is negative and we want \bar{E} to be positive.
            Real g_con_tt = G.gcon(Loci::center, j, i, 0, 0);
            
            Real discr = R_t_con[0]*R_t_con[0] + 3.0 * g_con_tt * invariant_scalar;
            GReal E_bar = 0.0;
            if (discr < 0.0) {
                // Unphysical state (Flux > Energy), likely due to numerical error
                // Reset to valid small number
                 E_bar = 1e-300;
            } else {
                // Solve quadratic: E = (R^{tt} - sqrt((R^{tt})^2 + 3 g^{tt} S)) / g^{tt}
                E_bar = (R_t_con[0] - m::sqrt(discr)) / g_con_tt;
            }

            if (E_bar <= 0.0 || !m::isfinite(E_bar)) {
                E_bar = 1e-300;
            }

            // Recover u^t_R and u^i_R
            Real u_R_t = 0.0;
            //then u^t_R = sqrt(1/8 g^{tt} - 9/(8 E_bar^2) * invariant_scalar)
            Real val_ut = 0.125 * g_con_tt - 1.125/(E_bar * E_bar) * invariant_scalar;
            
            if (val_ut < 0.0) {
                // If velocity solution is invalid, assume static in lab frame or fluid frame
                // Fallback: simple approximation or fluid velocity
                u_R_t = m::sqrt(-1.0 / g_con_tt); // u^i = 0
            } else {
                u_R_t = m::sqrt(val_ut);
            }

            // Calculate spatial components u^i_R
            // Eq: u^\mu = (R^{t\mu} - E/3 g^{t\mu}) / (4/3 E u^t)
            Real u_R_con[GR_DIM];
            // Safety for division
            Real denom = (4.0/3.0) * E_bar * u_R_t;
            //if (m::abs(denom) < 1e-25) denom = 1e-25;

            u_R_con[0] = u_R_t; // We already found component 0
            for(int nu=1; nu<4; ++nu) {
                u_R_con[nu] = (R_t_con[nu] - (E_bar / 3.0) * G.gcon(Loci::center, j, i, 0, nu)) / denom;
            }

            // Now the entire R^{\mu\nu} can be reconstructed from E_bar and u_R_con using Equation 27 in Sadowski et al. 2013

            Real R_con[GR_DIM][GR_DIM];
            for(int mu=0; mu<4; ++mu) {
                for(int nu=0; nu<4; ++nu) {
                    R_con[mu][nu] = (4.0/3.0) * E_bar * u_R_con[mu] * u_R_con[nu] + (E_bar / 3.0) * G.gcon(Loci::center, j, i, mu, nu);
                }
            }
 
            // We need to project R^munu onto the fluid 4-velocity u_gas
            // to get the fluid-frame energy E_hat and flux F_hat.

            //Get the Fluid 4-Velocity (u_gas) from Primitives
            Real uvec_gas[3] = {
                P(m_p.U1, k, j, i), 
                P(m_p.U2, k, j, i), 
                P(m_p.U3, k, j, i)
            };
            
            Real ucon_gas[GR_DIM], ucov_gas[GR_DIM];
            // Calculate u^mu (contravariant) and u_mu (covariant) for the gas
            GRMHD::calc_ucon(G, uvec_gas, k, j, i, Loci::center, ucon_gas);
            G.lower(ucon_gas, ucov_gas, k, j, i, Loci::center);

            //Calculate Fluid-Frame Energy Density (E_hat)
            // Projection: E_hat = R^{ab} * u_a * u_b
            Real E_hat = 0.0;
            for(int mu=0; mu<4; ++mu) {
                for(int nu=0; nu<4; ++nu) {
                    E_hat += R_con[mu][nu] * ucov_gas[mu] * ucov_gas[nu];
                }
            }


            // Calculate Fluid-Frame Flux (F_hat)
            // First, compute the Momentum Density J^mu = - R^{mu nu} u_nu
            Real J_con[GR_DIM] = {0.0};
            for(int mu=0; mu<4; ++mu) {
                for(int nu=0; nu<4; ++nu) {
                    J_con[mu] -= R_con[mu][nu] * ucov_gas[nu];
                }
            }

            // Project J^mu orthogonal to u_gas to get Flux F^mu
            // Formula: F^mu = (g^{mu nu} + u^mu u^nu) J_nu
            // Simplifies to: F^mu = J^mu - (J dot u) * u^mu
            // Note: (J dot u) is simply -E_hat
            
            Real F_con[GR_DIM]; 
            for(int mu=0; mu<4; ++mu) {
                F_con[mu] = J_con[mu] - (E_hat * ucon_gas[mu]);
            }

            if(isnan(E_hat) || isinf(E_hat)) {
                P(m_p.UU_RAD, k, j, i) =  0.001 * P(m_p.UU, k, j, i);
                P(m_p.U1_RAD, k, j, i) = 0.0;
                P(m_p.U2_RAD, k, j, i) = 0.0;
                P(m_p.U3_RAD, k, j, i) = 0.0; 
            }else{
                P(m_p.UU_RAD, k, j, i) = m::max(E_hat, 0.001 * P(m_p.UU, k, j, i)); // Convert back to conserved form and ensure positivity
                P(m_p.U1_RAD, k, j, i) = F_con[1];
                P(m_p.U2_RAD, k, j, i) = F_con[2];
                P(m_p.U3_RAD, k, j, i) = F_con[3];
            }

        }
    );
    return TaskStatus::complete;
}


KOKKOS_INLINE_FUNCTION void calc_Gnu(const Real rho, const Real ug, GReal * Gnu_gdet){

    Gnu_gdet[0] = 0.0; // Energy exchange term
    Gnu_gdet[1] = 0.0; // Momentum exchange term in x
    Gnu_gdet[2] = 0.0; // Momentum exchange term in
    Gnu_gdet[3] = 0.0; // Momentum exchange term in z
    
}


KOKKOS_INLINE_FUNCTION int SetImplicitInitialGuess(
    const GRCoordinates& G, const VarMap& m_p, const VarMap& m_u, 
    VariablePack<Real>& U, VariablePack<Real>& U_guess, VariablePack<Real>& P_guess, VariablePack<Real>& dU, 
    const Real dt, const int k, const int j, const int i, 
    GReal* Gnu_gdet, const int iter_max, const Real err_tol, const Real gam)
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
    calc_Gnu(0.0, 0.0, Gnu_gdet);
    
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
    GReal* Gnu_gdet, const int iter_max, const Real err_tol, const Real gam)
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
        
        calc_Gnu(0.0, 0.0, Gnu_gdet); // Recalculate Gnu with the new guess. Here we are just using a placeholder since we don't have the actual physics implemented yet.

        //Compute New Residuals
        for(int d = 0; d < 4; ++d) {
            Residual_curr[d] = R_curr[d] - U(rad_indices[d], k, j, i) + dt * Gnu_gdet[d];
            if (fabs(Residual_curr[d]) > max_res) max_res = fabs(Residual_curr[d]);
        }
        
        iter++;
    }
}

// This function will be a wrapper for the different implicit methods used to calculate the radiation four-force.
// For now I'm following Ben's advice to solely implement the secant method. It has no safe guards.
TaskStatus RadM1::AddImplicitRadiationSourceTerms(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{

    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
    const auto& G = pmb0->coords;

    const Real gam = pmb0->packages.Get("GRMHD")->Param<Real>("gamma");
    // Get from inverter package. Errors and tolerances
    auto &pars_inverter = pmb0->packages.Get("Inverter")->AllParams();
    const int iter_max = pmb0->packages.Get("Inverter")->Param<int>("iter_max");
    const Real err_tol = pmb0->packages.Get("Inverter")->Param<Real>("err_tol");

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
                G, m_p, m_u, U(b), U_guess(b), P_guess(b), dU(b), dt, k, j, i, local_Gnu_gdet, iter_max, err_tol, gam
            );

            if(proceed_with_implicit) {
                // Refine Gnu using the implicit solver.
                ImplicitSolver(G, m_p, m_u, P_guess(b), U_guess(b), U(b), dU(b), dt, k, j, i, local_Gnu_gdet, iter_max, err_tol, gam);
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

