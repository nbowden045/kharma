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
    // (ASK BEN)
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


    Real u_rad_floor = pin->GetOrAddReal("RadM1", "u_rad_floor", 1.e-50);
    pkg->AllParams().Add("u_rad_floor", u_rad_floor);

    // New Ratio Parameters (matching your legacy code)
    // Default values are placeholders; you should set reasonable defaults or require them in input.
    pkg->AllParams().Add("rad_rho_min", pin->GetOrAddReal("radM1", "rad_rho_min", 1.e-20));
    pkg->AllParams().Add("rad_rho_max", pin->GetOrAddReal("radM1", "rad_rho_max", 1.e20));
    pkg->AllParams().Add("rad_u_min", pin->GetOrAddReal("radM1", "rad_u_min", 1.e-20));
    pkg->AllParams().Add("rad_u_max", pin->GetOrAddReal("radM1", "rad_u_max", 1.e20));

    // Print all the parameters
    if (MPIRank0()){
        printf("RadM1 floor Parameters:\n");
        printf("u_rad_floor: %e\n", u_rad_floor);
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
    pkg->AddSource = RadM1::AddSource;

    //Add inversion to the tasks
    pkg->BlockUtoP = RadM1::BlockUtoP;

    pkg->BlockApplyFloors = RadM1::ApplyRadM1Floors;


    return pkg;
}

void RadM1::ApplyRadM1Floors(MeshBlockData<Real> *rc, IndexDomain domain)
{
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;
    auto& params = pmb->packages.Get("RadM1")->AllParams();

    // 1. Retrieve all parameters
    const Real erad_floor   = params.Get<Real>("u_rad_floor");
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

            //Radiation vs Density
            Real rho = P(m_p.RHO, k, j, i);
            
            // Radiation too small compared to mass
            if (ehat < erad_rho_min * rho) {
                // Boost Radiation
                ehat = erad_rho_min * rho;
                P(m_p.UU_RAD, k, j, i) = ehat;
            }

            // Radiation dominates mass too much
            if (ehat > erad_rho_max * rho) {
                // Boost Density -> modifying fluid var from Rad package? I'm blindly following Korals floors checks, should talk to ben
                // P(m_p.RHO, k, j, i) = ehat / erad_rho_max;
                P(m_p.UU_RAD, k, j, i) = erad_rho_max * rho;
            }

            //Radiation vs Internal Energy
            Real u_gas = P(m_p.UU, k, j, i);

            if (ehat < erad_u_min * u_gas) {
                // Boost Radiation
                ehat = erad_u_min * u_gas;
                P(m_p.UU_RAD, k, j, i) = ehat;
            }

            // if (ehat > erad_u_max * u_gas) {
            //     // Boost internal energy -> modifying fluid var from Rad package? I'm blindly following Korals floors checks, should talk to ben
            //     // P(m_p.UU, k, j, i) = ehat / erad_u_max;
            //     P(m_p.UU_RAD, k, j, i) = erad_u_max * u_gas;
            // }

            //Radiation and Magnetic Pressure
            if (has_b_field) {
                FourVectors Dtmp;
                GRMHD::calc_4vecs(G, P, m_p, k, j, i, Loci::center, Dtmp);

                GReal bsq = 0;
                DLOOP2 bsq += G.gcov(Loci::center, j, i, mu, nu) * Dtmp.bcon[mu] * Dtmp.bcon[nu];
                GReal mag_pressure = 0.5 * bsq;

                //Apply Magnetic Floor
                // Simplified: Boost radiation to match magnetic floor
                if (mag_pressure > erad_b_max * ehat) {
                    P(m_p.UU_RAD, k, j, i) = mag_pressure / erad_b_max;
                }
            }
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


// Is this function called before or after the conversion of U to P? 
// I think it should be before, since we need the radiation four-force to update the conserved variables.
TaskStatus RadM1::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{
    auto pmb0 = mdudt->GetBlockData(0)->GetBlockPointer();

    // Pack Conserved Variables
    PackIndexMap cons_map;
    auto dUdt = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    auto U = md->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    const VarMap m_u(cons_map, true);

    // Pack Primitive Variables
    PackIndexMap prim_map;
    auto P = md->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    // Get Loop Bounds
    const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, dUdt.GetDim(5) - 1};

    // KOKKOS loop over the meshblock, applying the radiation four-force source term to both the fluid and radiation conserved variables.
    pmb0->par_for("add_radM1_source", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            
            const auto& G = dUdt.GetCoords(b);
            Real gdet = G.gdet(Loci::center, j, i); 
            
            // Read necessary primitive variables for calculating absorption terms

            // We are gonna follow Mckinney et al. 2014 paper Eq. 16
            // For this we need R^\mu_\nu, u^\nu, and the opacities. The conserved radiation variables are R^t_t and R^t_i
            // To reconstruct the full R^\mu_\nu tensor, we need to reconstruct the primitives (E,F^i) and then use the M1 closure to get R^\mu_\nu by using

    
            // Reconstruct Radiation Stress-Energy Tensor R^{mu, nu} 
            // For the M1 scheme, we assume the radiation is isotropic and satisfies the Eddington approximation (P^{ij} = (1/3) E delta^{ij} in the fluid frame)
            // but in the radiation frame. Therefore, \bar{R}^{tt} = E, \bar{R}^{ii} = \bar{E}/3, and every other component is zero. 
            // So we have Equation 27 in Sadowski et al. 2013 in the radiation rest frame, but since it's covariant, it's valid for every other frame (including lab frame).
            // The equation goes as follows:
            // R^{mu nu} = (4/3) \bar{E} u_R^mu u_R^nu + (1/3) \bar{E} g^{mu nu}, \bar{E} is always in the radiation rest frame.

            // Get R^{t\mu} from R^t_\mu by doing R^{t\mu} = g^{t\nu} R_{\nu\mu} (using the raise function from Kharma)
            // (ASK BEN), are these terms multiplying sqrt(-g) and is this the right way to access them?
            // This is R^t_\mu
            Real R_t_cov[GR_DIM] ={ U(b, m_u.UU_RAD, k, j, i)/gdet, 
                                     U(b, m_u.U1_RAD, k, j, i)/gdet, 
                                     U(b, m_u.U2_RAD, k, j, i)/gdet, 
                                     U(b, m_u.U3_RAD, k, j, i)/gdet}; 

            // This is R^{t\mu}
            Real R_t_con[GR_DIM];
            G.raise(R_t_cov, R_t_con, k, j, i, Loci::center);

            //Do g_\mu\nu R^{t\nu} R^{t\mu} and call it invariant_scalar
            Real invariant_scalar = 0.0;
            for(int mu=0; mu<4; ++mu) {
                for(int nu=0; nu<4; ++nu) {
                    invariant_scalar += G.gcov(Loci::center, j, i, mu, nu) * R_t_con[mu] * R_t_con[nu];
                }
            }


            // Isolaring u^t_R^2 in Equation 33 to find u^t_R in Equation 32 from Sadowski et al. 2013.
            // It gives g^{tt}\bar{E}^2 - 2 R^{tt} \bar{E} - 3 invariant_scalar = 0
            // It yields the solution \bar{E} = R^{tt} +- sqrt((R^{tt})^2 + 3 g^tt invariant_scalar) / g^{tt}
            // We are gonna take the negative root since g^{tt} is negative and we want \bar{E} to be positive.
            Real g_con_tt = G.gcon(Loci::center, j, i, 0, 0);
            Real E_bar = (R_t_con[0] - m::sqrt(R_t_con[0]*R_t_con[0] + 3.0 * g_con_tt * invariant_scalar)) / g_con_tt;
            
            //then u^t_R = sqrt(1/8 g^{tt} - 9/(8 E_bar^2) * invariant_scalar)
            Real u_R_t = m::sqrt(0.125 * g_con_tt - 1.125/(E_bar * E_bar) * invariant_scalar);

            // check if any of the values are nan, if so, apply floor to both
            if(!isfinite(E_bar)) E_bar = 1e-30;
            if(!isfinite(u_R_t)) u_R_t = 1e-30;
            // Now calculate the other components of u_R^\mu using Equation 27 from Sadowski et al. 2013, which goes as 
            // R^{\mu\nu} = 4/3 \bar{E} u_R^\mu u_R^\nu + 1/3 \bar{E} g^{\mu\nu}
            // I do have R^{t\mu} and R^{tt}, so I can rearrange the equation to find u_R^\mu as follows:
            // u_R^\mu = (R^{t \nu} - 1/3 \bar{E} g^{t\nu}) / (4/3 \bar{E} u_R^t)
            Real u_R_con[GR_DIM]; 
            for(int nu=0; nu<4; ++nu) {
                u_R_con[nu] = (R_t_con[nu] - (E_bar / 3.0) * G.gcon(Loci::center, j, i, 0, nu)) / ((4.0/3.0) * E_bar * u_R_t);
            }


            // Now calculating the other terms of R^\mu\nu using the same equation
            // However, we don't need to calculate R^ij, since we have R^tt and R^ti 
            Real R_uu[GR_DIM - 1][GR_DIM - 1]; // R^ij
            for(int mu=0; mu<GR_DIM-1; ++mu) {
                for(int nu=0; nu<GR_DIM-1; ++nu) {
                    R_uu[mu][nu] = (4.0/3.0) * E_bar * u_R_con[mu+1] * u_R_con[nu+1] + (E_bar / 3.0) * G.gcon(Loci::center, j, i, mu+1, nu+1);
                }
            } 

            Real R_uu_complete[GR_DIM][GR_DIM] = {
                {R_t_con[0], R_t_con[1], R_t_con[2], R_t_con[3]},
                {R_t_con[1], R_uu[0][0], R_uu[0][1], R_uu[0][2]},
                {R_t_con[2], R_uu[1][0], R_uu[1][1], R_uu[1][2]},
                {R_t_con[3], R_uu[2][0], R_uu[2][1], R_uu[2][2]}
            };

            // Calculate the mixed index R^\mu_\nu
            Real R_mixed[GR_DIM][GR_DIM] = {0.0}; 
            for(int mu=0; mu<GR_DIM; ++mu) {
                for(int nu=0; nu<GR_DIM; ++nu) {
                    for(int sigma=0; sigma<GR_DIM; ++sigma) {
                        R_mixed[mu][nu] +=  R_uu_complete[mu][sigma] * G.gcov(Loci::center, j, i, sigma, nu);
                    }
                }
            }


            // Here we'll need to calculate Kabs depending on the physical processes involved in the radiation field
            Real kappa_a = calc_kabs(0.0, 0.0);
            Real kappa_s = 0.4; // Thomson scattering opacity for electron scattering, for example.
            Real lambda      = calc_lambda(1.e4);

            // Now we are gonna use equation 16 from Mckinney et al. 2014 to calculate the radiation four-force G_\nu, which goes as
            // G^\nu = -kappa_a (R^\mu_\nu u^\nu + \lambda u^\mu) - kappa_s (R^\mu_\alpha u^\alpha + R^\alpha_\beta u_\alpha u^\beta u^\mu)

            Real uvec[3];
            uvec[0] = P(b, m_p.U1, k, j, i);
            uvec[1] = P(b, m_p.U2, k, j, i);
            uvec[2] = P(b, m_p.U3, k, j, i);

            Real ucon[GR_DIM];
            Real ucov[GR_DIM];
            GRMHD::calc_ucon(G, uvec, k, j, i, Loci::center, ucon);
            G.lower(ucon, ucov, k, j, i, Loci::center);

            // Calculating the scattering term

            Real Ruu_scalar = 0.0;
            for(int alpha=0; alpha<GR_DIM; ++alpha) {
                for(int beta=0; beta<GR_DIM; ++beta) {
                    Ruu_scalar += R_mixed[alpha][beta] * ucov[alpha] * ucon[beta];
                }
            }

            Real G_con[GR_DIM] = {0.0};

            for(int mu=0; mu<GR_DIM; ++mu) {
                //Calculate the vector projection V^mu = R^mu_nu * u^nu
                Real Ru_vec = 0.0;
                for(int nu=0; nu<GR_DIM; ++nu) {
                    Ru_vec += R_mixed[mu][nu] * ucon[nu];
                }
                // Absorption Term: G_abs = -kappa_a * ( V^mu + lambda * u^mu )
                Real force_abs = -kappa_a * (Ru_vec) + lambda * ucon[mu];
                //Scattering Term: G_scat = -kappa_s * ( V^mu + S * u^mu )
                Real force_scatt = -kappa_s * (Ru_vec + Ruu_scalar * ucon[mu]);
                // Combine
                G_con[mu] = force_abs + force_scatt;
            }
       
            // I think we need to lower the index here? (ASK BEN)
            Real G_cov[GR_DIM];

            G.lower(G_con, G_cov, k, j, i, Loci::center);
            //printf("G_cov: %e %e %e %e\n", G_cov[0], G_cov[1], G_cov[2], G_cov[3]);
            // G_cov[0]  = 1e-8;
            // G_cov[1]  = 1e-8;
            // G_cov[2]  = 1e-8;
            // G_cov[3]  = 1e-8;

            // // FLUID
            // dUdt(b, m_u.UU, k, j, i) += gdet * G_cov[0];
            // dUdt(b, m_u.U1, k, j, i) += gdet * G_cov[1];
            // dUdt(b, m_u.U2, k, j, i) += gdet * G_cov[2];
            // dUdt(b, m_u.U3, k, j, i) += gdet * G_cov[3];

            // //RADIATION
            // dUdt(b, m_u.UU_RAD, k, j, i) -= gdet * G_cov[0];
            // dUdt(b, m_u.U1_RAD, k, j, i) -= gdet * G_cov[1];
            // dUdt(b, m_u.U2_RAD, k, j, i) -= gdet * G_cov[2];
            // dUdt(b, m_u.U3_RAD, k, j, i) -= gdet * G_cov[3];
        }
    );

    return TaskStatus::complete;
}


// TaskStatus RadM1::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
// {
//     auto pmesh = mdudt->GetMeshPointer();
//     auto pmb0 = mdudt->GetBlockData(0)->GetBlockPointer();

//     const auto& pars = pmb0->packages.Get("RadM1")->AllParams();
//     const Real G0 = pars.Get<Real>("const_G0");
//     const Real G1 = pars.Get<Real>("const_G1");
//     const Real G2 = pars.Get<Real>("const_G2");
//     const Real G3 = pars.Get<Real>("const_G3");

//     // Pack Variables
//     PackIndexMap cons_map;
//     auto dUdt = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
//     const VarMap m_u(cons_map, true);

//     // Get Loop Bounds
//     const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
//     const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
//     const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
//     const IndexRange block = IndexRange{0, dUdt.GetDim(5) - 1};

//     // 4. Parallel Loop
//     pmb0->par_for("add_radM1_source", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
//         KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            
//             const auto& G = dUdt.GetCoords(b);
//             Real gdet = G.gdet(Loci::center, j, i); 

//             // Calculate G_nu (Covariant Components)
//             Real Gnu_lower[GR_DIM] = {0.0};
//             //calc_Gnu(G0, G1, G2, G3, Gnu_lower);

//             // Apply Source to Gas
            
//             // Energy (Component 0)
//             dUdt(b, m_u.UU, k, j, i) += gdet * Gnu_lower[0];

//             // Momentum (Components 1, 2, 3)
//             dUdt(b, m_u.U1, k, j, i) += gdet * Gnu_lower[1];
//             dUdt(b, m_u.U2, k, j, i) += gdet * Gnu_lower[2];
//             dUdt(b, m_u.U3, k, j, i) += gdet * Gnu_lower[3];
//         }
//     );

//     return TaskStatus::complete;
// }