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


std::shared_ptr<KHARMAPackage> RadM1::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    auto pkg = std::make_shared<KHARMAPackage>("RadM1");
    Params &params = pkg->AllParams();

    //Place holder parameters for the radiation four-force, for testing purposes.  
    Real const_G0 = pin->GetOrAddReal("RadM1", "const_G0", 0.0);
    params.Add("const_G0", const_G0);
    Real const_G1 = pin->GetOrAddReal("RadM1", "const_G1", 0.0);
    params.Add("const_G1", const_G1);
    Real const_G2 = pin->GetOrAddReal("RadM1", "const_G2", 0.0);
    params.Add("const_G2", const_G2);
    Real const_G3 = pin->GetOrAddReal("RadM1", "const_G3", 0.0);
    params.Add("const_G3", const_G3);

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

    auto flags_prim_vec = flags_prim;
    flags_prim_vec.push_back(Metadata::Vector);

    auto flags_cons_vec = flags_cons;
    flags_cons_vec.push_back(Metadata::Vector);

    std::vector<int> s_vector({NVEC});

    auto m_prim_vector = Metadata(flags_prim_vec, s_vector);
    pkg->AddField("prims.uvec_rad", m_prim_vector);

    auto m_cons_vector = Metadata(flags_cons_vec, s_vector);
    pkg->AddField("cons.uvec_rad", m_cons_vector);
    
    //Right now, to execute the torus problem with radM1, we need to initialize the radiation primitives in fm_torus.cpp (this is stupid) (ASK BEN)
    //I think this should be moved to RadM1 method, maybe call an initialization method like a task straight after initializing the torus?
    //Especially because we'll want to initialize the radiation field for other problems. 

    //This method should allow you to add source terms to both plasma and radiation variables separately.
    pkg->AddSource = RadM1::AddSource;

    return pkg;
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
            if(isnan(E_bar)){
                printf("R_tcon[0]: %e, g_con_tt: %e, invariant_scalar: %e\n", R_t_con[0], g_con_tt, invariant_scalar);
            }
            //Floor applied in Mckinney et al 2013? (ASK BEN ABOUT FLOORS)
            if(E_bar < 1e-150){
                E_bar = 1e-150;
            }

            //then u^t_R = sqrt(1/8 g^{tt} - 9/(8 E_bar^2) * invariant_scalar)
            Real u_R_t = m::sqrt(0.125 * g_con_tt - 1.125/(E_bar * E_bar) * invariant_scalar);
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

            // Calculate the mixed index R^\mu_\nu
            Real R_mixed[GR_DIM][GR_DIM]; 
            for(int mu=0; mu<GR_DIM; ++mu) {
                for(int nu=0; nu<GR_DIM; ++nu) {
                    R_mixed[mu][nu] = 0.0;
                    for(int sigma=0; sigma<GR_DIM; ++sigma) {
                        R_mixed[mu][nu] += G.gcov(Loci::center, j, i, mu, sigma) * ((sigma == 0) ? R_t_con[nu] : R_uu[sigma-1][nu-1]);
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
            Real scattering_term[GR_DIM] = {0};
            for(int mu=0; mu<GR_DIM; ++mu) {
                for(int alpha=0; alpha<GR_DIM; ++alpha) {
                    for (int beta=0; beta<GR_DIM; ++beta) {
                        scattering_term[mu] += R_mixed[mu][alpha] * ucov[alpha] + R_mixed[alpha][beta] * ucov[alpha] * ucov[beta] * ucon[mu];
                    }
                }
                scattering_term[mu] *= - kappa_s;
            }
            // Calculating the absorption term
            Real absorption_term[GR_DIM] = {0};
            for(int mu=0; mu<GR_DIM; ++mu) {
                for(int nu=0; nu<GR_DIM; ++nu) {
                    absorption_term[mu] += kappa_a * R_mixed[mu][nu] * ucon[nu] + lambda * ucon[mu];
                }
                absorption_term[mu] *= -1;
            }

            Real G_con[GR_DIM]; // G^\nu

            for (int mu=0; mu<GR_DIM; ++mu) {
                G_con[mu] = absorption_term[mu] + scattering_term[mu];
            }
       
            // I think we need to lower the index here? (ASK BEN)
            Real G_cov[GR_DIM];

            
            for(int nu=0; nu<4; ++nu) {
                for(int mu=0; mu<4; ++mu) {
                    G_cov[nu] += G.gcov(Loci::center, j, i, nu, mu) * G_con[mu];
                }
            }

            
            // FLUID
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