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
    auto m = Metadata(flags_prim);
    pkg->AddField("prims.u_rad", m);
    pkg->AddField("prims.uvec_rad", m);
    m = Metadata(flags_cons);
    pkg->AddField("cons.u_rad", m);
    pkg->AddField("cons.uvec_rad", m);
    
    //Right now, to execute the torus problem with radM1, we need to initialize the radiation primitives in fm_torus.cpp (this is stupid)
    //I think this should be moved to RadM1 method, maybe call an initialization method like a task straight after initializing the torus?
    //Especially because we'll want to initialize the radiation field for other problems. 

    //This method should allow you to add source terms to both plasma and radiation variables separately.
    pkg->AddSource = RadM1::AddSource;

    return pkg;
}


TaskStatus RadM1::AddSource(MeshData<Real> *md, MeshData<Real> *mdudt, IndexDomain domain)
{
    auto pmesh = mdudt->GetMeshPointer();
    auto pmb0 = mdudt->GetBlockData(0)->GetBlockPointer();

    const auto& pars = pmb0->packages.Get("RadM1")->AllParams();
    const Real G0 = pars.Get<Real>("const_G0");
    const Real G1 = pars.Get<Real>("const_G1");
    const Real G2 = pars.Get<Real>("const_G2");
    const Real G3 = pars.Get<Real>("const_G3");

    // Pack Variables
    PackIndexMap cons_map;
    auto dUdt = mdudt->PackVariables(std::vector<MetadataFlag>{Metadata::Conserved}, cons_map);
    const VarMap m_u(cons_map, true);

    // Get Loop Bounds
    const IndexRange ib = mdudt->GetBoundsI(IndexDomain::interior);
    const IndexRange jb = mdudt->GetBoundsJ(IndexDomain::interior);
    const IndexRange kb = mdudt->GetBoundsK(IndexDomain::interior);
    const IndexRange block = IndexRange{0, dUdt.GetDim(5) - 1};

    // 4. Parallel Loop
    pmb0->par_for("add_radM1_source", block.s, block.e, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int& b, const int &k, const int &j, const int &i) {
            
            const auto& G = dUdt.GetCoords(b);
            Real gdet = G.gdet(Loci::center, j, i); 

            // Calculate G_nu (Covariant Components)
            Real Gnu_lower[GR_DIM];
            calc_Gnu(G0, G1, G2, G3, Gnu_lower);

            // Apply Source to Gas
            
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