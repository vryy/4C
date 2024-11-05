// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FS3I_BIOFILM_FSI_HPP
#define FOUR_C_FS3I_BIOFILM_FSI_HPP


#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_fs3i_partitioned_1wc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Adapter
{
  class AleFsiWrapper;
  class StructureBio;
  class FSIStructureWrapper;
}  // namespace Adapter

namespace ALE
{
  class AleBaseAlgorithm;
}

namespace FS3I
{
  class BiofilmFSI : public PartFS3I1Wc
  {
   public:
    BiofilmFSI(const Epetra_Comm& comm);

    void init() override;

    void setup() override;

    void timeloop() override;

    void inner_timeloop();

    //! information transfer FSI -> ScaTra
    void set_fsi_solution();

    void compute_interface_vectors(Core::LinAlg::Vector<double>& idispnp_,
        Core::LinAlg::Vector<double>& iveln_,
        std::shared_ptr<Core::LinAlg::Vector<double>> struidispnp_,
        Core::LinAlg::Vector<double>& struiveln_);

    std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_ale(
        Core::LinAlg::Vector<double>& iv) const;

    std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        Core::LinAlg::Vector<double>& iv) const;

    /// field transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_struct_field(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv) const;

    /// field transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_struct_field(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv) const;

    /// interface transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> struct_to_ale(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv) const;

    /// interface transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> struct_to_ale(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv) const;

    /// solve fluid-ale
    virtual void fluid_ale_solve();

    /// solve structure-ale
    virtual void struct_ale_solve();

    void update_and_output();

    const Epetra_Comm& comm() { return comm_; }

    void vec_to_scatravec(Core::FE::Discretization& scatradis, Core::LinAlg::Vector<double>& vec,
        Core::LinAlg::MultiVector<double>& scatravec);

    void struct_gmsh_output();

    void fluid_gmsh_output();

   private:
    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    /// coupling of fluid and ale (interface only)
    std::shared_ptr<Coupling::Adapter::Coupling> icoupfa_;

    /// coupling of fluid and ale (whole field)
    std::shared_ptr<Coupling::Adapter::Coupling> coupfa_;

    /// coupling of structure and ale (interface only)
    std::shared_ptr<Coupling::Adapter::Coupling> icoupsa_;

    /// coupling of structure and ale (whole field)
    std::shared_ptr<Coupling::Adapter::Coupling> coupsa_;

    // std::shared_ptr< ::Adapter::FSIStructureWrapper> structure_;

    // std::shared_ptr<ALE::AleBaseAlgorithm> ale_;
    std::shared_ptr<Adapter::AleFsiWrapper> ale_;

    //    // total flux at the interface overall the InnerTimeloop
    //    std::shared_ptr<Core::LinAlg::MultiVector<double>> flux;
    //
    //    // total flux at the structure interface overall the InnerTimeloop
    //    std::shared_ptr<Core::LinAlg::MultiVector<double>> struflux;

    std::shared_ptr<Core::LinAlg::Vector<double>> norminflux_;

    std::shared_ptr<Core::LinAlg::Vector<double>> lambda_;
    std::shared_ptr<Core::LinAlg::Vector<double>> normtraction_;
    std::shared_ptr<Core::LinAlg::Vector<double>> tangtractionone_;
    std::shared_ptr<Core::LinAlg::Vector<double>> tangtractiontwo_;

    std::vector<double> nvector_;

    // coefficients used in the calculation of the displacement due to growth
    // fluxcoef_ multiply the scalar influx at the interface,
    // while normforcecoef_, tangoneforcecoef_ and tangtwoforcecoef_  multiply forces
    // in the normal and in the two tangential directions at the interface
    double fluxcoef_;
    double normforceposcoef_;
    double normforcenegcoef_;
    double tangoneforcecoef_;
    double tangtwoforcecoef_;

    //// growth time parameters

    // number of steps
    int nstep_bio_;

    // current step
    int step_bio_;

    // time step size
    double dt_bio_;

    // total time of the outer loop
    double time_bio_;


    //// scatra and fsi time parameters

    // number of steps
    int nstep_fsi_;

    // current step
    int step_fsi_;

    // time step size
    double dt_fsi_;

    // total time of the inner loop
    double time_fsi_;

    // maximum time
    double maxtime_fsi_;

    // total time
    double time_;

    /// fluid interface displacement at time t^{n}
    std::shared_ptr<Core::LinAlg::Vector<double>> idispn_;

    /// fluid interface displacement at time t^{n+1}
    std::shared_ptr<Core::LinAlg::Vector<double>> idispnp_;

    /// fluid velocity at interface (always zero!)
    std::shared_ptr<Core::LinAlg::Vector<double>> iveln_;

    /// structure interface displacement at time t^{n}
    std::shared_ptr<Core::LinAlg::Vector<double>> struidispn_;

    /// structure interface displacement at time t^{n+1}
    std::shared_ptr<Core::LinAlg::Vector<double>> struidispnp_;

    /// structure velocity at interface (always zero!)
    std::shared_ptr<Core::LinAlg::Vector<double>> struiveln_;

    /// total structure displacement due to growth
    std::shared_ptr<Core::LinAlg::Vector<double>> struct_growth_disp_;

    /// total fluid displacement due to growth
    std::shared_ptr<Core::LinAlg::Vector<double>> fluid_growth_disp_;

    /// total scatra structure displacement due to growth
    std::shared_ptr<Core::LinAlg::MultiVector<double>> scatra_struct_growth_disp_;

    /// total scatra fluid displacement due to growth
    std::shared_ptr<Core::LinAlg::MultiVector<double>> scatra_fluid_growth_disp_;
  };

}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
