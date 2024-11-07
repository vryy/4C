// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_HPP
#define FOUR_C_ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_HPP

#include "4C_config.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_adapter_scatra_fluid_coupling_algorithm.hpp"
#include "4C_coupling_adapter.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Adapter
{
  /// basis coupling algorithm for scalar transport with Navier-Stokes on moving meshes
  /*!

    Base class for scalar transport problems coupled to Navier-Stokes velocity field on
    deforming meshes.
    Derives from ScaTraFluidCouplingAlgorithm and AleBaseAlgorithm and establishes
    the fluid-ale coupling.
    Different application coupling algorithms will inherit from this base class
    (at the moment only electrochemistry applications).

    \author gjb
    \date 05/09
   */
  class ScaTraFluidAleCouplingAlgorithm : public Adapter::ScaTraFluidCouplingAlgorithm,
                                          public Adapter::AleBaseAlgorithm
  {
   public:
    /// constructor using a Epetra_Comm
    ScaTraFluidAleCouplingAlgorithm(const Epetra_Comm& comm,  ///< communicator
        const Teuchos::ParameterList& prbdyn,                 ///< problem-specific parameters
        const std::string condname,  ///< name of condition that defines fluid-ale coupling
        const Teuchos::ParameterList& solverparams);


    /// setup
    void setup() override;

    /// init
    void init() override;

    /// read restart data (pure virtual)
    void read_restart(int step  ///< step number where the calculation is continued
        ) override = 0;

    /// solve fluid-ale
    virtual void fluid_ale_nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel, const bool pseudotransient);

    /// access to ale field
    const std::shared_ptr<Adapter::AleFluidWrapper>& ale_field() { return ale_; }

   protected:
    //! @name Transfer helpers

    /// field transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv) const;

    /// field transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv) const;

    /// interface transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_ale(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv) const;

    /// interface transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_ale(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv) const;

   private:
    /// ALE-fluid wrapper
    std::shared_ptr<AleFluidWrapper> ale_;

    /// coupling of fluid and ale (whole field)
    std::shared_ptr<Coupling::Adapter::Coupling> coupfa_;

    /// coupling of fluid and ale (interface only)
    std::shared_ptr<Coupling::Adapter::Coupling> icoupfa_;

    /// condition name
    const std::string condname_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
