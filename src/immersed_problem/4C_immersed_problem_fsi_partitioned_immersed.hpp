// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IMMERSED_PROBLEM_FSI_PARTITIONED_IMMERSED_HPP
#define FOUR_C_IMMERSED_PROBLEM_FSI_PARTITIONED_IMMERSED_HPP

#include "4C_config.hpp"

#include "4C_fsi_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  class PartitionedImmersed : public Partitioned
  {
   public:
    //! constructor
    explicit PartitionedImmersed(const Epetra_Comm& comm);

    //! setup this object
    void setup() override;

    //! overrides method of base class.
    void setup_coupling(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm) override;

    //! call the time loop of the base class
    void timeloop(const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface) override
    {
      FSI::Partitioned::timeloop(interface);
    };

    //! override version of fsi partitioned
    void extract_previous_interface_solution() override;

    //! Implement pure virtual functions (again overloaded by corresponding partitioned subclass in
    //! immersed_problem)
    void fsi_op(const Core::LinAlg::Vector<double>& x, Core::LinAlg::Vector<double>& F,
        const FillType fillFlag) override
    {
      return;
    };

    //! empty; overridden in sub class
    std::shared_ptr<Core::LinAlg::Vector<double>> fluid_op(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp, const FillType fillFlag) override
    {
      return nullptr;
    };

    //! empty; overridden in sub class
    std::shared_ptr<Core::LinAlg::Vector<double>> struct_op(
        std::shared_ptr<Core::LinAlg::Vector<double>> iforce, const FillType fillFlag) override
    {
      return nullptr;
    };

    //! empty; overridden in sub class
    std::shared_ptr<Core::LinAlg::Vector<double>> initial_guess() override { return nullptr; };


  };  // class PartitionedImmersed
}  // namespace FSI


FOUR_C_NAMESPACE_CLOSE

#endif
