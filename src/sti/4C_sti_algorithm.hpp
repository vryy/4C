// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STI_ALGORITHM_HPP
#define FOUR_C_STI_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  template <typename T>
  class Vector;
}

namespace Adapter
{
  class ScaTraBaseAlgorithm;
}  // namespace Adapter

namespace ScaTra
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace ScaTra

namespace STI
{
  //! monolithic algorithm for scatra-thermo interaction
  class Algorithm : public Adapter::AlgorithmBase
  {
   public:
    //! return counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    const unsigned& iter() const { return iter_; };

    //! read restart data
    void read_restart(int step  //! time step for restart
        ) override;

    //! access scatra time integrator
    std::shared_ptr<ScaTra::ScaTraTimIntImpl> scatra_field() const
    {
      return scatra_->scatra_field();
    };

    //! access thermo time integrator
    std::shared_ptr<ScaTra::ScaTraTimIntImpl> thermo_field() const
    {
      return thermo_->scatra_field();
    };

    //! time loop
    void time_loop();

   protected:
    //! constructor
    explicit Algorithm(const Epetra_Comm& comm,  //! communicator
        const Teuchos::ParameterList& stidyn,    //! parameter list for scatra-thermo interaction
        const Teuchos::ParameterList&
            scatradyn,  //! scalar transport parameter list for scatra and thermo fields
        const Teuchos::ParameterList&
            solverparams_scatra,  //! solver parameter list for scatra field
        const Teuchos::ParameterList&
            solverparams_thermo  //! solver parameter list for thermo field
    );

    //! prepare time step
    void prepare_time_step() override;

    //! pass scatra degrees of freedom to thermo discretization
    void transfer_scatra_to_thermo(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> scatra  //!< scatra state vector
    ) const;

    //! pass thermo degrees of freedom to scatra discretization
    void transfer_thermo_to_scatra(
        const std::shared_ptr<const Core::LinAlg::Vector<double>> thermo  //!< thermo state vector
    ) const;

    //! scatra time integrator
    std::shared_ptr<Adapter::ScaTraBaseAlgorithm> scatra_;

    //! thermo time integrator
    std::shared_ptr<Adapter::ScaTraBaseAlgorithm> thermo_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    std::shared_ptr<ScaTra::MeshtyingStrategyS2I> strategyscatra_;

    //! meshtying strategy for scatra-scatra interface coupling on thermo discretization
    std::shared_ptr<ScaTra::MeshtyingStrategyS2I> strategythermo_;

    //! input parameters for scatra and thermo fields
    std::shared_ptr<Teuchos::ParameterList> fieldparameters_;

    //! counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling iterations
    //! (partitioned algorithm)
    unsigned int iter_;

    //! maximum number of Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    unsigned int itermax_;

    //! tolerance for Newton-Raphson iteration (monolithic algorithm) or outer coupling iteration
    //! (partitioned algorithm)
    double itertol_;

    //! input parameters for scatra-thermo interaction
    std::shared_ptr<Teuchos::ParameterList> stiparameters_;

    //! timer for Newton-Raphson iteration (monolithic algorithm) or outer coupling iteration
    //! (partitioned algorithm)
    std::shared_ptr<Teuchos::Time> timer_;

   private:
    //! modify field parameters for thermo field
    void modify_field_parameters_for_thermo_field();

    //! output solution to screen and files
    void output() override;

    //! evaluate time step using Newton-Raphson iteration (monolithic algorithm) or outer coupling
    //! iteration (partitioned algorithm)
    virtual void solve() = 0;

    //! update scatra and thermo fields after time step evaluation
    void update() override;
  };  // class Algorithm : public Adapter::AlgorithmBase
}  // namespace STI
FOUR_C_NAMESPACE_CLOSE

#endif
