// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_POROMULTI_HPP
#define FOUR_C_SCATRA_TIMINT_POROMULTI_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_bdf2.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_ost.hpp"
#include "4C_scatra_timint_stat.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class ScaTraTimIntPoroMulti : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    ScaTraTimIntPoroMulti(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    // -----------------------------------------------------------------
    // general methods
    // -----------------------------------------------------------------

    /// initialize algorithm
    void init() override;

    //! update the solution after convergence of the nonlinear iteration.
    void update() override {};

    //! set the nodal L2-flux
    virtual void set_l2_flux_of_multi_fluid(
        std::shared_ptr<const Core::LinAlg::MultiVector<double>> multiflux);

    //! set solution field of the multiphase fluid
    virtual void set_solution_field_of_multi_fluid(
        std::shared_ptr<const Core::LinAlg::Vector<double>> phinp_fluid,
        std::shared_ptr<const Core::LinAlg::Vector<double>> phin_fluid);

    //! set the velocity field (zero or field by function)
    virtual void set_velocity_field(const int nds)
    {
      FOUR_C_THROW(
          "set_velocity_field(...) cannot be used for transport within a multiphase porous medium!"
          " Use SetSolutionFields(...) instead!");
    };

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void set_velocity_field(std::shared_ptr<const Core::LinAlg::Vector<double>>
                                        convvel,  //!< convective velocity/press. vector
        std::shared_ptr<const Core::LinAlg::Vector<double>> acc,    //!< acceleration vector
        std::shared_ptr<const Core::LinAlg::Vector<double>> vel,    //!< velocity vector
        std::shared_ptr<const Core::LinAlg::Vector<double>> fsvel,  //!< fine-scale velocity vector
        const int nds,  //!< number of the dofset the velocity/pressure state belongs to
        const bool setpressure =
            false  //!< flag whether the fluid pressure needs to be known for the scatra
    )
    {
      FOUR_C_THROW(
          "set_velocity_field(...) cannot be used for transport within a multiphase porous medium!"
          " Use SetSolutionFields(...) instead!");
    };

    void collect_runtime_output_data() override;

    //! add parameters depending on the problem
    void add_problem_specific_parameters_and_vectors(
        Teuchos::ParameterList& params  //!< parameter list
        ) override;

   protected:
    //! do we employ L2-projection for reconstruction of velocity field
    bool L2_projection_;
  };


  class ScaTraTimIntPoroMultiOST : public ScaTraTimIntPoroMulti, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiOST(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

  };  // class TimIntPoroMultiOST


  class ScaTraTimIntPoroMultiBDF2 : public ScaTraTimIntPoroMulti, public TimIntBDF2
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiBDF2(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

  };  // class TimIntPoroMultiBDF2


  class ScaTraTimIntPoroMultiGenAlpha : public ScaTraTimIntPoroMulti, public TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiGenAlpha(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

  };  // class TimIntPoroMultiGenAlpha


  class ScaTraTimIntPoroMultiStationary : public ScaTraTimIntPoroMulti, public TimIntStationary
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiStationary(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

  };  // class TimIntPoroMultiStationary
}  // namespace ScaTra



FOUR_C_NAMESPACE_CLOSE

#endif
