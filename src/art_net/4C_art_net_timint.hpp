// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_TIMINT_HPP
#define FOUR_C_ART_NET_TIMINT_HPP

#include "4C_config.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace Core::LinAlg
{
  class Solver;
}

namespace Discret
{
  class ResultTest;
}

namespace Arteries
{
  /*!
   * \brief time integration for artery network problems
   */

  class TimInt : public Adapter::ArtNet
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimInt(std::shared_ptr<Core::FE::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        Core::IO::DiscretizationWriter& output);


    //! initialize time integration
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) override;

    //! get discretization
    std::shared_ptr<Core::FE::Discretization> discretization() override { return discret_; }

    double dt() const override { return dta_; }

    double time() const { return time_; }
    int step() const { return step_; }

    int itemax() const { return params_.get<int>("max nonlin iter steps"); }

    void output(
        bool CoupledTo3D, std::shared_ptr<Teuchos::ParameterList> CouplingParams) override = 0;

    /*!
    \brief start time loop for startingalgo, normal problems and restarts

    */
    void integrate(
        bool CoupledTo3D, std::shared_ptr<Teuchos::ParameterList> CouplingParams) override;

    /*!
    \brief prepare the loop

    */
    void prepare_time_loop() override;

    /*!
    \brief Do time integration (time loop)

    */
    void time_loop(bool CoupledTo3D, std::shared_ptr<Teuchos::ParameterList> CouplingTo3DParams);

    //! set the initial field on the artery discretization
    virtual void set_initial_field(
        const Inpar::ArtDyn::InitialField init,  //!< type of initial field
        const int startfuncno                    //!< number of spatial function
    )
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
    }


    /// setup the variables to do a new time step
    void prepare_time_step() override;

    /// setup the variables to do a new time step
    void prepare_linear_solve() override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
    }

    /// setup the variables to do a new time step
    void assemble_mat_and_rhs() override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
    }

    /// direct access to system matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      return std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(sysmat_);
    };

    //! right-hand side alias the dynamic force residual
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() const override { return rhs_; }

    //! iterative update of primary variable
    void update_iter(const std::shared_ptr<const Core::LinAlg::Vector<double>> inc) override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
      return;
    }

    // get solution vector
    std::shared_ptr<const Core::LinAlg::Vector<double>> pressurenp() const override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
      return nullptr;
    }
    /*!
    \brief solve linearised artery and bifurcation

    */
    void solve(std::shared_ptr<Teuchos::ParameterList> CouplingTo3DParams) override = 0;

    void solve_scatra() override = 0;

    //! is output needed for the current time step?
    bool do_output() { return ((step_ % upres_ == 0) or (step_ % uprestart_ == 0)); };

    // set solve scatra flag
    void set_solve_scatra(const bool solvescatra) override
    {
      solvescatra_ = solvescatra;
      return;
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    std::shared_ptr<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() const override
    {
      return dbcmaps_;
    }

    // create field test
    std::shared_ptr<Core::Utils::ResultTest> create_field_test() override = 0;

   protected:
    //! @name general algorithm parameters
    //! arterial network discretization
    std::shared_ptr<Core::FE::Discretization> discret_;
    //! linear solver
    std::shared_ptr<Core::LinAlg::Solver> solver_;
    const Teuchos::ParameterList& params_;
    Core::IO::DiscretizationWriter& output_;
    //! the processor ID from the communicator
    int myrank_;

    /// (standard) system matrix
    std::shared_ptr<Core::LinAlg::SparseOperator> sysmat_;

    /// maps for extracting Dirichlet and free DOF sets
    std::shared_ptr<Core::LinAlg::MapExtractor> dbcmaps_;

    /// rhs: right hand side vector
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs_;

    /// (scatra) system matrix
    std::shared_ptr<Core::LinAlg::SparseOperator> scatra_sysmat_;

    /// rhs: right hand side vector of scatra
    std::shared_ptr<Core::LinAlg::Vector<double>> scatra_rhs_;

    //! @name time step sizes
    double dta_;
    double dtp_;

    //! @name cpu-time measures
    double dtele_;
    double dtsolve_;
    //@}

    //! @name time stepping variables
    double time_;     ///< physical time
    int step_;        ///< timestep
    int stepmax_;     ///< maximal number of timesteps
    double maxtime_;  ///< maximal physical computation time
    //@}

    //! @name restart variables
    int uprestart_;
    int upres_;
    //@}

    bool solvescatra_;
    const int linsolvernumber_;

    //!
    bool coupledTo3D_;
    //@}


  };  // class TimInt
}  // namespace Arteries



FOUR_C_NAMESPACE_CLOSE

#endif
