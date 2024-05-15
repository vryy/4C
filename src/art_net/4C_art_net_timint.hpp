/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for arterial network (time) integration.


\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_TIMINT_HPP
#define FOUR_C_ART_NET_TIMINT_HPP

#include "4C_config.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_lib_discret.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace CORE::LINALG
{
  class Solver;
}

namespace DRT
{
  class ResultTest;
}

namespace ART
{
  /*!
   * \brief time integration for artery network problems
   */

  class TimInt : public ADAPTER::ArtNet
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimInt(Teuchos::RCP<DRT::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
        IO::DiscretizationWriter& output);


    //! initialize time integration
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) override;

    //! get discretization
    Teuchos::RCP<DRT::Discretization> Discretization() override { return discret_; }

    double Dt() const override { return dta_; }

    double Time() const { return time_; }
    int Step() const { return step_; }

    int Itemax() const { return params_.get<int>("max nonlin iter steps"); }

    void Output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) override = 0;

    /*!
    \brief start time loop for startingalgo, normal problems and restarts

    */
    void Integrate(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) override;

    /*!
    \brief prepare the loop

    */
    void PrepareTimeLoop() override;

    /*!
    \brief Do time integration (time loop)

    */
    void TimeLoop(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams);

    //! set the initial field on the artery discretization
    virtual void SetInitialField(const INPAR::ARTDYN::InitialField init,  //!< type of initial field
        const int startfuncno  //!< number of spatial function
    )
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
    }


    /// setup the variables to do a new time step
    void PrepareTimeStep() override;

    /// setup the variables to do a new time step
    void PrepareLinearSolve() override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
    }

    /// setup the variables to do a new time step
    void AssembleMatAndRHS() override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
    }

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override
    {
      return Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(sysmat_);
    };

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> RHS() const override { return rhs_; }

    //! iterative update of primary variable
    void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc) override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
      return;
    }

    // get solution vector
    Teuchos::RCP<const Epetra_Vector> Pressurenp() const override
    {
      // each artery integration should overwrite this if used
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    }
    /*!
    \brief solve linearised artery and bifurcation

    */
    void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams) override = 0;

    void SolveScatra() override = 0;

    //! is output needed for the current time step?
    bool DoOutput() { return ((step_ % upres_ == 0) or (step_ % uprestart_ == 0)); };

    // set solve scatra flag
    void SetSolveScatra(const bool solvescatra) override
    {
      solvescatra_ = solvescatra;
      return;
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() const override
    {
      return dbcmaps_;
    }

    // create field test
    Teuchos::RCP<CORE::UTILS::ResultTest> CreateFieldTest() override = 0;

   protected:
    //! @name general algorithm parameters
    //! arterial network discretization
    Teuchos::RCP<DRT::Discretization> discret_;
    //! linear solver
    Teuchos::RCP<CORE::LINALG::Solver> solver_;
    const Teuchos::ParameterList& params_;
    IO::DiscretizationWriter& output_;
    //! the processor ID from the communicator
    int myrank_;

    /// (standard) system matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> sysmat_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps_;

    /// rhs: right hand side vector
    Teuchos::RCP<Epetra_Vector> rhs_;

    /// (scatra) system matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_sysmat_;

    /// rhs: right hand side vector of scatra
    Teuchos::RCP<Epetra_Vector> scatra_rhs_;

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
}  // namespace ART



FOUR_C_NAMESPACE_CLOSE

#endif
