/*----------------------------------------------------------------------*/
/*! \file

\brief  Time integration for variational formulation of chemical diffusion
\level 2

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_VARIATIONAL_H
#define SCATRA_TIMINT_VARIATIONAL_H

#include "scatra_timint_ost.H"
//!#include "scatra_timint_bdf2.H"
//!#include "scatra_timint_genalpha.H"
//!#include "scatra_timint_stat.H"

//
// To save error files
#include "../drt_io/io_control.H"

namespace SCATRA
{
  class TimIntVariational : public virtual ScaTraTimIntImpl
  {
   public:
    //! Standard Constructor
    TimIntVariational(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<LINALG::Solver> solver,                     //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,             //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,        //!< supplementary parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output,           //!< output writer
        const int probnum = 0                                    //!< global problem number
    );

    //! Destructor
    virtual ~TimIntVariational();

    //! initialize time integration scheme
    virtual void Init();

    //! setup time integration scheme
    virtual void Setup();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

    //! set the initial scalar field phi
    virtual void SetInitialField(const INPAR::SCATRA::InitialField init,  //!< type of initial field
        const int startfuncno  //!< number of spatial function
    );

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/
    //! calculate error compared to analytical solution
    void EvaluateErrorComparedToAnalyticalSol();

    //! return elch parameter list
    const Teuchos::RCP<const Teuchos::ParameterList> VarParameterList() const { return varparams_; }

    /*========================================================================*/
    //! @name Preconditioning
    /*========================================================================*/

    void SetupSplitter();


    /*========================================================================*/
    //! @name general framework //TODO try to make them static to access them also from the element
    //! and simplify code
    /*========================================================================*/
    void AnalyticSolution_SeriesErrorFnt(double& Conc,  //!< Concentration
        const int series_end,                           //!< stopping number for the series
        const double x,                                 //!< evaluation position in x
        const double t,                                 //!< evaluation time
        const double L,                                 //!< length of the domain
        const double D,                                 //!< Diffusion parameter
        const double c_0,                               //!< initial/reference concentration
        const double Cext                               //!< concentration at Dirichlet nodes
    );

    void PostProcess_ChemPot(double& ChemPot,  //!< chemical potential
        const double conct,                    //!< concentration
        DRT::Node* node  //!< Node: to extract corresponding post-process relation
    );
    /*========================================================================*/
    //! @name scalar degrees of freedom and related
    /*========================================================================*/

    //! return scalar field phi at time 0
    Teuchos::RCP<Epetra_Vector> Phi0() const { return phi0_; }

   protected:
    //! explicit predictor for nonlinear solver
    void ExplicitPredictor() const;

    /// add parameters specific for time-integration scheme
    virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);

    //! Apply Dirichlet boundary conditions on provided state vector
    void ApplyDirichletBC(const double time,  //!< evaluation time
        Teuchos::RCP<Epetra_Vector> phinp,    //!< transported scalar(s) (may be = null)
        Teuchos::RCP<Epetra_Vector> phidt     //!< first time derivative (may be = null)
    );

    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    //! add parameters depending on the problem, i.e., loma, level-set, ...
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params  //!< parameter list
    );

    /*========================================================================*/
    //! @name query and output
    /*========================================================================*/
    //! write extra output vectors
    virtual void OutputState();

    virtual void SaveError2File(const double RelError_Conc,  // !< Relative error for Concentration
        const double RelError_ChemPot  // !< Relative error for Chemical potential
    );

    /*========================================================================*/
    //! @name Preconditioning
    /*========================================================================*/
    //! Parameter check for diffusion-conduction formulation
    void ValidParameterDiffCond();
    //! initialize meshtying strategy (including standard case without meshtying)
    void CreateMeshtyingStrategy();

    /*========================================================================*/
    //! @name Variational variables
    /*========================================================================*/

    //! the parameter list for elch problems
    Teuchos::RCP<const Teuchos::ParameterList> varparams_;

    /*========================================================================*/
    //! @name general framework variables
    /*========================================================================*/

    //! initialize meshtying strategy (including standard case without meshtying)
    virtual void CreateScalarHandler();

    //! phi at time 0
    Teuchos::RCP<Epetra_Vector> phi0_;

    //! To be used when analytic solution is known to export to paraview
    Teuchos::RCP<Epetra_Vector> phiAnalytic_;

  };  // class TimIntVariational

  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Helper class for managing different number of degrees of freedom per node
   */
  class ScalarHandlerVar : public ScalarHandler
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScalarHandlerVar();

    //! Destructor
    virtual ~ScalarHandlerVar(){};

    //! initialize time integration
    virtual void Setup(const ScaTraTimIntImpl* const scatratimint);

    /*========================================================================*/
    //! @name Access and Query methods
    /*========================================================================*/

    //! return maximum number of transported scalars per node (not including potential and current
    //! density)
    virtual int NumScal() const
    {
      CheckIsSetup();
      return *(numscal_.rbegin());
    }

    //! return maximum number of transported scalars per node (not including potential and current
    //! density)
    virtual int NumScalInCondition(const DRT::Condition& condition,
        const Teuchos::RCP<const DRT::Discretization>& discret) const;

    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
   protected:
    //! number of transported scalars (without potential and current density)
    std::set<int> numscal_;

  };  // class ScalarHandlerVar

}  // namespace SCATRA
#endif  // #ifndef SCATRA_TIMINT_VARIATIONAL_H
