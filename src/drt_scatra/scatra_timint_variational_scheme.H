/*----------------------------------------------------------------------*/
/*! \file

\brief  Time integration for variational formulation problems using different schemes
\level 2

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_VARIATIONAL_SCHEME_H
#define SCATRA_TIMINT_VARIATIONAL_SCHEME_H

#include "scatra_timint_variational.H"
#include "scatra_timint_ost.H"
//!#include "scatra_timint_bdf2.H"
//!#include "scatra_timint_genalpha.H"
//!#include "scatra_timint_stat.H"

namespace SCATRA
{
  class TimIntVariationalOST : public TimIntVariational, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    TimIntVariationalOST(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<LINALG::Solver> solver,                        //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,                //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
        const int probnum = 0                              //!< global problem number
    );

    //! Destructor
    virtual ~TimIntVariationalOST();

    //! initialize time integration scheme
    virtual void Init();

    //! setup time integration scheme
    virtual void Setup();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

    //! Indicates if the functional is evaluated fully implicit or semi implicit.
    bool IsSemImplicitFunctional() { return semImplicitFunctional_; }

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! initialization procedure prior to evaluation of first time step
    virtual void PrepareFirstTimeStep();

   protected:
    //! explicit predictor for nonlinear solver
    void ExplicitPredictor() const;

    /// add parameters specific for time-integration scheme
    virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);


    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    //! add parameters depending on the problem, i.e., loma, level-set, ...
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params  //!< parameter list
    );

   private:
    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for potentially skipping computation of initial time derivative
    bool semImplicitFunctional_;

  };  // class TimIntVariationalOST
}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_VARIATIONAL_SCHEME_H
