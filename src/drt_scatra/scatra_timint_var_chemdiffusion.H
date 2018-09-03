/*!----------------------------------------------------------------------
\file scatra_timint_var_chemdiffusion.H

\brief  Time integration for variational formulation of chemical diffusion,
    it inherits directly from scatra_timint_ost.
\level 2

<pre>
\maintainer Jorge De Anda Salazar
            deanda@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_VAR_CHEMDIFFUSIONOST_H
#define SCATRA_TIMINT_VAR_CHEMDIFFUSIONOST_H

#include "scatra_timint_ost.H"
//!#include "scatra_timint_bdf2.H"
//!#include "scatra_timint_genalpha.H"
//!#include "scatra_timint_stat.H"

namespace SCATRA
{
  class TimIntVarChemDiffusionOST : public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    TimIntVarChemDiffusionOST(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<LINALG::Solver> solver,                             //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,                     //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
        const int probnum = 0                              //!< global problem number
    );

    //! Destructor
    virtual ~TimIntVarChemDiffusionOST();

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

    //! Indicates if the functional is evaluated fully implicit or semi implicit.
    bool IsSemImplicitFunctional() { return semImplicitFunctional_; }

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
    //! @name general framework variables
    /*========================================================================*/
    //! phi at time 0
    Teuchos::RCP<Epetra_Vector> phi0_;

   private:
    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for potentially skipping computation of initial time derivative
    bool semImplicitFunctional_;

  };  // class TimIntVarChemDiffusionOST
}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_VAR_CHEMDIFFUSIONOST_H
