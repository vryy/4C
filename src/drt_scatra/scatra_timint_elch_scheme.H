/*----------------------------------------------------------------------*/
/*! \file

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)
\level 2



*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_ELCH_SCHEME_H
#define SCATRA_TIMINT_ELCH_SCHEME_H

#include "scatra_timint_ost.H"
#include "scatra_timint_bdf2.H"
#include "scatra_timint_genalpha.H"
#include "scatra_timint_stat.H"
#include "scatra_timint_elch.H"

namespace SCATRA
{
  class ScaTraTimIntElchOST : public ScaTraTimIntElch, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchOST(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! initialize time integration scheme
    void Init() override;

    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   protected:
    //! calculate initial electric potential field
    void CalcInitialPotentialField() override;

    //! write additional data required for restart
    void OutputRestart() const override;

    //! update time-dependent electrode state variables at the end of an time step
    void ElectrodeKineticsTimeUpdate() override;

    //! explicit predictor for nonlinear solver
    void ExplicitPredictor() const override;

    //!  compute scalar time derivate parameters of the input voltage to compute double layer
    //!  current densities
    void ComputeTimeDerivPot0(const bool init) override;

    //! set part of residual vector belonging to previous time step
    void SetOldPartOfRighthandside() override;
  };  // class TimIntElchOST


  class ScaTraTimIntElchBDF2 : public ScaTraTimIntElch, public TimIntBDF2
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchBDF2(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! initialize time integration scheme
    void Init() override;

    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   protected:
    //! calculate initial electric potential field
    void CalcInitialPotentialField() override;

    //! write additional data required for restart
    void OutputRestart() const override;

    //! update time-dependent electrode state variables at the end of an time step
    void ElectrodeKineticsTimeUpdate() override;

    //!  compute scalar time derivate parameters of the input voltage to compute double layer
    //!  current densities
    void ComputeTimeDerivPot0(const bool init) override;

    //! set part of residual vector belonging to previous time step
    void SetOldPartOfRighthandside() override;
  };  // class TimIntElchBDF2


  class ScaTraTimIntElchGenAlpha : public ScaTraTimIntElch, public TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! initialize time integration scheme
    void Init() override;

    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   protected:
    //! calculate initial electric potential field
    void CalcInitialPotentialField() override;

    //! write additional data required for restart
    void OutputRestart() const override;

    //! update time-dependent electrode state variables at the end of an time step
    void ElectrodeKineticsTimeUpdate() override;

    //!  compute scalar time derivate parameters of the input voltage to compute double layer
    //!  current densities
    void ComputeTimeDerivPot0(const bool init) override;
  };  // class TimIntElchGenAlpha


  class ScaTraTimIntElchStationary : public ScaTraTimIntElch, public TimIntStationary
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! initialize time integration scheme
    void Init() override;

    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   protected:
    //! calculate initial electric potential field
    void CalcInitialPotentialField() override;

    //! write additional data required for restart
    void OutputRestart() const override;

    //! update time-dependent electrode state variables at the end of an time step
    void ElectrodeKineticsTimeUpdate() override
    {
      dserror("Galvanostatic-BC is not implemented for the stationary time-integration scheme");
    };

    //!  compute scalar time derivate parameters of the input voltage to compute double layer
    //!  current densities
    void ComputeTimeDerivPot0(const bool init) override;
  };  // class TimIntElchStatationary
}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_ELCH_SCHEME_H
