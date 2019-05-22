/*!----------------------------------------------------------------------
\brief One-Step-Theta time-integration scheme

\level 1

\maintainer Anh-Tu Vuong

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_OST_H
#define SCATRA_TIMINT_OST_H

#include "scatra_timint_implicit.H"

#include "../linalg/linalg_utils.H"

namespace SCATRA
{
  class TimIntOneStepTheta : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntOneStepTheta(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<LINALG::Solver> solver,                      //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,              //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,         //!< supplementary parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output,            //!< output writer
        const int probnum = 0                                     //!< global problem number
    );

    /// Destructor
    virtual ~TimIntOneStepTheta();

    /// initialize time integration scheme
    virtual void Init();

    /// setup time integration scheme
    virtual void Setup();

    //! operator for manipulations before call to \ref Solve() ; May be overridden by subclass.
    virtual void PreSolve() { return; };

    //! operator for manipulations after call to \ref Solve() ; May be overridden by subclass.
    virtual void PostSolve() { return; };

    /// Print information about current time step to screen (reimplementation for OST)
    virtual void PrintTimeStepInfo();

    /// compute values at intermediate time steps (required for generalized-alpha)
    virtual void ComputeIntermediateValues() { return; };

    /// compute values at the interior of the elements (required for hdg)
    virtual void ComputeInteriorValues() { return; };

    ///  compute scalar time derivative
    virtual void ComputeTimeDerivative();

    ///  compute scalar time derivate parameters of the input voltage
    virtual void ComputeTimeDerivPot0(const bool init) { return; };

    /// update the solution after convergence of the nonlinear iteration.
    /// current solution becomes old solution of next timestep.
    virtual void Update(const int num = 0);

    /// read restart data
    virtual void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

    // routine to return scalar field phi at time step n-1
    virtual Teuchos::RCP<Epetra_Vector> Phinm() { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_F
    virtual Teuchos::RCP<Epetra_Vector> Phiaf() { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_M
    virtual Teuchos::RCP<Epetra_Vector> Phiam() { return Teuchos::null; }

    /// routine to return time derivative of scalar field phi at time step n+alpha_M
    virtual Teuchos::RCP<Epetra_Vector> Phidtam() { return Teuchos::null; }

    /// routine to return fine-scale scalar field fsphi at time step n+1
    virtual Teuchos::RCP<Epetra_Vector> FsPhi()
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    };

    /// routine to return time integration specific parameters
    virtual Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
    {
      Teuchos::RCP<Teuchos::ParameterList> timeparams;
      timeparams = Teuchos::rcp(new Teuchos::ParameterList());
      timeparams->set("using stationary formulation", false);
      timeparams->set("using generalized-alpha time integration", false);
      timeparams->set("total time", time_);
      timeparams->set("time factor", theta_ * dta_);
      timeparams->set("alpha_F", 1.0);
      return timeparams;
    }

    //! set state on micro scale in multi-scale simulations
    void SetState(Teuchos::RCP<Epetra_Vector> phin,  //!< micro-scale state vector at old time step
        Teuchos::RCP<Epetra_Vector> phinp,           //!< micro-scale state vector at new time step
        Teuchos::RCP<Epetra_Vector>
            phidtn,  //!< time derivative of micro-scale state vector at old time step
        Teuchos::RCP<Epetra_Vector>
            phidtnp,  //!< time derivative of micro-scale state vector at new time step
        Teuchos::RCP<Epetra_Vector> hist,               //!< micro-scale history vector
        Teuchos::RCP<IO::DiscretizationWriter> output,  //!< micro-scale discretization writer
        const std::vector<double>&
            phinp_macro,   //!< values of state variables at macro-scale Gauss point
        const int step,    //!< time step
        const double time  //!< time
    );

    //! clear state on micro scale in multi-scale simulations
    void ClearState();

   protected:
    /// don't want = operator and cctor
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old);

    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    /// set time parameter for element evaluation
    virtual void SetElementTimeParameter(bool forcedincrementalsolver = false) const;

    //! set time for evaluation of Neumann boundary conditions
    virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

    //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
    virtual void CalcInitialTimeDerivative();

    /// set part of residual vector belonging to previous time step
    virtual void SetOldPartOfRighthandside();

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    virtual void ExplicitPredictor() const;

    /// add actual Neumann loads with time factor
    virtual void AddNeumannToResidual();

    /// AVM3-based scale separation
    virtual void AVM3Separation();

    /// dynamic Smagorinsky model
    virtual void DynamicComputationOfCs();

    // dynamic Vreman model
    virtual void DynamicComputationOfCv();

    /// add parameters specific for time-integration scheme
    virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);

    /// write additional data required for restart
    virtual void OutputRestart() const;

    /// return the right time-scaling-factor for the true residual
    virtual double ResidualScaling() const { return 1.0 / (dta_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;

  };  // class TimIntOneStepTheta

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_OST_H
