/*----------------------------------------------------------------------*/
/*! \file
 \brief One-step-theta time integration scheme for porous multiphase flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROFLUIDMULTIPHASE_TIMINT_OST_H
#define BACI_POROFLUIDMULTIPHASE_TIMINT_OST_H

#include "baci_linalg_utils_densematrix_communication.H"
#include "baci_porofluidmultiphase_timint_implicit.H"

namespace POROFLUIDMULTIPHASE
{
  class TimIntOneStepTheta : public TimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntOneStepTheta(Teuchos::RCP<DRT::Discretization> dis,  //!< discretization
        const int linsolvernumber,                             //!< number of linear solver
        const Teuchos::ParameterList& probparams,              //!< parameter list of global problem
        const Teuchos::ParameterList& poroparams,              //!< paramter list of poro problem
        Teuchos::RCP<IO::DiscretizationWriter> output          //!< output writer
    );


    /// Print information about current time step to screen (reimplementation for OST)
    void PrintTimeStepInfo() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() override { return; };

    ///  compute scalar time derivative
    void ComputeTimeDerivative() override;

    /// update the solution after convergence of the nonlinear iteration.
    /// current solution becomes old solution of next timestep.
    void Update() override;

    /// read restart data
    void ReadRestart(const int step) override;

   protected:
    /// don't want = operator and cctor
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old);

    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    /// set time parameter for element evaluation (called before every time step)
    void SetElementTimeStepParameter() const override;

    //! set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    //! initialization procedure prior to evaluation of first time step
    void CalcInitialTimeDerivative() override;

    /// set part of residual vector belonging to previous time step
    void SetOldPartOfRighthandside() override;

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    virtual void ExplicitPredictor();

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual() override;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors() override;

    /// write additional data required for restart
    void OutputRestart() override;

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const override { return 1.0 / (dt_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

  };  // class TimIntOneStepTheta

}  // namespace POROFLUIDMULTIPHASE



#endif  // POROFLUIDMULTIPHASE_TIMINT_OST_H
