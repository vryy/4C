/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_timint_ost.H

 \brief One-step-theta time integration scheme for porous multiphase flow

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_TIMINT_OST_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_TIMINT_OST_H_

#include "../linalg/linalg_utils.H"
#include "porofluidmultiphase_timint_implicit.H"

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
        FILE* errfile,                                         //!< error file
        Teuchos::RCP<IO::DiscretizationWriter> output          //!< output writer
    );

    /// Destructor
    virtual ~TimIntOneStepTheta();

    /// Print information about current time step to screen (reimplementation for OST)
    virtual void PrintTimeStepInfo();

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() { return; };

    ///  compute scalar time derivative
    void ComputeTimeDerivative();

    /// update the solution after convergence of the nonlinear iteration.
    /// current solution becomes old solution of next timestep.
    virtual void Update();

    /// read restart data
    void ReadRestart(const int step);

   protected:
    /// don't want = operator and cctor
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old);

    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    /// set time parameter for element evaluation (called before every time step)
    virtual void SetElementTimeStepParameter() const;

    //! set time for evaluation of Neumann boundary conditions
    virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

    //! initialization procedure prior to evaluation of first time step
    virtual void CalcInitialTimeDerivative();

    /// set part of residual vector belonging to previous time step
    virtual void SetOldPartOfRighthandside();

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    virtual void ExplicitPredictor();

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual();

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors();

    /// write additional data required for restart
    void OutputRestart();

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const { return 1.0 / (dt_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

  };  // class TimIntOneStepTheta

}  // namespace POROFLUIDMULTIPHASE



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_TIMINT_OST_H_ */
