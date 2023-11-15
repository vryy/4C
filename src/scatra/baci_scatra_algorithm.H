/*----------------------------------------------------------------------*/
/*! \file

\brief Transport of passive scalars in Navier-Stokes velocity field

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ALGORITHM_H
#define BACI_SCATRA_ALGORITHM_H

#include "baci_adapter_scatra_fluid_coupling_algorithm.H"

namespace SCATRA
{
  /// Basic algorithm for the transport of passive scalars in the fluid velocity field
  /*!

    \author gjb
    \date 07/08
   */
  class ScaTraAlgorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
   public:
    /// constructor
    ScaTraAlgorithm(const Epetra_Comm& comm,        ///< communicator
        const Teuchos::ParameterList& scatradyn,    ///< scatra parameter list
        const Teuchos::ParameterList& fdyn,         ///< fluid parameter list
        const std::string scatra_disname,           ///< scatra discretization name
        const Teuchos::ParameterList& solverparams  ///< solver parameter list
    );


    void Setup() override;

    void Init() override;

    /// outer level time loop
    void TimeLoop() override;

    /// read restart for preceding turbulent inflow simulation
    void ReadInflowRestart(int restart);

    /// Add tests to global problem and start tests
    void TestResults();

   protected:
    /// time loop for one-way coupling
    void TimeLoopOneWay();

    /// time loop for two-way coupling (natural convection)
    void TimeLoopTwoWay();

    /// provide information about initial field
    virtual void PrepareTimeLoop(){};

    /// initial calculations for two-way coupling time loop
    void PrepareTimeLoopTwoWay();

    /// start a new time step
    void PrepareTimeStep() override;

    /// start a new time step
    void PrepareTimeStepConvection();

    /// print scatra solver type to screen
    virtual void PrintScaTraSolver();

    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();

    /// solve transport (convection-diffusion) equations for current time step
    void DoTransportStep();

    /// set fluid velocites in scatra field
    void SetVelocityField();

    /// Outer iteration loop for natural convection
    void OuterIterationConvection();

    /// take current results for converged and save for next time step
    void Update() override;

    /// take current results for converged and save for next time step
    void UpdateConvection();

    /// write output
    void Output() override;

    /// convergence check for natural convection solver
    virtual bool ConvergenceCheck(int natconvitnum, int natconvitmax, double natconvittol);

    /// flag for natural convection effects
    int natconv_;

    /// maximum iteration steps of outer loop for natural convection
    const int natconvitmax_;

    /// convection tolerance for outer loop for natural convection
    const double natconvittol_;

    /// outer loop velocity increment for natural convection
    Teuchos::RCP<Epetra_Vector> velincnp_;

    /// outer loop phi increment for natural convection
    Teuchos::RCP<Epetra_Vector> phiincnp_;

    /// start step for sampling of statistical data
    const int samstart_;

    /// end step for sampling of statistical data
    const int samstop_;
  };

}  // namespace SCATRA

#endif  // SCATRA_ALGORITHM_H
