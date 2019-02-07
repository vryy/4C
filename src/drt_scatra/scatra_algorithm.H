/*----------------------------------------------------------------------*/
/*!
\file scatra_algorithm.H

\brief Transport of passive scalars in Navier-Stokes velocity field

\level 1

\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ALGORITHM_H
#define SCATRA_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"

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

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraAlgorithm();

    /// setup
    virtual void Setup();

    /// init
    virtual void Init(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            scatradyn,  ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
        const std::string& disname = "scatra",       ///< name of scalar transport discretization
        const bool isale = false                     ///< ALE flag
    );

    /// outer level time loop
    void TimeLoop();

    /// read restart for preceding turbulent inflow simulation
    void ReadInflowRestart(int restart);

   protected:
    /// time loop for one-way coupling
    void TimeLoopOneWay();

    /// time loop for two-way coupling (natural convection)
    void TimeLoopTwoWay();

    /// provide information about initial field
    virtual void PrepareTimeLoop()
    {
      return;
    };  // do nothing in base algorithm (only relevant for derived algorithms)

    /// initial calculations for two-way coupling time loop
    void PrepareTimeLoopTwoWay();

    /// start a new time step
    void PrepareTimeStep();

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
    void Update(const int num = 0);

    /// take current results for converged and save for next time step
    void UpdateConvection();

    /// write output
    void Output();

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
