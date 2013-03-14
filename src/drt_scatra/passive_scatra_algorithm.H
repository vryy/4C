/*----------------------------------------------------------------------*/
/*!
\file passive_scatra_algorithm.H

\brief Transport of passive scalars in Navier-Stokes velocity field
       (one-way coupling)

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef PASSIVE_SCATRA_ALGORITHM_H
#define PASSIVE_SCATRA_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"

namespace SCATRA{

  /// Basic algorithm for the transport of passive scalars in the fluid velocity field
  /*!

    \author gjb
    \date 07/08
   */
  class PassiveScaTraAlgorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor
    explicit PassiveScaTraAlgorithm(
        const Epetra_Comm& comm,              ///< communicator
        const Teuchos::ParameterList& prbdyn, ///< problem-specific parameter list
        const std::string disname,            ///< scatra discretization name
        const Teuchos::ParameterList& solverparams
        );

    /// virtual destructor to support polymorph destruction
    virtual ~PassiveScaTraAlgorithm();

    /// outer level time loop
    void TimeLoop();

    /// start a new time step
    void PrepareTimeStep();

    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();

    /// solve transport (convection-diffusion) equations for current time step
    void DoTransportStep();

    /// take current results for converged and save for next time step
    void Update(const int num=0);

    /// write output
    void Output();
    
    /// read restart for preceding turbulent inflow simulation
    void ReadInflowRestart(int restart);

  private:

    /// start step for sampling of statistical data
    const int samstart_;

    /// end step for sampling of statistical data
    const int samstop_;

  };

} // namespace SCATRA

#endif // PASSIVE_SCATRA_ALGORITHM_H
