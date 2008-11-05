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

#ifdef CCADISCRET

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
        Epetra_Comm& comm,                   ///< communicator
        const Teuchos::ParameterList& prbdyn ///< problem-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~PassiveScaTraAlgorithm();

    /// outer level time loop
    void TimeLoop();

  protected:

    /// start a new time step
    void PrepareTimeStep();
    
    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();
    
    /// solve transport (convection-diffusion) equations for current time step
    void DoTransportStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

  private:

  };

} // namespace SCATRA

#endif // PASSIVE_SCATRA_ALGORITHM_H
#endif // CCADISCRET
