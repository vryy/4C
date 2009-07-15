/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.H

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ELCH_ALGORITHM_H
#define ELCH_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"


/*!
\brief ELCH: namespace of the electrochemistry module in BACI

*/
namespace ELCH
{

  /// ELCH algorithm base
  /*!

    Base class of ELCH algorithms. Derives from FluidBaseAlgorithm
    and ConDifBaseAlgorithm.
    There can (and will) be different subclasses that implement
    different coupling schemes.

    \author gjb
    \date 03/08
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor
    explicit Algorithm(
        Epetra_Comm& comm,                   ///< communicator
        const Teuchos::ParameterList& prbdyn ///< problem-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level ELCH time loop
    void TimeLoop();

  protected:

    /// start a new time step
    void PrepareTimeStep();

    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();

    /// solve transport equations for current time step
    void DoTransportStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

  private:

    /// flag for printing out mean values of concentrations and density
    const bool outmean_;

  };

} // namespace ELCH

#endif // ELCH_ALGORITHM_H
#endif // CCADISCRET
