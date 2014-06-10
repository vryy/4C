/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.H

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ELCH_ALGORITHM_H
#define ELCH_ALGORITHM_H

#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
//#include "../drt_inpar/inpar_elch.H"


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
        const Epetra_Comm& comm,                   ///< communicator
        const Teuchos::ParameterList& prbdyn, ///< problem-specific parameter list
        const Teuchos::ParameterList& solverparams
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level ELCH time loop
    void TimeLoop();

  protected:

    /// time loop for elch algorithm without natural convection
    void TimeLoopElch();

    /// time loop for elch algorithm with natural convection
    void TimeLoopConvection();

    /// initial calculations
    void InitialCalculations();

    /// start a new time step
    void PrepareTimeStep();

    /// start a new time step
    void PrepareTimeStepConvection();

    /// solve Navier-Stokes equations for current time step
    void DoFluidStep();

    /// solve transport equations for current time step
    void DoTransportStep();

    /// take current results for converged and save for next time step
    void Update();

    /// take current results for converged and save for next time step
    void UpdateConvection();

    /// write output
    void Output();

    /// Outer Iteration loop convection
    void OuterIterationConvection();

    /// convergence check for natural convection solver
    bool ConvergenceCheck(int itnum, int itmax, double ittol);

  private:

    /// flag for natural convection effects
    int natconv_;

    /// maximum iteration steps
    const int itmax_;

    /// convection tolerance
    const double ittol_;

    /// velocity increment of the outer loop
    Teuchos::RCP<Epetra_Vector> velincnp_;

    /// concentration and potential increment of the outer loop
    Teuchos::RCP<Epetra_Vector> conpotincnp_;

    /// start step for sampling of statistical data
    const int samstart_;

    /// end step for sampling of statistical data
    const int samstop_;

  };

} // namespace ELCH

#endif // ELCH_ALGORITHM_H
