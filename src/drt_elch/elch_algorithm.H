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

    /// get initial fluid density
    double GetInitialFluidDensity();

  private:

    /// flag for natural convection effects
    const string natconv_;

    /// maximum iteration steps
    const int itmax_;

    /// convection tolerance
    const double ittol_;

    /// velocity increment of the outer loop
    Teuchos::RCP<Epetra_Vector> velincnp_;

    /// concentration and potential increment of the outer loop
    Teuchos::RCP<Epetra_Vector> conpotincnp_;

    /// initial fluid density
    const double density0_;

  };

} // namespace ELCH

#endif // ELCH_ALGORITHM_H
#endif // CCADISCRET
