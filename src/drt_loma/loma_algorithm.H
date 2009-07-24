/*----------------------------------------------------------------------*/
/*!
\file loma_algorithm.H

\brief Basis of all LOMA algorithms

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef LOMA_ALGORITHM_H
#define LOMA_ALGORITHM_H

#include <Epetra_Comm.h>
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"


/*!
\brief LOMA: namespace of the module for low-Mach-number flow problems in BACI

*/
namespace LOMA
{

  /// LOMA algorithm base
  /*!

    Base class of LOMA algorithms. Derives from FluidBaseAlgorithm
    and ScatraBaseAlgorithm.
    There can (and will) be different subclasses that implement
    different coupling schemes.

    \author vg
    \date 08/08
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

    /// just the general entry to LOMA algorithm where particular algorithm
    /// will be decided on
    void TimeLoop();

  protected:

    /// LOMA stationary algorithm
    void SolveStationaryProblem();

    /// LOMA generalized-alpha time loop
    void GenAlphaTimeLoop();

    /// LOMA one-step-theta/BDF2 time loop
    void OSTBDF2TimeLoop();

    /// prepare time step
    void InitialCalculations();

    /// prepare generalized-alpha time step
    void GenAlphaPrepareTimeStep();

    /// prepare one-step-theta/BDF2 time step
    void OSTBDF2PrepareTimeStep();

    /// do generalized-alpha outer iteration loop
    void GenAlphaOuterLoop();

    /// do one-step-theta/BDF2 outer iteration loop
    void OSTBDF2OuterLoop();

    /// do generalized-alpha update for next time step
    void GenAlphaUpdate();

    /// do one-step-theta/BDF2 update for next time step
    void OSTBDF2Update();

    /// write output
    void Output();

    /// (preliminary) maximum number of iterations and tolerance for outer iteration
    int     itmax_;
    int     itmaxpre_;
    double  ittol_;

    /// flag for constant thermodynamic pressure
    string consthermpress_;

    /// thermodynamic pressure (default: 98100.0)
    double thermpress_;

    /// gas constant (default: 287.0)
    double gasconstant_;

    /// factor for equation of state (i.e., therm. press. / gas constant)
    double eosfac_;

    /// initial total mass (required if thermodynamic pressure is
    /// calculated from mass conservation)
    double initialmass_;

    /// flag for special flow
    string special_flow_;

    /// start of sampling period
    int samstart_;

  private:

  };

} // namespace LOMA

#endif // LOMA_ALGORITHM_H
#endif // CCADISCRET
