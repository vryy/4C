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

    /// time loop
    void TimeLoop();

  protected:

    /// initial calculations
    void InitialCalculations();

    /// prepare time step
    void PrepareTimeStep();

    /// outer iteration loop
    void OuterLoop();

    /// update for next time step
    void Update();

    /// write output
    void Output();

    /// (preliminary) maximum number of iterations and tolerance for outer iteration
    int     itmax_;
    int     itmaxpre_;
    double  ittol_;

    /// flag for constant thermodynamic pressure
    string consthermpress_;

    /// flag for special flow
    string special_flow_;

    /// start of sampling period
    int samstart_;

  private:

  };

} // namespace LOMA

#endif // LOMA_ALGORITHM_H
#endif // CCADISCRET
