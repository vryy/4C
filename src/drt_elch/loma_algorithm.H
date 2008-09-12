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
        RCP<DRT::Discretization>  scatradis,
        const Teuchos::ParameterList& prbdyn ///< problem-specific parameter list
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer LOMA time loop
    void TimeLoop();

  protected:

    /// prepare time step
    void PrepareTimeStep();

    /// prepare time step
    void PrepareIterationStep();

    /// do outer iteration loop
    void OuterLoop();

    /// compute current density based on current temperature + therm. pres.
    void ComputeDensity();

    /// save current (converged) results for next time step
    void Update();

    /// write output
    void Output();

    /// thermodynamic pressure at time step n+1
    double thermpressnp_;

    /// scatra discretization
    RCP<DRT::Discretization>      scatradiscret_;

    /// node-based densities at time steps n+1, n and n-1
    RCP<Epetra_Vector>    noddensnp_;
    RCP<Epetra_Vector>    noddensn_;
    RCP<Epetra_Vector>    noddensnm_;

    /// node-based density increment for convergence check
    RCP<Epetra_Vector>    densincnp_;

    /// node-based temperatures at time step n+1
    RCP<Epetra_Vector>    nodtempnp_;

  private:

  };

} // namespace LOMA

#endif // LOMA_ALGORITHM_H
#endif // CCADISCRET
