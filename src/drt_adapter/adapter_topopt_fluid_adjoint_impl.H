/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_fluid_adjoint_impl.H

\brief adapter for element routines of fluid adjoint equations in topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_
#define ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_


#include "adapter_topopt_fluid_adjoint.H"


/// forward declarations
namespace DRT
{
  class Discretization;
  class ResultTest;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class Solver;
}

namespace TOPOPT
{
  class Optimizer;
  namespace ADJOINT
  {
    class ImplicitTimeInt;
  }
}

// Algorithm adapter classes for use in topology optimization of fluid fields
namespace ADAPTER {

  /// adapter to adjoint fluid algorithm
  /*!
   *
   * This adapter communicates the element level of the fluid adjoint equations
   * when an implicit time integration scheme is used.

    \author winklmaier
    \date 01/12
   */
  class FluidAdjointImpl : public FluidAdjoint
  {
  public:

    /// constructor
    FluidAdjointImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! @name Vector access

    /// velocities (and pressures) at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();

    /// velocities (and pressures) at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Veln();

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    /// set initial flow field
    virtual void SetInitialAdjointField(const INPAR::TOPOPT::InitialAdjointField initfield,const int startfuncno);

    //@}

    virtual void TimeLoop();
    virtual void Output();
    virtual IO::DiscretizationWriter& DiscWriter();

    virtual void ReadRestart(int step);

    /*!
    \brief return type of time integration scheme

    */
    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const;

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

    virtual void SetTopOptData(
        Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > fluidvelocities,
        Teuchos::RCP<Epetra_Vector> porosity,
        Teuchos::RCP<TOPOPT::Optimizer> optimizer);


  private:

    /// the actual fluid algorithm
    Teuchos::RCP<TOPOPT::ADJOINT::ImplicitTimeInt> adjointTimeInt_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}
  };

}

#endif /* ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_ */
