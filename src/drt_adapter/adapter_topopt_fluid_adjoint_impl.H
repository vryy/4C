/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_fluid_adjoint_impl.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_
#define ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_


#include "adapter_topopt_fluid_adjoint.H"


//forward declarations
namespace DRT
{
  class Discretization;
  class DofSet;
  class ResultTest;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MapExtractor;
  class BlockSparseMatrixBase;
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


    \author winklmaier
    \date 01/12
   */
  class FluidAdjointImpl : public FluidAdjoint
  {
  public:
    FluidAdjointImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    virtual Teuchos::RCP<const Epetra_Vector> Velnp();
    virtual Teuchos::RCP<const Epetra_Vector> Veln();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();

    virtual void TimeLoop();
    virtual void Output();
    virtual IO::DiscretizationWriter& DiscWriter();

    virtual void ReadRestart(int step);
    virtual void SetRestart(const int step,
                    const double time,
                    Teuchos::RCP<const Epetra_Vector> readvelnp,
                    Teuchos::RCP<const Epetra_Vector> readveln,
                    Teuchos::RCP<const Epetra_Vector> readvelnm,
                    Teuchos::RCP<const Epetra_Vector> readaccnp,
                    Teuchos::RCP<const Epetra_Vector> readaccn);

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

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno);

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
#endif  // #ifdef CCADISCRET
