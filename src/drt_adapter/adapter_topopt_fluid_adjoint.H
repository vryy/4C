/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_fluid_adjoint.H

\brief Fluid adjoint field adapter

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_TOPOPT_FLUID_ADJOINT_H_
#define ADAPTER_TOPOPT_FLUID_ADJOINT_H_


#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

namespace LINALG
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class Solver;
}

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

namespace FLD
{
  class TurbulenceStatisticManager;
  namespace UTILS
  {
    class MapExtractor;
  }
}

namespace COMBUST
{
  class InterfaceHandleCombust;
  class FlameFront;
}

namespace TOPOPT
{
  class Optimizer;
}

/// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general fluid adjoint field interface for topology optimization
  /*!

  This is the general interface for the fluid adjoints equations. Different
  approaches exist for the computation of the adjoints of the fluid in the
  context of topology optimization (f.e. discrete and continous)

  This pure virtual class shall define an interface of the required functions.

  \author winklmaier
  \date 01/12
  */
  class FluidAdjoint
  {
  public:

    /// virtual destructor to get polymorph destruction
    virtual ~FluidAdjoint(){};

    //! @name Vector access

    /// velocities (and pressures) at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;

    /// velocities (and pressures) at t(n+alpha_F)
    virtual Teuchos::RCP<const Epetra_Vector> Velaf() = 0;

    /// velocities (and pressures) at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;

    /// velocities (and pressures) at t(n-1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno) = 0;
    //@}

    //! @name Time step helpers

    /// run a complete simulation (for fluid stand-alone simulations)
    virtual void TimeLoop() = 0;

    /// update fluid unknowns and evaluate elements
    ///
    /// there are two increments possible
    /// x^n+1_i+1 = x^n+1_i + iterinc, and
    ///
    /// x^n+1_i+1 = x^n     + stepinc
    ///
    /// with n and i being time and Newton iteration step
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc ///< increment between time step n and n+1
        ) = 0;

    /// output results
    virtual void Output() = 0;

    /// access to output
    virtual IO::DiscretizationWriter& DiscWriter() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;
    /// set restart
    virtual void SetRestart(const int step,
                    const double time,
                    Teuchos::RCP<const Epetra_Vector> readvelnp,
                    Teuchos::RCP<const Epetra_Vector> readveln,
                    Teuchos::RCP<const Epetra_Vector> readvelnm,
                    Teuchos::RCP<const Epetra_Vector> readaccnp,
                    Teuchos::RCP<const Epetra_Vector> readaccn) = 0;
    //@}

    /// return type of time integration scheme
    virtual const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const = 0;

    //! @name Transfer the porosity field used in the topology optimization to the fluid

    /// set pseudo-porosity field in the fluid
    virtual void SetTopOptData(RCP<Epetra_Vector> porosity,RCP<TOPOPT::Optimizer> optimizer) = 0;
    //@}

    /// create result test for encapulated fluid algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };
}

#endif /* ADAPTER_TOPOPT_FLUID_ADJOINT_H_ */
#endif  // #ifdef CCADISCRET
