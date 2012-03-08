/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_fluid_adjoint_base.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_TOPOPT_FLUID_ADJOINT_BASE_H_
#define ADAPTER_TOPOPT_FLUID_ADJOINT_BASE_H_


//forward declarations
namespace DRT
{
  class Discretization;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

namespace TOPOPT
{
  class Optimizer;
}

/// fluid algorithms
namespace FLD
{
  namespace UTILS
  {
    class MapExtractor;
    class KSPMapExtractor;
    class FluidImpedanceWrapper;
    class FluidWkOptimizationWrapper;
    class FluidVolumetricSurfaceFlowWrapper;
    class TotalTractionCorrector;
  }
}


namespace ADAPTER
{
  class FluidAdjoint;

  /// adapter for the coupling of the fluid adjoint algorithm in topology optimization
  /*!

    couples either the discrete or continuous adjoint approach to the other parts
    of a topology optimization due to a density approach

    \author winklmaier
    \date 01/12
   */
  class TopOptFluidAdjointAlgorithm
  {
  public:

    /// constructor
    explicit TopOptFluidAdjointAlgorithm(
        const Teuchos::ParameterList& prbdyn
        );

    /// virtual desctructor
    virtual ~TopOptFluidAdjointAlgorithm();

    virtual void ReadRestart(int step);

    /// fluid adjoint field solver
    Teuchos::RCP<FluidAdjoint> AdjointFluidField();
    /// return instance of fluid adjoint time integration class
    const Teuchos::RCP<const FluidAdjoint> AdjointFluidField() const;

  protected:

    /// adjoint fluid field
    Teuchos::RCP<FluidAdjoint> adjoint_;

  private:

    void SetupAdjointFluid(const Teuchos::ParameterList& prbdyn);
  };
} // namespace ADAPTER

#endif /* ADAPTER_TOPOPT_FLUID_ADJOINT_H_ */
#endif  // #ifdef CCADISCRET
