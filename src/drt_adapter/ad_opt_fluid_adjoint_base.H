/*----------------------------------------------------------------------*/
/*! \file

\brief base adapter of adjoint fluid equations for topology optimization

\level 3

\maintainer Martin Kronbichler

 *------------------------------------------------------------------------------------------------*/


#ifndef AD_OPT_FLUID_ADJOINT_BASE_H_
#define AD_OPT_FLUID_ADJOINT_BASE_H_


#include <Teuchos_RCP.hpp>
#include "../drt_opti/topopt_fluidAdjoint_timeint.H"


namespace Teuchos
{
  class ParameterList;
}

namespace TOPOPT
{
  namespace ADJOINT
  {
    class FluidAdjointTimeInt;
  }
}  // namespace TOPOPT
namespace ADAPTER
{
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
    explicit TopOptFluidAdjointAlgorithm(const Teuchos::ParameterList& prbdyn);

    /// virtual desctructor
    virtual ~TopOptFluidAdjointAlgorithm();

    virtual void ReadRestart(int step);

    /// fluid adjoint field solver
    Teuchos::RCP<TOPOPT::ADJOINT::FluidAdjointTimeInt> AdjointFluidField()
    {
      return adjointTimeInt_;
    }
    /// return instance of fluid adjoint time integration class
    const Teuchos::RCP<const TOPOPT::ADJOINT::FluidAdjointTimeInt> AdjointFluidField() const
    {
      return adjointTimeInt_;
    }

   protected:
    /// adjoint fluid field
    Teuchos::RCP<TOPOPT::ADJOINT::FluidAdjointTimeInt> adjointTimeInt_;

   private:
    /// setup the adjoint fluid time integration data
    void SetupAdjointFluid(const Teuchos::ParameterList& prbdyn);
  };
}  // namespace ADAPTER

#endif /* ADAPTER_TOPOPT_FLUID_ADJOINT_H_ */
