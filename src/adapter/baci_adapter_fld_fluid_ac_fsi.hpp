/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for AC-FS3I problems



\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_ADAPTER_FLD_FLUID_AC_FSI_HPP
#define BACI_ADAPTER_FLD_FLUID_AC_FSI_HPP

#include "baci_config.hpp"

#include "baci_adapter_fld_fluid_fsi.hpp"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class FluidACFSI : public FluidFSI
  {
   public:
    /// Constructor
    FluidACFSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// Get vector of current relative pressures errors from the Windkessels
    std::vector<double> GetWindkesselErrors();

    //@}
  };
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
