/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for AC-FS3I problems



\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_FLD_FLUID_AC_FSI_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_AC_FSI_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FluidACFSI : public FluidFSI
  {
   public:
    /// Constructor
    FluidACFSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Discret::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// Get vector of current relative pressures errors from the Windkessels
    std::vector<double> GetWindkesselErrors();

    //@}
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
