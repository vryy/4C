/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for AC-FS3I problems



\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_AC_FSI_H
#define ADAPTER_FLUID_AC_FSI_H

#include "ad_fld_fluid_fsi.H"

namespace ADAPTER
{
  class FluidACFSI : public FluidFSI
  {
   public:
    /// Constructor
    FluidACFSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// Get vector of current relative pressures errors from the Windkessels
    std::vector<double> GetWindkesselErrors();

    //@}
  };
}  // namespace ADAPTER

#endif
