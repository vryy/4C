/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_ac_fsi.H

\brief Fluid field adapter for AC-FS3I problems

\date 2015-07-29

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            089/289-10364

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
