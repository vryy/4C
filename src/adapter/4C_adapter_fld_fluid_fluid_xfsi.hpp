/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for XFSI allowing multiple fluid discretizations. Can only be used in
conjunction with XFluidFluid!

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_FLD_FLUID_FLUID_XFSI_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_FLUID_XFSI_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_xfsi.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class XFluidFluid;
}

namespace Adapter
{
  class FluidFluidXFSI : public XFluidFSI
  {
   public:
    /// Constructor
    FluidFluidXFSI(Teuchos::RCP<Fluid> fluid,
        const std::string coupling_name_xfsi,  // name of the FSI coupling condition
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    /// initialize algorithm
    void Init() override;

   protected:
    /// A casted pointer to a fluid with multiple discretizations
    Teuchos::RCP<FLD::XFluidFluid> xfluidfluid_;
  };
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
