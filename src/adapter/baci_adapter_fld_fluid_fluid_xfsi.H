/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for XFSI allowing multiple fluid discretizations. Can only be used in
conjunction with XFluidFluid!

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_FLD_FLUID_FLUID_XFSI_H
#define BACI_ADAPTER_FLD_FLUID_FLUID_XFSI_H

#include "baci_config.H"

#include "baci_adapter_fld_fluid_xfsi.H"

BACI_NAMESPACE_OPEN

namespace FLD
{
  class XFluidFluid;
}

namespace ADAPTER
{
  class FluidFluidXFSI : public XFluidFSI
  {
   public:
    /// Constructor
    FluidFluidXFSI(Teuchos::RCP<Fluid> fluid,
        const std::string coupling_name_xfsi,  // name of the FSI coupling condition
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// initialize algorithm
    void Init() override;

   protected:
    /// A casted pointer to a fluid with multiple discretizations
    Teuchos::RCP<FLD::XFluidFluid> xfluidfluid_;
  };
}  // namespace ADAPTER


BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_FLD_FLUID_FLUID_XFSI_H
