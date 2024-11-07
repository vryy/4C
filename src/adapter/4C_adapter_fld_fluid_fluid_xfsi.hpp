// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    FluidFluidXFSI(std::shared_ptr<Fluid> fluid,
        const std::string coupling_name_xfsi,  // name of the FSI coupling condition
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    /// initialize algorithm
    void init() override;

   protected:
    /// A casted pointer to a fluid with multiple discretizations
    std::shared_ptr<FLD::XFluidFluid> xfluidfluid_;
  };
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
