// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_FLUID_FSI_MSHT_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_FSI_MSHT_HPP


#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  namespace Utils
  {
    class FsiMapExtractor;
  }
}  // namespace FLD

namespace Adapter
{
  /*! \brief Fluid field adapter for fsi with internal mesh tying or mesh sliding
   *
   *
   *  Can only be used in conjunction with #FLD::FluidImplicitTimeInt
   */
  class FluidFSIMsht : public FluidFSI
  {
   public:
    /// Constructor
    FluidFSIMsht(std::shared_ptr<Fluid> fluid, std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    /// initialize algorithm
    void init() override;

    /// communication object at the interface
    virtual std::shared_ptr<FLD::Utils::FsiMapExtractor> const& fsi_interface() const
    {
      return fsiinterface_;
    }


   protected:
    /// create conditioned dof-map extractor for the fluid
    virtual void setup_fsi_interface();

    //! \brief interface map setup for fsi interface and other
    //!
    //! Note: full map contains velocity AND pressure DOFs
    std::shared_ptr<FLD::Utils::FsiMapExtractor> fsiinterface_;
  };
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
