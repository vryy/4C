// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_ALE_FPSI_HPP
#define FOUR_C_ADAPTER_ALE_FPSI_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace Adapter
{
  class AleFpsiWrapper : public AleWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFpsiWrapper(std::shared_ptr<Ale> ale);

    //! specialized method to apply displacements to fpsi interface
    void apply_interface_displacements(std::shared_ptr<const Core::LinAlg::Vector<double>> idisp);

    //! specialized method to apply displacements to fsi interface
    void apply_fsi_interface_displacements(
        std::shared_ptr<const Core::LinAlg::Vector<double>> idisp);

    //! communicate object at the interface
    std::shared_ptr<const ALE::Utils::MapExtractor> interface() const;

    //@}

   private:
    //! interface map extractor
    std::shared_ptr<ALE::Utils::MapExtractor> interface_;


  };  // class AleFpsiWrapper
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
