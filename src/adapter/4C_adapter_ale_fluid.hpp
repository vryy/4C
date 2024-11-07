// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_ALE_FLUID_HPP
#define FOUR_C_ADAPTER_ALE_FLUID_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/* forward declarations */
namespace ALE
{
  namespace Utils
  {
    class MapExtractor;
  }
}  // namespace ALE
/*----------------------------------------------------------------------------*/
/* class definitions */
namespace Adapter
{
  /*! \brief ALE Wrapper for fluid problems with moving boundaries
   *
   *  Provide ALE functionalities for problems involving moving boundaries
   *  here by overloading the respective routines from Adapter::AleWrapper
   *
   *  \sa Adapter::Ale, Adapter::AleWrapper
   *
   *  \author kruse \date 10/2014
   */
  class AleFluidWrapper : public AleWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFluidWrapper(std::shared_ptr<Ale> ale);

    //@}

    //! communicate object at the interface
    std::shared_ptr<const ALE::Utils::MapExtractor> interface() const;

    //! solve
    int solve() override;

    //! apply displacements at the free surface nodes
    void apply_free_surface_displacements(
        std::shared_ptr<const Core::LinAlg::Vector<double>> fsdisp);

    //! apply displacements at the ale update condition nodes
    void apply_ale_update_displacements(std::shared_ptr<const Core::LinAlg::Vector<double>> audisp);

    //! apply FSI interface displacements (required for staggered FSI)
    void apply_interface_displacements(std::shared_ptr<const Core::LinAlg::Vector<double>> idisp);

   private:
    std::shared_ptr<ALE::Utils::MapExtractor> interface_;

  };  // class AleFluidWrapper
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
