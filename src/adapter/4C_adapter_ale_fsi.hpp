// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_ALE_FSI_HPP
#define FOUR_C_ADAPTER_ALE_FSI_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
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
  /*! \brief ALE Wrapper for FSI Problems
   *
   *  Provide FSI specific ALE functionalities here by overloading the respective
   *  routines from Adapter::AleWrapper
   *
   *  \sa Adapter::Ale, Adapter::AleWrapper
   *
   *  \author mayr.mt \date 10/2014
   */
  class AleFsiWrapper : public AleWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFsiWrapper(std::shared_ptr<Ale> ale);

    //@}

    //! communicate object at the interface
    std::shared_ptr<const ALE::Utils::MapExtractor> interface() const;

    //! apply interface displacements
    void apply_interface_displacements(std::shared_ptr<const Core::LinAlg::Vector<double>> idisp)
    {
      interface_->insert_fsi_cond_vector(*idisp, *write_access_dispnp());
    }

    //! get Dirichlet map extractor
    std::shared_ptr<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() override
    {
      return AleWrapper::get_dbc_map_extractor();
    }

   private:
    std::shared_ptr<ALE::Utils::MapExtractor> interface_;

  };  // class AleFsiWrapper
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
