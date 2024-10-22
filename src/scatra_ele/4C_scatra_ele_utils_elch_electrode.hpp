// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_UTILS_ELCH_ELECTRODE_HPP
#define FOUR_C_SCATRA_ELE_UTILS_ELCH_ELECTRODE_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_utils_elch.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerElchElectrode;

    // class implementation
    template <Core::FE::CellType distype>
    class ScaTraEleUtilsElchElectrode : public ScaTraEleUtilsElch<distype>
    {
      //! abbreviation
      typedef ScaTraEleUtilsElch<distype> myelch;

     public:
      //! singleton access method
      static ScaTraEleUtilsElchElectrode<distype>* instance(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );



      //! evaluate electrode material
      void mat_electrode(Teuchos::RCP<const Core::Mat::Material> material,  //!< electrode material
          double concentration,                                             //!< concentration
          double temperature,                                               //!< temperature
          Teuchos::RCP<ScaTraEleDiffManagerElchElectrode> diffmanager       //!< diffusion manager
      );

     protected:
      //! protected constructor for singletons
      ScaTraEleUtilsElchElectrode(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );
    };  // class ScaTraEleUtilsElchElectrode
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
