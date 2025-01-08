// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_UTILS_ELCH_DIFFCOND_HPP
#define FOUR_C_SCATRA_ELE_UTILS_ELCH_DIFFCOND_HPP

#include "4C_config.hpp"

#include "4C_inpar_elch.hpp"
#include "4C_scatra_ele_utils_elch_electrode.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    // forward declarations
    template <int nsd, int nen>
    class ScaTraEleInternalVariableManagerElchDiffCond;
    class ScaTraEleDiffManagerElchDiffCond;

    // class implementation
    template <Core::FE::CellType distype>
    class ScaTraEleUtilsElchDiffCond : public ScaTraEleUtilsElchElectrode<distype>
    {
      //! abbreviations
      typedef ScaTraEleUtilsElch<distype> myelch;
      typedef ScaTraEleUtilsElchElectrode<distype> myelectrode;

     public:
      //! singleton access method
      static ScaTraEleUtilsElchDiffCond<distype>* instance(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

      //! evaluate electrolyte material
      void mat_elch_mat(
          std::shared_ptr<const Core::Mat::Material> material,  //!< electrolyte material
          const std::vector<double>& concentrations,            //!< local concentration values
          double temperature,                                   //!< temperature
          Inpar::ElCh::EquPot equpot,  //!< type of closing equation for electric potential
          double ffrt,                 //!< factor F^2/RT
          std::shared_ptr<ScaTraEleDiffManagerElchDiffCond> diffmanager,  //!< diffusion manager
          Inpar::ElCh::DiffCondMat& diffcondmat                           //!< ion type
      );

      //! evaluate electrolyte phase
      void mat_elch_phase(
          std::shared_ptr<const Core::Mat::Material> material,  //!< electrolyte phase
          const std::vector<double>& concentrations,            //!< local concentration values
          double temperature,                                   //!< temperature
          const Inpar::ElCh::EquPot& equpot,  //!< type of closing equation for electric potential
          const double& ffrt,                 //!< factor F^2/RT
          std::shared_ptr<ScaTraEleDiffManagerElchDiffCond> diffmanager,  //!< diffusion manager
          Inpar::ElCh::DiffCondMat& diffcondmat                           //!< ion type
      );

      //! evaluate standard Newman material
      void mat_newman(std::shared_ptr<const Core::Mat::Material> material,  //!< Newman material
          double concentration,  //!< local concentration value
          double temperature,    //!< temperature
          std::shared_ptr<ScaTraEleDiffManagerElchDiffCond> diffmanager  //!< diffusion manager
      );

     protected:
      //! private constructor for singletons
      ScaTraEleUtilsElchDiffCond(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

    };  // class ScaTraEleUtilsElchDiffCond
  }  // namespace Elements
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
