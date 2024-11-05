// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_UTILS_ELCH_SCL_HPP
#define FOUR_C_SCATRA_ELE_UTILS_ELCH_SCL_HPP
#include "4C_config.hpp"

#include "4C_inpar_elch.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"
#include "4C_scatra_ele_utils_elch_diffcond.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    // forward declarations
    template <int nsd, int nen>
    class ScaTraEleInternalVariableManagerScl;
    class ScaTraEleDiffManagerElchScl;

    // class implementation
    template <Core::FE::CellType distype>
    class ScaTraEleUtilsElchScl : public ScaTraEleUtilsElchDiffCond<distype>
    {
      //! abbreviations
      using myelch = ScaTraEleUtilsElch<distype>;
      using mydiffcond = ScaTraEleUtilsElchDiffCond<distype>;

     public:
      //! singleton access method
      static ScaTraEleUtilsElchScl<distype>* instance(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

      //! evaluate electrolyte material
      void mat_elch_mat(
          std::shared_ptr<const Core::Mat::Material> material,       //!< electrolyte material
          const std::vector<double>& concentrations,                 //!< local concentration values
          double temperature,                                        //!< temperature
          std::shared_ptr<ScaTraEleDiffManagerElchScl> diffmanager,  //!< diffusion manager
          Inpar::ElCh::DiffCondMat& diffcondmat                      //!< ion type
      );
      //! evaluate electrolyte phase
      void mat_elch_phase(
          std::shared_ptr<const Core::Mat::Material> material,       //!< electrolyte phase
          const std::vector<double>& concentrations,                 //!< local concentration values
          double temperature,                                        //!< temperature
          std::shared_ptr<ScaTraEleDiffManagerElchScl> diffmanager,  //!< diffusion manager
          Inpar::ElCh::DiffCondMat& diffcondmat                      //!< ion type
      );

      //! evaluate Space Charge Layer Material
      void mat_scl(std::shared_ptr<const Core::Mat::Material> material, const double concentration,
          const double temperature, std::shared_ptr<ScaTraEleDiffManagerElchScl> diffmanager);

     protected:
      //! private constructor for singletons
      ScaTraEleUtilsElchScl(const int numdofpernode,  ///< number of degrees of freedom per node
          const int numscal,                          ///< number of transported scalars per node
          const std::string& disname                  ///< name of discretization
      );
    };
  }  // namespace Elements
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
