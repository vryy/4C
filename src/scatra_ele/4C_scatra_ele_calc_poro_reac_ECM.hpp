// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_PORO_REAC_ECM_HPP
#define FOUR_C_SCATRA_ELE_CALC_PORO_REAC_ECM_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc_poro_reac.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class StructPoroReactionECM;
}


namespace Discret
{
  namespace Elements
  {
    template <Core::FE::CellType distype>
    class ScaTraEleCalcPoroReacECM : public ScaTraEleCalcPoroReac<distype>
    {
     private:
      /// private constructor, since we are a Singleton.
      ScaTraEleCalcPoroReacECM(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcPoroReac<distype> pororeac;
      typedef ScaTraEleCalcPoro<distype> poro;
      typedef ScaTraEleCalcAdvReac<distype> advreac;
      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static ScaTraEleCalcPoroReacECM<distype>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! get the material parameters
      void get_material_params(
          const Core::Elements::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,          //!< density at t_(n)
          std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,         //!< density at t_(n+alpha_M)
          double& visc,                        //!< fluid viscosity
          const int iquad = -1                 //!< id of current gauss point (default = -1)
          ) override;

      //! evaluate material
      void materials(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! compute_struct_chem_potential
      double compute_struct_chem_potential(Mat::StructPoroReactionECM& structmat, const int gp);
    };

  }  // namespace Elements
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
