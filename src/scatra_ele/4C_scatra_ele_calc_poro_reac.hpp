// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_PORO_REAC_HPP
#define FOUR_C_SCATRA_ELE_CALC_PORO_REAC_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_calc_advanced_reaction.hpp"
#include "4C_scatra_ele_calc_poro.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    template <Core::FE::CellType distype>
    class ScaTraEleCalcPoroReac : virtual public ScaTraEleCalcPoro<distype>,
                                  virtual public ScaTraEleCalcAdvReac<distype>
    {
     protected:
      /// protected constructor, since we are a Singleton.
      ScaTraEleCalcPoroReac(const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcPoro<distype> poro;
      typedef ScaTraEleCalcAdvReac<distype> advreac;

     public:
      /// Singleton access method
      static ScaTraEleCalcPoroReac<distype>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     protected:
      //! evaluate material
      void materials(const std::shared_ptr<const Core::Mat::Material>
                         material,  //!< pointer to current material
          const int k,              //!< id of current scalar
          double& densn,            //!< density at t_(n)
          double& densnp,           //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,           //!< density at t_(n+alpha_M)
          double& visc,             //!< fluid viscosity
          const int iquad = -1      //!< id of current gauss point (default = -1)
          ) override;

      //! get the material parameters
      void get_material_params(
          const Core::Elements::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,          //!< density at t_(n)
          std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,         //!< density at t_(n+alpha_M)
          double& visc,                        //!< fluid viscosity
          const int iquad = -1                 //!< id of current gauss point (default = -1)
          ) override;

      //! material ScaTra
      void mat_scatra(const std::shared_ptr<const Core::Mat::Material>
                          material,  //!< pointer to current material
          const int k,               //!< id of current scalar
          double& densn,             //!< density at t_(n)
          double& densnp,            //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,            //!< density at t_(n+alpha_M)
          double& visc,              //!< fluid viscosity
          const int iquad = -1       //!< id of current gauss point (default = -1)
          ) override;

      //! extract element based or nodal values
      //  return extracted values of phinp
      void extract_element_and_node_values(Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::LocationArray& la) override;
    };
  }  // namespace Elements
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
