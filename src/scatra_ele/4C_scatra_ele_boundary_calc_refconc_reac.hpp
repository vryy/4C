// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_REFCONC_REAC_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_REFCONC_REAC_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcRefConcReac : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef Discret::Elements::ScaTraEleBoundaryCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcRefConcReac<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! Factor needed for the calculation of reference concentrations
      double fac_for_ref_conc(const int iquad,      ///< current boundary integration point
          const Core::Elements::FaceElement* bele,  ///< current boundary element
          Teuchos::ParameterList& params,           ///< parameter list
          Core::FE::Discretization& discretization  ///< discretization
          ) override;

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcRefConcReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      template <Core::FE::CellType bdistype,
          Core::FE::CellType pdistype>
      double calc_jat_int_point(const int iquad,    ///< current boundary integration point
          const Core::Elements::FaceElement* bele,  ///< current boundary element
          const Core::Elements::Element* pele,      ///< current parent element
          Teuchos::ParameterList& params,           ///< parameter list
          Core::FE::Discretization& discretization  ///< discretization
      );

    };  // class ScaTraEleBoundaryCalcRefConcReac

  }  // namespace Elements

}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
