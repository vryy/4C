/*----------------------------------------------------------------------*/
/*! \file
\brief main file containing routines for calculation of scatra element formulated in reference
concentrations and with advanced reaction terms

\level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_REFCONC_REAC_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_REFCONC_REAC_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcRefConcReac : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcRefConcReac<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! Factor needed for the calculation of reference concentrations
      double FacForRefConc(const int iquad,    ///< current boundary integration point
          const DRT::FaceElement* bele,        ///< current boundary element
          Teuchos::ParameterList& params,      ///< parameter list
          DRT::Discretization& discretization  ///< discretization
          ) override;

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcRefConcReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      template <CORE::FE::CellType bdistype,
          CORE::FE::CellType pdistype>
      double CalcJatIntPoint(const int iquad,  ///< current boundary integration point
          const DRT::FaceElement* bele,        ///< current boundary element
          const DRT::Element* pele,            ///< current parent element
          Teuchos::ParameterList& params,      ///< parameter list
          DRT::Discretization& discretization  ///< discretization
      );

    };  // class ScaTraEleBoundaryCalcRefConcReac

  }  // namespace ELEMENTS

}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
