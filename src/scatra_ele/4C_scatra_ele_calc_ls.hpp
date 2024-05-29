/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluations for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_LS_HPP
#define FOUR_C_SCATRA_ELE_CALC_LS_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcLS : public ScaTraEleCalc<distype>
    {
     private:
      //! private constructor for singletons
      ScaTraEleCalcLS(const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static ScaTraEleCalcLS<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      int evaluate_action(CORE::Elements::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          CORE::Elements::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     private:
      // calculate error compared to analytical solution
      void cal_error_compared_to_analyt_solution(const CORE::Elements::Element* ele,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephizero,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& errors);

      // smoothed heaviside function
      void smooth_heaviside_function(const double charelelength, const double phi, double& smoothH);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
