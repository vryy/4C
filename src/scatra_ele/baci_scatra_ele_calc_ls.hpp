/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluations for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_LS_HPP
#define FOUR_C_SCATRA_ELE_CALC_LS_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN

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

      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     private:
      // calculate error compared to analytical solution
      void CalErrorComparedToAnalytSolution(const DRT::Element* ele,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephizero,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& errors);

      // smoothed heaviside function
      void SmoothHeavisideFunction(const double charelelength, const double phi, double& smoothH);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
