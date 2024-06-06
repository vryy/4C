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

namespace Discret
{
  namespace ELEMENTS
  {
    template <Core::FE::CellType distype>
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

      int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, const ScaTra::Action& action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

     private:
      // calculate error compared to analytical solution
      void cal_error_compared_to_analyt_solution(const Core::Elements::Element* ele,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephizero,
          Teuchos::ParameterList& params, Core::LinAlg::SerialDenseVector& errors);

      // smoothed heaviside function
      void smooth_heaviside_function(const double charelelength, const double phi, double& smoothH);
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
