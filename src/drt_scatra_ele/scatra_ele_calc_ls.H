/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluations for level sets

\level 2

\maintainer Christoph Ager
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_LS_H
#define SCATRA_ELE_CALC_LS_H

#include "scatra_ele_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcLS : public ScaTraEleCalc<distype>
    {
     private:
      //! private constructor for singletons
      ScaTraEleCalcLS(const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;

     public:
      virtual ~ScaTraEleCalcLS() {}

      /// Singleton access method
      static ScaTraEleCalcLS<distype>* Instance(const int numdofpernode, const int numscal,
          const std::string& disname, const ScaTraEleCalcLS* delete_me = NULL);

      /// called upon destruction
      virtual void Done();

      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     private:
      // calculate error compared to analytical solution
      void CalErrorComparedToAnalytSolution(const DRT::Element* ele,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& ephizero, Teuchos::ParameterList& params,
          Epetra_SerialDenseVector& errors);

      // smoothed heaviside function
      void SmoothHeavisideFunction(const double charelelength, const double phi, double& smoothH);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
