/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_ls.H

\brief evaluations for level sets

<pre>
Maintainer: Ursula Rasthofer
            erasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_LS_H
#define SCATRA_ELE_CALC_LS_H

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcLS : public ScaTraEleCalc<distype>
{
private:
    /// private constructor, since we are a Singleton.
    ScaTraEleCalcLS(const int numdofpernode,const int numscal);

    typedef ScaTraEleCalc<distype> my;

public:
    virtual ~ScaTraEleCalcLS() {}

    /// Singleton access method
    static ScaTraEleCalcLS<distype> * Instance(
      const int numdofpernode,
      const int numscal,
      bool create=true );

    /// called upon destruction
    virtual void Done();

    virtual int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization &         discretization,
                         const std::vector<int> &      lm,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra);

private:

    // calculate error compared to analytical solution
    void CalErrorComparedToAnalytSolution(
      const DRT::Element*                             ele,
      const std::vector<LINALG::Matrix<my::nen_,1> >& ephizero,
      Teuchos::ParameterList&                         params,
      Epetra_SerialDenseVector&                       errors);

    // smoothed heaviside function
    void SmoothHeavisideFunction(
      const double                      charelelength,
      const double                      phi,
      double&                           smoothH);

};
}
}

#endif

