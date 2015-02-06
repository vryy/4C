/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_std.H

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_STD_H
#define SCATRA_ELE_BOUNDARY_CALC_STD_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{

namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleBoundaryCalcStd : public ScaTraEleBoundaryCalc<distype>
{
  typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

public:

  /// Singleton access method
  static ScaTraEleBoundaryCalcStd<distype> * Instance(const int numdofpernode, const int numscal, bool create=true );

  /// Empty destructor
  virtual ~ScaTraEleBoundaryCalcStd(){return;};

  /// called upon destruction
  void Done();

private:
  /// private constructor since we are singleton
  ScaTraEleBoundaryCalcStd(const int numdofpernode, const int numscal);

  int EvaluateAction(
      DRT::ELEMENTS::TransportBoundary* ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm,
      Epetra_SerialDenseMatrix&         elemat1_epetra,
      Epetra_SerialDenseMatrix&         elemat2_epetra,
      Epetra_SerialDenseVector&         elevec1_epetra,
      Epetra_SerialDenseVector&         elevec2_epetra,
      Epetra_SerialDenseVector&         elevec3_epetra
  );

}; // class ScaTraEleBoundaryCalcStd

} //namsepace ELEMENTS
} // namespace DRT

#endif

