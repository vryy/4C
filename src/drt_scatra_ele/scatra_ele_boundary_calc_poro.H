/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_std.H

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_PORO_H
#define SCATRA_ELE_BOUNDARY_CALC_PORO_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{

namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleBoundaryCalcPoro : public ScaTraBoundaryImpl<distype>
{
  typedef DRT::ELEMENTS::ScaTraBoundaryImpl<distype> my;

public:

  /// Singleton access method
  static ScaTraEleBoundaryCalcPoro<distype> * Instance(const int numdofpernode, const int numscal, bool create=true );

  /// Empty destructor
  virtual ~ScaTraEleBoundaryCalcPoro() {};

  /// called upon destruction
  void Done();

private:
  /// private constructor since we are singleton
  ScaTraEleBoundaryCalcPoro(const int numdofpernode, const int numscal);

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

  void EvaluateFPSSICoupling(
      DRT::ELEMENTS::TransportBoundary* ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm,
      Epetra_SerialDenseMatrix&         elemat1_epetra,
      Epetra_SerialDenseMatrix&         elemat2_epetra,
      Epetra_SerialDenseVector&         elevec1_epetra,
      Epetra_SerialDenseVector&         elevec2_epetra,
      Epetra_SerialDenseVector&         elevec3_epetra);

}; // class ScaTraEleBoundaryCalcPoro

} //namsepace ELEMENTS
} // namespace DRT

#endif

