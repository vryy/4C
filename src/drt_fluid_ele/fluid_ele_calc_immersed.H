/*!----------------------------------------------------------------------
\file fluid_ele_calc_immersed.H

\brief calc class for immersed problems

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240
</pre>
*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_IMMERSED_H
#define FLUID_ELE_CALC_IMMERSED_H

#include "fluid_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{


template<DRT::Element::DiscretizationType distype>
class FluidEleCalcImmersed : public FluidEleCalc<distype>
{
  typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

protected:
  /// private Constructor since we are a Singleton.
  FluidEleCalcImmersed();

public:
  virtual ~FluidEleCalcImmersed() {}

  /// Singleton access method
  static FluidEleCalcImmersed<distype>* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /*!
    Evaluate

    \param eid              (i) element id
    \param discretization   (i) fluid discretization the element belongs to
    \param lm               (i) location matrix of element
    \param params           (i) element parameter list
    \param mat              (i) material
    \param elemat1_epetra   (o) element matrix to calculate
    \param elemat2_epetra   (o) element matrix to calculate
    \param elevec1_epetra   (o) element vector to calculate
    \param elevec2_epetra   (o) element vector to calculate
    \param elevec3_epetra   (o) element vector to calculate
    \param intpoints        (i) Gaussian integration points

   */
  virtual int Evaluate(DRT::ELEMENTS::Fluid*                    ele,
      DRT::Discretization &                discretization,
      const std::vector<int> &             lm,
      Teuchos::ParameterList&              params,
      Teuchos::RCP<MAT::Material> &        mat,
      Epetra_SerialDenseMatrix&            elemat1_epetra,
      Epetra_SerialDenseMatrix&            elemat2_epetra,
      Epetra_SerialDenseVector&            elevec1_epetra,
      Epetra_SerialDenseVector&            elevec2_epetra,
      Epetra_SerialDenseVector&            elevec3_epetra,
      const DRT::UTILS::GaussIntegration & intpoints);

};
}
}

#endif

