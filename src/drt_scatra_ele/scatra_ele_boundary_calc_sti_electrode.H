/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_sti_electrode.H

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_H
#define SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_H

#include "../drt_scatra_ele/scatra_ele_boundary_calc.H"

namespace DRT
{
  // forward declaration
  class Discretization;

  namespace ELEMENTS
  {
    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleBoundaryCalcSTIElectrode : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

    public:

      //! singleton access method
      static ScaTraEleBoundaryCalcSTIElectrode<distype>* Instance(
          const int numdofpernode,
          const int numscal,
          const std::string& disname,
          const ScaTraEleBoundaryCalcSTIElectrode* delete_me = NULL
          );

      //! destructor
      virtual ~ScaTraEleBoundaryCalcSTIElectrode(){return;};

      //! singleton destruction
      virtual void Done();

    private:

      //! private constructor for singletons
      ScaTraEleBoundaryCalcSTIElectrode(
          const int numdofpernode,
          const int numscal,
          const std::string& disname
          );

      //! evaluate main-diagonal system matrix contributions associated with scatra-scatra interface coupling condition
      void EvaluateS2ICoupling(
          const DRT::FaceElement*        ele,              ///< current boundary element
          Teuchos::ParameterList&        params,           ///< parameter list
          DRT::Discretization&           discretization,   ///< discretization
          DRT::Element::LocationArray&   la,               ///< location array
          Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
          Epetra_SerialDenseVector&      eslaveresidual    ///< element residual for slave side
          );

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface coupling condition
      void EvaluateS2ICouplingOD(
          const DRT::FaceElement*        ele,              ///< current boundary element
          Teuchos::ParameterList&        params,           ///< parameter list
          DRT::Discretization&           discretization,   ///< discretization
          DRT::Element::LocationArray&   la,               ///< location array
          Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
          Epetra_SerialDenseMatrix&      emastermatrix     ///< element matrix for master side
          );

      //! evaluate action
      int EvaluateAction(
          DRT::FaceElement*              ele,              //!< boundary element
          Teuchos::ParameterList&        params,           //!< parameter list
          DRT::Discretization&           discretization,   //!< discretization
          SCATRA::BoundaryAction         action,           //!< action
          DRT::Element::LocationArray&   la,               //!< location array
          Epetra_SerialDenseMatrix&      elemat1_epetra,   //!< element matrix 1
          Epetra_SerialDenseMatrix&      elemat2_epetra,   //!< element matrix 2
          Epetra_SerialDenseVector&      elevec1_epetra,   //!< element right-hand side vector 1
          Epetra_SerialDenseVector&      elevec2_epetra,   //!< element right-hand side vector 2
          Epetra_SerialDenseVector&      elevec3_epetra    //!< element right-hand side vector 3
          );
    }; // class ScaTraEleBoundaryCalcSTIElectrode
  } // namespace ELEMENTS
} // namespace DRT
#endif
