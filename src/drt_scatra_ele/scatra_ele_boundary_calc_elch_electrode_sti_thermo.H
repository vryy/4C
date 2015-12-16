/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode_sti_thermo.H

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_STI_THERMO_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_STI_THERMO_H

#include "../drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode.H"

namespace DRT
{
  // forward declaration
  class Discretization;

  namespace ELEMENTS
  {
    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleBoundaryCalcElchElectrodeSTIThermo : public ScaTraEleBoundaryCalcElchElectrode<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype> myelectrode;

    public:

      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>* Instance(
          const int numdofpernode,
          const int numscal,
          const std::string& disname,
          const ScaTraEleBoundaryCalcElchElectrodeSTIThermo* delete_me = NULL
          );

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchElectrodeSTIThermo(){return;};

      //! singleton destruction
      virtual void Done();

    private:

      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchElectrodeSTIThermo(
          const int numdofpernode,
          const int numscal,
          const std::string& disname
          );

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface coupling condition
      void EvaluateS2ICouplingOD(
          const DRT::FaceElement*        ele,              ///< current boundary element
          Teuchos::ParameterList&        params,           ///< parameter list
          DRT::Discretization&           discretization,   ///< discretization
          DRT::Element::LocationArray&   la,               ///< location array
          Epetra_SerialDenseMatrix&      eslavematrix      ///< element matrix for slave side
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

      //! evaluate factor F/RT
      const double GetFRT(
          DRT::Discretization&           discretization,   ///< discretization
          DRT::Element::LocationArray&   la                ///< location array
          ) const;

      //! evaluate factor F/RT²
      const double GetFRTT(
          DRT::Discretization&           discretization,   ///< discretization
          DRT::Element::LocationArray&   la                ///< location array
          ) const;

      //! evaluate local temperature value
      const double GetTemperature(
          DRT::Discretization&           discretization,   ///< discretization
          DRT::Element::LocationArray&   la                ///< location array
          ) const;
    }; // class ScaTraEleBoundaryCalcElchElectrodeSTIThermo
  } // namespace ELEMENTS
} // namespace DRT
#endif
