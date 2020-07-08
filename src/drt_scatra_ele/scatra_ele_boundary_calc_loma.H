/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for low Mach number problems


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_LOMA_H
#define SCATRA_ELE_BOUNDARY_CALC_LOMA_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcLoma : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcLoma<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcLoma* delete_me = NULL);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcLoma() { return; };

      //! singleton destruction
      virtual void Done();

      //! evaluate action
      int EvaluateAction(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, SCATRA::BoundaryAction action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcLoma(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate loma thermal press
      void CalcLomaThermPress(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! calculate Neumann inflow boundary conditions
      void NeumannInflow(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs);

      //! get density at integration point
      double GetDensity(Teuchos::RCP<const MAT::Material> material,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& ephinp, const int k);

      //! integral of normal diffusive flux and velocity over boundary surface
      void NormDiffFluxAndVelIntegral(const DRT::Element* ele, Teuchos::ParameterList& params,
          const std::vector<double>& enormdiffflux, const std::vector<double>& enormvel);

      //! thermodynamic pressure
      double thermpress_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
