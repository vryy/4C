/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for low Mach number problems


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_LOMA_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_LOMA_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_boundary_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcLoma : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcLoma<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


      //! evaluate action
      int EvaluateAction(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, SCATRA::BoundaryAction action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

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
          CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs) override;

      //! get density at integration point
      double GetDensity(Teuchos::RCP<const MAT::Material> material,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp, const int k) override;

      //! integral of normal diffusive flux and velocity over boundary surface
      void NormDiffFluxAndVelIntegral(const DRT::Element* ele, Teuchos::ParameterList& params,
          const std::vector<double>& enormdiffflux, const std::vector<double>& enormvel);

      //! thermodynamic pressure
      double thermpress_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif
