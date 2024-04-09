/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_boundary_calc_elch_electrode.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchDiffCond;
    class ScaTraEleParameterElchDiffCond;

    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchDiffCond
        : public ScaTraEleBoundaryCalcElchElectrode<distype, probdim>
    {
      using my = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelch = DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>;
      using myelectrode = DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>;

     protected:
      using my::nen_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);



     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchDiffCond(
          const int numdofpernode, const int numscal, const std::string& disname);

      int EvaluateAction(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, SCATRA::BoundaryAction action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

      int EvaluateNeumann(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseVector& elevec1,
          const double scalar) override;

      void EvaluateElchBoundaryKinetics(const DRT::Element* ele,
          CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist, double timefac,
          Teuchos::RCP<const MAT::Material> material, Teuchos::RCP<DRT::Condition> cond,
          const int nume, const std::vector<int> stoich, const int kinetics, const double pot0,
          const double frt, const double scalar) override;

      void EvaluateS2ICoupling(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& eslavematrix,
          CORE::LINALG::SerialDenseMatrix& emastermatrix,
          CORE::LINALG::SerialDenseVector& eslaveresidual) override;

      void EvaluateS2ICouplingOD(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& eslavematrix) override;

      double GetValence(
          const Teuchos::RCP<const MAT::Material>& material, const int k) const override;

      //! diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc_;
    };  // class ScaTraEleBoundaryCalcElchDiffCond
  }     // namespace ELEMENTS
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif
