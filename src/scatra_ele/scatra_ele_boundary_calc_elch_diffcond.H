/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation


\level 2
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_H

#include "scatra_ele_boundary_calc_elch_electrode.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchDiffCond;
    class ScaTraEleParameterElchDiffCond;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcElchDiffCond : public ScaTraEleBoundaryCalcElchElectrode<distype>
    {
      using my = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>;
      using myelch = DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>;
      using myelectrode = DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>;

     protected:
      using my::nen_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchDiffCond<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);



     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchDiffCond(
          const int numdofpernode, const int numscal, const std::string& disname);

      int EvaluateAction(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, SCATRA::BoundaryAction action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra) override;

      int EvaluateNeumann(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition,
          DRT::Element::LocationArray& la, Epetra_SerialDenseVector& elevec1,
          const double scalar) override;

      void EvaluateElchBoundaryKinetics(const DRT::Element* ele, Epetra_SerialDenseMatrix& emat,
          Epetra_SerialDenseVector& erhs, const std::vector<LINALG::Matrix<nen_, 1>>& ephinp,
          const std::vector<LINALG::Matrix<nen_, 1>>& ehist, double timefac,
          Teuchos::RCP<const MAT::Material> material, Teuchos::RCP<DRT::Condition> cond,
          const int nume, const std::vector<int> stoich, const int kinetics, const double pot0,
          const double frt, const double scalar) override;

      void EvaluateS2ICoupling(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix,
          Epetra_SerialDenseVector& eslaveresidual) override;

      void EvaluateS2ICouplingOD(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& eslavematrix) override;

      double GetValence(
          const Teuchos::RCP<const MAT::Material>& material, const int k) const override;

      //! diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc_;
    };  // class ScaTraEleBoundaryCalcElchDiffCond
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
