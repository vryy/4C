/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       with EAS element technology
\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_CALC_EAS_H
#define SOLID_ELE_CALC_EAS_H

#include "solid_ele_calc.H"
#include "element.H"
#include "inpar_structure.H"
#include "solid_ele_eas_utils.H"

namespace DRT::ELEMENTS
{
  template <DRT::Element::DiscretizationType distype, int neas>
  class SolidEleCalcEas : public virtual SolidEleCalc<distype>
  {
   public:
    /// private constructor, since we are a Singleton.
    SolidEleCalcEas();
    /// evaluate element
    int nln_force_stiff_mass(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
        const std::vector<int>& lm, Teuchos::ParameterList& params,
        Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
        Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
        Epetra_SerialDenseVector* elevec3_epetra,
        DRT::UTILS::GaussIntegration* intpoints = nullptr) override;

    /// evaluate with GEMM time integration
    int nln_force_stiff_mass_gemm(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
        const std::vector<int>& lm, Teuchos::ParameterList& params,
        Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
        Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
        Epetra_SerialDenseVector* elevec3_epetra,
        DRT::UTILS::GaussIntegration* intpoints = nullptr) override;

    /// Update the element at converged time step
    int UpdateElement(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
        const std::vector<int>& lm, Teuchos::ParameterList& params,
        Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
        Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
        Epetra_SerialDenseVector* elevec3_epetra) override;

    /// recover variables condensed at element level
    int RecoverCondensed(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
        const std::vector<int>& lm, Teuchos::ParameterList& params,
        Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
        Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
        Epetra_SerialDenseVector* elevec3_epetra) override;

    /// setup element
    virtual void Setup(DRT::ELEMENTS::Solid* ele, DRT::INPUT::LineDefinition* linedef);

   protected:
    /// Prepare terms for EAS evaluation at before GP loop
    void PrepareEAS(DRT::ELEMENTS::Solid* ele);

    /// integrate the EAS terms
    void IntegrateEAS(DRT::ELEMENTS::Solid* ele);

    /// integrate the EAS terms
    void EvaluateEASshape(const ::STR::ELEMENTS::EASType eastype);

    /// integrate the EAS terms
    void EnhanceEASstrains(DRT::ELEMENTS::Solid* ele);

    using my = SolidEleCalc<distype>;

   protected:
    LINALG::Matrix<my::numstr_, my::numstr_> T0invT_;
    LINALG::Matrix<my::numstr_, neas> shape_eas_;
    /// deformation gradient consistent to enhanced strains
    LINALG::Matrix<my::nsd_, my::nsd_> defgrd_enh_;

  };  // class SolidEleCalcEas
}  // namespace DRT::ELEMENTS

#endif  // SOLID_ELE_CALC_EAS_H
