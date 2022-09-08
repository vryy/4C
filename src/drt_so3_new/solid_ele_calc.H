/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       simple displacement based
\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_CALC_H
#define SOLID_ELE_CALC_H

#include "solid_ele_interface.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "solid_utils.H"
#include "../drt_lib/drt_element_integration_select.H"

namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace DRT
{
  namespace UTILS
  {
    template <DRT::Element::DiscretizationType>
    struct DisTypeToNumNodePerEle;
    template <DRT::Element::DiscretizationType>
    struct DisTypeToDim;
  }  // namespace UTILS

  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class SolidEleCalc : public virtual SolidEleInterface
    {
     public:
      SolidEleCalc();

      virtual ~SolidEleCalc() = default;
      /// Singleton access method
      static SolidEleCalc<distype>* Instance(bool create = true);

      /// lump mass matrix
      void LumpMass(Epetra_SerialDenseMatrix* mass) override;

      /// setup element
      void Setup(DRT::ELEMENTS::Solid* ele, DRT::INPUT::LineDefinition* linedef) override;

      /*!
      \brief Evaluate an element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(DRT::ELEMENTS::Solid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3) override;

      /// Evaluate Neumann condition on this element
      int EvaluateNeumann(DRT::ELEMENTS::Solid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          const DRT::UTILS::GaussIntegration* intpoints, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = NULL) override;

      /// post-process stress
      int PostProcessStress(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override;

      /// calculate stress for output
      int CalcStress(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override;

      /// Update the element at converged time step
      int UpdateElement(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override;

      /// Evaluate the element for GEMM time integration
      /*!
        This class does not provide a definition for this function; it
        must be defined in SolidEleCalc.
       */
      int nln_force_stiff_mass_gemm(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra,
          DRT::UTILS::GaussIntegration* intpoints = nullptr) override
      {
        dserror("todo");
        return 0;
      }

      int ResetStep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override
      {
        dserror("todo");
        return 0;
      }

      int StoreStep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override
      {
        dserror("todo");
        return 0;
      }

      int RecoverStep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override
      {
        dserror("todo");
        return 0;
      }

      int ResetIstep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override
      {
        dserror("todo");
        return 0;
      }

      int GlobalGPstressMap(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override
      {
        dserror("todo");
        return 0;
      }

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in SolidEleCalc.
       */
      int nln_force_stiff_mass(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra,
          DRT::UTILS::GaussIntegration* intpoints = nullptr) override;

      /// recover variables condensed at element level
      int RecoverCondensed(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) override
      {
        return 0;
      };



     protected:
      /// static values for matrix sizes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
      static const int numdofperelement_ = nen_ * nsd_;
      static const int numstr_ = nsd_ * (nsd_ + 1) / 2;

      Teuchos::RCP<DRT::UTILS::GaussIntegration> default_integration_;

      virtual void EvaluateShapeDeriv(const DRT::UTILS::GaussIntegration& intpoints, const int gp);

      virtual void JacobianMapping(const DRT::UTILS::GaussIntegration& intpoints, const int gp);

      virtual void Strains();

      virtual void StrainGradient();

      virtual void AddInternalForce(LINALG::Matrix<nsd_ * nen_, 1>& force, const double& fac);

      virtual void AddElasticStiff(
          LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& stiff, const double& fac);

      virtual void AddGeometricStiff(
          LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& stiff, const double& fac);

      virtual void AddMass(LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& mass, const double& fac);

      virtual void Material(LINALG::Matrix<nsd_, nsd_>* defgrd,
          LINALG::Matrix<numstr_, 1>* glstrain, DRT::ELEMENTS::Solid* ele,
          Teuchos::ParameterList& params, int gp);

      virtual void FillPositionArray(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm);

      virtual void CalcStressStrainOutput(DRT::ELEMENTS::Solid* ele, const int gp);
      virtual void WriteStressStrainOutput(DRT::ELEMENTS::Solid* ele);

      /// initialize default GP locations, shape functions and derivs
      // todo: replace this with template "distype to default rule"
      //       once the DRT::UTILS::GaussIntegration does no longer use Intrepid
      virtual void InitializeDefaultQuadrature();


      //! interface ptr for data exchange between the element and the time integrator.
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      /// Integration factor: detJ * w(gp)
      double fac_;

      /// jacobian determinant
      double detJ_;

      /// matrices refilled for every element or gp, to avoid alloc/dealloc
      // X, material coord. of element
      LINALG::Matrix<nen_, nsd_> xrefe_;
      // x, current  coord. of element
      LINALG::Matrix<nen_, nsd_> xcurr_;
      /// shape function
      LINALG::Matrix<nen_, 1> shapefunction_;
      /// shape function deriv
      LINALG::Matrix<nsd_, nen_> deriv_;
      /// deformation gradient
      LINALG::Matrix<nsd_, nsd_> defgrd_;
      /// shape function derivative w.r.t. XYZ
      LINALG::Matrix<nsd_, nen_> n_xyz_;
      /// inverse jacobian
      LINALG::Matrix<nsd_, nsd_> invJ_;
      /// right cauchy-green
      LINALG::Matrix<nsd_, nsd_> rcg_;
      /// Green-Lagrange strain
      LINALG::Matrix<numstr_, 1> gl_strain_;
      /// gp coords
      LINALG::Matrix<nsd_, 1> xi_;
      /// B-Operator
      LINALG::Matrix<numstr_, nsd_ * nen_> bop_;
      /// material tangent
      LINALG::Matrix<numstr_, numstr_> cmat_;
      /// 2pk stress
      LINALG::Matrix<numstr_, 1> pk2_;
      /// matrix cmat_*bop_
      LINALG::Matrix<numstr_, nen_ * nsd_> cb_;



      /// stress output
      LINALG::Matrix<
          STR::UTILS::IntRuleToNquad<DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule>::ngp,
          numstr_>
          stress_output_;
      /// strain output
      LINALG::Matrix<
          STR::UTILS::IntRuleToNquad<DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule>::ngp,
          numstr_>
          strain_output_;

    };  // class SolidEleCalc
  }     // namespace ELEMENTS
}  // namespace DRT

#endif  // SOLID_ELE_CALC_STD_H
