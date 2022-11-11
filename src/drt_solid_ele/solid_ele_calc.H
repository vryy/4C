/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       simple displacement based
\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_CALC_H
#define SOLID_ELE_CALC_H

#include <memory>
#include <string>
#include <unordered_map>
#include "solid_ele_interface.H"
#include "drt_element.H"
#include "drt_utils_gausspoints.H"
#include "solid_utils.H"
#include "drt_element_integration_select.H"

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

      /// Singleton access method
      static SolidEleCalc<distype>* Instance(bool create = true);

      void Setup(const DRT::ELEMENTS::Solid& ele, DRT::INPUT::LineDefinition* linedef) override;
      void MaterialPostSetup(const DRT::ELEMENTS::Solid& ele) override;

      void EvaluateNonlinearForceStiffnessMass(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseVector* force_vector,
          Epetra_SerialDenseMatrix* stiffness_matrix,
          Epetra_SerialDenseMatrix* mass_matrix) override;

      void Recover(const DRT::ELEMENTS::Solid& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params) override;

      void CalculateStress(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) override;

      void PostProcessStressStrain(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) override;

      void Update(const DRT::ELEMENTS::Solid& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params) override;

      void InitializeGaussPointDataOutput(const DRT::ELEMENTS::Solid& ele) const override;

      void EvaluateGaussPointDataOutput(const DRT::ELEMENTS::Solid& ele) const override;


     private:
      /// static values for matrix sizes
      static constexpr int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static constexpr int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
      static constexpr int numdofperelement_ = nen_ * nsd_;
      static constexpr int numstr_ = nsd_ * (nsd_ + 1) / 2;

      std::unique_ptr<DRT::UTILS::GaussIntegration> stiffness_matrix_integration_;
      std::unique_ptr<DRT::UTILS::GaussIntegration> mass_matrix_integration_;



      /// initialize default GP locations, shape functions and derivs
      // todo: replace this with template "distype to default rule"
      //       once the DRT::UTILS::GaussIntegration does no longer use Intrepid
      virtual void InitializeDefaultQuadrature();


      Teuchos::RCP<DRT::UTILS::GaussIntegration> default_integration_ = Teuchos::null;

    };  // class SolidEleCalc
  }     // namespace ELEMENTS
}  // namespace DRT

#endif  // SOLID_ELE_CALC_STD_H
