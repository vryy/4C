/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
*-----------------------------------------------------------------------*/

#ifndef SO_UTILS_H
#define SO_UTILS_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"  // namespace INPAR
#include "../drt_inpar/inpar_structure.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class PreStress;

    namespace UTILS
    {
      template <DRT::Element::DiscretizationType distype>
      void CalcR(const DRT::Element* ele, const std::vector<double>& disp,
          LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& R);

      template <DRT::Element::DiscretizationType distype>
      void NodalFiber(DRT::Node** nodes,
          const std::vector<
              LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>>&
              shapefcts,
          std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1>>& fibergp,
          std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1>>& sheetgp);

      template <DRT::Element::DiscretizationType distype>
      void GetTemperatureForStructuralMaterial(
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>&
              shapefctsGP,
          Teuchos::ParameterList& params);

      template <DRT::Element::DiscretizationType distype>
      void ComputeDeformationGradient(LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
                                          DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
          INPAR::STR::KinemType kinemType,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          const INPAR::STR::PreStress prestressType,
          const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, int gp);

      template <DRT::Element::DiscretizationType distype>
      void ComputeDeformationGradientMulf(LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
                                              DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, int gp);

      template <DRT::Element::DiscretizationType distype>
      void ComputeDeformationGradientStandard(LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
                                                  DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian);

      template <DRT::Element::DiscretizationType distype>
      void EvaluateNodalCoordinates(DRT::Node** nodes,
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe);

      template <DRT::Element::DiscretizationType distype>
      void EvaluateNodalDisplacements(const std::vector<double>& disp,
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp);

      template <DRT::Element::DiscretizationType distype>
      void EvaluateCurrentNodalCoordinates(
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr);


      template <DRT::Element::DiscretizationType distype>
      void EvaluateInverseJacobian(
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian);
    }  // namespace UTILS
  }    // namespace ELEMENTS
}  // namespace DRT

#endif
