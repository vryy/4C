/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
*-----------------------------------------------------------------------*/

#ifndef SO_UTILS_H
#define SO_UTILS_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
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

      /*!
       * \brief Compute the deformation gradient of the element at a specific Gauss point (including
       * the MULF switch)
       *
       * \tparam distype Shape of the element
       * \param defgrd [out] : Deformation gradient
       * \param kinemType [in] : Type of kinematics
       * \param xdisp [in] : Nodal displacements
       * \param xcurr [in] : Current nodal coordinates
       * \param inverseJacobian [in] : Inverse jacobian at the point of evaluation
       * \param derivs [in] Derivatives of the shape functions evaluated at the evaluation point
       * \param prestressType [in] :
       * \param mulfHistory [in] : Internal MULF history variables
       * \param gp [in] : Gauss point
       */
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

      /*!
       * \brief Compute the deformation gradient in the case of MULF
       *
       * \tparam distype Shape of the element
       * \param defgrd [out] : Deformation gradient
       * \param xdisp [in] : Nodal displacements of the element
       * \param derivs [in] : Derivatives of the shape functions with respect to the reference
       * coordinates
       *
       * \param mulfHistory [in] : Internal MULF history variables
       * \param gp [in] : Gauss point
       */
      template <DRT::Element::DiscretizationType distype>
      void ComputeDeformationGradientMulf(LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
                                              DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, int gp);

      /*!
       * \brief Compute the deformation gradient in case of no prestressing
       *
       * \tparam distype Shape of the element
       * \param defgrd [out] : Deformation gradient
       * \param xcurr[in] : Current nodal coordinates
       * \param derivs : Derivatives of the shape functions with respect to the reference
       * \param inverseJacobian [in] : Inverse jacobian at the point of evaluation
       */
      template <DRT::Element::DiscretizationType distype>
      void ComputeDeformationGradientStandard(LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
                                                  DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian);

      /*!
       * \brief Evaluate the nodal coordinates of an element
       *
       * \tparam distype Shape of the element
       * \param nodes [in] : List of nodes of the element
       * \param xrefe [out] : reference coordinates of the element
       */
      template <DRT::Element::DiscretizationType distype>
      void EvaluateNodalCoordinates(DRT::Node** nodes,
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe);

      /*!
       * \brief Evaluate the nodal displacement of the element
       *
       * \tparam distype  Shape of the element
       * \param disp [in] : Local displacement vector
       * \param xdisp [out] : Nodal displacements
       */
      template <DRT::Element::DiscretizationType distype>
      void EvaluateNodalDisplacements(const std::vector<double>& disp,
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp);

      /*!
       * \brief Compute the current nodal coordinates of the element
       *
       * \tparam distype Shape of the element
       * \param xrefe [in] : Reference coordinates of the element
       * \param xdisp [in] : Nodal displacements of the element
       * \param xcurr [out] : Current coordinates of the element
       */
      template <DRT::Element::DiscretizationType distype>
      void EvaluateCurrentNodalCoordinates(
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
          LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr);

      /*!
       * \brief Evaluates the inverse jacobian of the element
       *
       * \tparam distype Shape of the element
       * \param xrefe [in] : reference coordinates of the nodes
       * \param derivs [in] : Derivatives of the shape functions w.r.t. the reference coordinates
       * \param inverseJacobian [out] : inverse jacobian
       */
      template <DRT::Element::DiscretizationType distype>
      void EvaluateInverseJacobian(
          const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
              DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
          const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
          LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian);

      /*!
       * \brief Check whether all nodes of the element have fibers so that the fibers
       *
       * \tparam distype Discretization type
       * \param nodes Pointer to the nodes of the element
       * \return true All nodes have fibers
       * \return false At least one one does not have a fiber
       */
      template <DRT::Element::DiscretizationType distype>
      bool HaveNodalFibers(DRT::Node** nodes);
    }  // namespace UTILS
  }    // namespace ELEMENTS
}  // namespace DRT

#endif
