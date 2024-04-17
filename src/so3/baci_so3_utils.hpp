/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
*-----------------------------------------------------------------------*/

#ifndef FOUR_C_SO3_UTILS_HPP
#define FOUR_C_SO3_UTILS_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class PreStress;

    namespace UTILS
    {
      template <CORE::FE::CellType distype>
      void CalcR(const DRT::Element* ele, const std::vector<double>& disp,
          CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::dim<distype>>& R);

      template <CORE::FE::CellType distype>
      void GetTemperatureForStructuralMaterial(
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, 1>& shapefctsGP,
          Teuchos::ParameterList& params);

      /*!
       * @brief compute the deformation gradient of the element at a specific Gauss point
       *
       * @tparam distype     shape of the element
       * @tparam probdim     dimension of the problem
       * @param[out] defgrd  deformation gradient
       * @param[in] nodes    list of nodes of the element
       * @param[in] xsi      position of the gauss point in parameter space
       * @param[in] xdisp    nodal displacements of the element
       */
      template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
      void ComputeDeformationGradient(CORE::LINALG::Matrix<probdim, probdim>& defgrd,
          DRT::Node** nodes, const CORE::LINALG::Matrix<probdim, 1>& xsi,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xdisp);

      /*!
       * @brief compute the deformation gradient of the element at a specific Gauss point
       *
       * @tparam distype          shape of the element
       * @tparam probdim          dimension of the problem
       * @param[out] defgrd       deformation gradient
       * @param[in] nodes         list of nodes of the element
       * @param[in] xsi           position of the gauss point in parameter space
       * @param[in] displacement  displacement vector of the element
       */
      template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
      void ComputeDeformationGradient(CORE::LINALG::Matrix<probdim, probdim>& defgrd,
          DRT::Node** nodes, const CORE::LINALG::Matrix<probdim, 1>& xsi,
          const std::vector<double>& displacement);

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
      template <CORE::FE::CellType distype>
      void ComputeDeformationGradient(
          CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::dim<distype>>& defgrd,
          INPAR::STR::KinemType kinemType,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, CORE::FE::dim<distype>>& xdisp,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, CORE::FE::dim<distype>>& xcurr,
          const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::dim<distype>>&
              inverseJacobian,
          const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>>& derivs,
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
      template <CORE::FE::CellType distype>
      void ComputeDeformationGradientMulf(
          CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::dim<distype>>& defgrd,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, CORE::FE::dim<distype>>& xdisp,
          const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>>& derivs,
          const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, int gp);

      /*!
       * \brief Compute the deformation gradient in case of no prestressing
       *
       * \tparam distype Shape of the element
       * \tparam probdim Dimension of the problem
       * \param defgrd [out] : Deformation gradient
       * \param xcurr[in] : Current nodal coordinates
       * \param derivs : Derivatives of the shape functions with respect to the reference
       * \param inverseJacobian [in] : Inverse jacobian at the point of evaluation
       */
      template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
      void ComputeDeformationGradientStandard(CORE::LINALG::Matrix<probdim, probdim>& defgrd,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xcurr,
          const CORE::LINALG::Matrix<probdim, CORE::FE::num_nodes<distype>>& derivs,
          const CORE::LINALG::Matrix<probdim, probdim>& inverseJacobian);

      /*!
       * \brief Evaluate the nodal coordinates of an element
       *
       * \tparam distype Shape of the element
       * \tparam probdim Dimension of the problem
       * \param nodes [in] : List of nodes of the element
       * \param xrefe [out] : reference coordinates of the element
       */
      template <CORE::FE::CellType distype, int probdim = 3>
      void EvaluateNodalCoordinates(
          DRT::Node** nodes, CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xrefe);

      /*!
       * \brief Evaluate the nodal displacement of the element
       *
       * \tparam distype  Shape of the element
       * \tparam probdim Dimension of the problem
       * \param disp [in] : Local displacement vector
       * \param xdisp [out] : Nodal displacements
       */
      template <CORE::FE::CellType distype, int probdim = 3>
      void EvaluateNodalDisplacements(const std::vector<double>& disp,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xdisp);

      /*!
       * \brief Compute the current nodal coordinates of the element
       *
       * \tparam distype Shape of the element
       * \tparam probdim Dimension of the problem
       * \param xrefe [in] : Reference coordinates of the element
       * \param xdisp [in] : Nodal displacements of the element
       * \param xcurr [out] : Current coordinates of the element
       */
      template <CORE::FE::CellType distype, int probdim = 3>
      void EvaluateCurrentNodalCoordinates(
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xdisp,
          CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, probdim>& xcurr);

      /*!
       * \brief Evaluates the inverse jacobian of the element
       *
       * \tparam distype Shape of the element
       * \param xrefe [in] : reference coordinates of the nodes
       * \param derivs [in] : Derivatives of the shape functions w.r.t. the reference coordinates
       * \param inverseJacobian [out] : inverse jacobian
       */
      template <CORE::FE::CellType distype>
      void EvaluateInverseJacobian(
          const CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, CORE::FE::dim<distype>>& xrefe,
          const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>>& derivs,
          CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::dim<distype>>& inverseJacobian);

      /*!
       * \brief Checks whether maerial tangent should be computed via finite difference. If so,
       * throws an error.
       *
       * \param sdyn [in] : Structural dynamics parameter list
       * \param eletype [out] : Element type string
       */
      void ThrowErrorFDMaterialTangent(
          const Teuchos::ParameterList& sdyn, const std::string& eletype);

    }  // namespace UTILS
  }    // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
