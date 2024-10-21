#ifndef FOUR_C_SO3_UTILS_HPP
#define FOUR_C_SO3_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class PreStress;

    namespace Utils
    {
      template <Core::FE::CellType distype>
      void calc_r(const Core::Elements::Element* ele, const std::vector<double>& disp,
          Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& R);

      template <Core::FE::CellType distype>
      void get_temperature_for_structural_material(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, 1>& shapefctsGP,
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
      template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>
      void compute_deformation_gradient(Core::LinAlg::Matrix<probdim, probdim>& defgrd,
          Core::Nodes::Node** nodes, const Core::LinAlg::Matrix<probdim, 1>& xsi,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xdisp);

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
      template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>
      void compute_deformation_gradient(Core::LinAlg::Matrix<probdim, probdim>& defgrd,
          Core::Nodes::Node** nodes, const Core::LinAlg::Matrix<probdim, 1>& xsi,
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
      template <Core::FE::CellType distype>
      void compute_deformation_gradient(
          Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& defgrd,
          Inpar::Solid::KinemType kinemType,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xdisp,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xcurr,
          const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>&
              inverseJacobian,
          const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& derivs,
          const Inpar::Solid::PreStress prestressType, Discret::ELEMENTS::PreStress& mulfHistory,
          int gp);

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
      template <Core::FE::CellType distype>
      void compute_deformation_gradient_mulf(
          Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& defgrd,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xdisp,
          const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& derivs,
          Discret::ELEMENTS::PreStress& mulfHistory, int gp);

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
      template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>
      void compute_deformation_gradient_standard(Core::LinAlg::Matrix<probdim, probdim>& defgrd,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xcurr,
          const Core::LinAlg::Matrix<probdim, Core::FE::num_nodes<distype>>& derivs,
          const Core::LinAlg::Matrix<probdim, probdim>& inverseJacobian);

      /*!
       * \brief Evaluate the nodal coordinates of an element
       *
       * \tparam distype Shape of the element
       * \tparam probdim Dimension of the problem
       * \param nodes [in] : List of nodes of the element
       * \param xrefe [out] : reference coordinates of the element
       */
      template <Core::FE::CellType distype, int probdim = 3>
      void evaluate_nodal_coordinates(Core::Nodes::Node** nodes,
          Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xrefe);

      /*!
       * \brief Evaluate the nodal displacement of the element
       *
       * \tparam distype  Shape of the element
       * \tparam probdim Dimension of the problem
       * \param disp [in] : Local displacement vector
       * \param xdisp [out] : Nodal displacements
       */
      template <Core::FE::CellType distype, int probdim = 3>
      void evaluate_nodal_displacements(const std::vector<double>& disp,
          Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xdisp);

      /*!
       * \brief Compute the current nodal coordinates of the element
       *
       * \tparam distype Shape of the element
       * \tparam probdim Dimension of the problem
       * \param xrefe [in] : Reference coordinates of the element
       * \param xdisp [in] : Nodal displacements of the element
       * \param xcurr [out] : Current coordinates of the element
       */
      template <Core::FE::CellType distype, int probdim = 3>
      void evaluate_current_nodal_coordinates(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xdisp,
          Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xcurr);

      /*!
       * \brief Evaluates the inverse jacobian of the element
       *
       * \tparam distype Shape of the element
       * \param xrefe [in] : reference coordinates of the nodes
       * \param derivs [in] : Derivatives of the shape functions w.r.t. the reference coordinates
       * \param inverseJacobian [out] : inverse jacobian
       */
      template <Core::FE::CellType distype>
      void evaluate_inverse_jacobian(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& derivs,
          Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& inverseJacobian);

      /*!
       * \brief Checks whether maerial tangent should be computed via finite difference. If so,
       * throws an error.
       *
       * \param sdyn [in] : Structural dynamics parameter list
       * \param eletype [out] : Element type string
       */
      void throw_error_fd_material_tangent(
          const Teuchos::ParameterList& sdyn, const std::string& eletype);

    }  // namespace Utils
  }    // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
