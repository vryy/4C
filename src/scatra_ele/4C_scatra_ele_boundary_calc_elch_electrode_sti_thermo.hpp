/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_STI_THERMO_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_STI_THERMO_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc_elch_electrode.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchElectrodeSTIThermo
        : public ScaTraEleBoundaryCalcElchElectrode<distype, probdim>
    {
      using my = Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelch = Discret::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>;
      using myelectrode = Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>;
      using myelectrodeutils = Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);


      /*!
       * \brief evaluate off-diagonal system matrix contributions associated with scatra-scatra
       *        interface coupling condition at integration point
       *
       * \remark This is a static method as it is also called from
       *         `scatra_timint_meshtying_strategy_s2i_elch.cpp` for the mortar implementation.
       *
       * @param[in] matelectrode     electrode material
       * @param[in] eslavephinp      scatra state variables at slave-side nodes
       * @param[in] eslavetempnp     thermo state variables at slave-side nodes
       * @param[in] emastertempnp    thermo state variables at master-side nodes
       * @param[in] emasterphinp     scatra state variables at master-side nodes
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] funct_slave      slave-side shape function values
       * @param[in] funct_master     master-side shape function values
       * @param[in] test_slave       slave-side test function values
       * @param[in] test_master      master-side test function values
       * @param[in] dsqrtdetg_dd     derivatives of the square root of the determinant of the metric
       *                             tensor w.r.t. the displacement dofs
       * @param[in] shape_spatial_derivatives  spatial derivatives of shape functions
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] differentiationtype        type of variable for linearization
       * @param[in] timefacfac       time-integration factor times domain-integration factor
       * @param[in] timefacwgt       time-integration factor times Gauss point weight
       * @param[in] detF             determinant of jacobian at current integration point
       * @param[in] num_dof_per_node number of dofs per node
       * @param[out] k_ss            linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_ms            linearizations of master-side residuals w.r.t. slave-side dofs
       *
       * \tparam distype_master  This method is templated on the master-side discretization type.
       */
      template <Core::FE::CellType distype_master>
      static void evaluate_s2_i_coupling_od_at_integration_point(
          const Teuchos::RCP<const Mat::Electrode>& matelectrode,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& eslavephinp,
          const Core::LinAlg::Matrix<nen_, 1>& eslavetempnp,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& emastertempnp,
          const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>>&
              emasterphinp,
          double pseudo_contact_fac, const Core::LinAlg::Matrix<nen_, 1>& funct_slave,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& funct_master,
          const Core::LinAlg::Matrix<nen_, 1>& test_slave,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& test_master,
          const Core::LinAlg::Matrix<nsd_, nen_>& dsqrtdetg_dd,
          const Core::LinAlg::Matrix<nsd_, nen_>& shape_spatial_derivatives,
          const Discret::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          ScaTra::DifferentiationType differentiationtype, double timefacfac, double timefacwgt,
          double detF, int num_dof_per_node, Core::LinAlg::SerialDenseMatrix& k_ss,
          Core::LinAlg::SerialDenseMatrix& k_ms);

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchElectrodeSTIThermo(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      void evaluate_s2_i_coupling_od(
          const Core::Elements::FaceElement* ele,        ///< current boundary element
          Teuchos::ParameterList& params,                ///< parameter list
          Core::FE::Discretization& discretization,      ///< discretization
          Core::Elements::Element::LocationArray& la,    ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
          ) override;

      //! evaluate action
      int evaluate_action(Core::Elements::FaceElement* ele,  //!< boundary element
          Teuchos::ParameterList& params,                    //!< parameter list
          Core::FE::Discretization& discretization,          //!< discretization
          ScaTra::BoundaryAction action,                     //!< action
          Core::Elements::Element::LocationArray& la,        //!< location array
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,   //!< element matrix 1
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,   //!< element matrix 2
          Core::LinAlg::SerialDenseVector& elevec1_epetra,   //!< element right-hand side vector 1
          Core::LinAlg::SerialDenseVector& elevec2_epetra,   //!< element right-hand side vector 2
          Core::LinAlg::SerialDenseVector& elevec3_epetra    //!< element right-hand side vector 3
          ) override;

      //! extract nodal state variables associated with boundary element
      void extract_node_values(const Core::FE::Discretization& discretization,  //!< discretization
          Core::Elements::Element::LocationArray& la                            //!< location array
          ) override;

      //! evaluate factor F/RT
      [[nodiscard]] double get_frt() const override;

      //! nodal temperature variables associated with time t_{n+1} or t_{n+alpha_f}
      Core::LinAlg::Matrix<nen_, 1> etempnp_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeSTIThermo
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
