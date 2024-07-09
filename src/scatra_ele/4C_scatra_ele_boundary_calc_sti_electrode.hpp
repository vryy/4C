/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Mat
{
  class Electrode;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class ScaTraEleBoundaryCalcElchElectrodeUtils;

    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcSTIElectrode : public ScaTraEleBoundaryCalc<distype, probdim>
    {
      using my = Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelectrodeutils = Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /*!
       * \brief evaluate scatra-scatra interface coupling condition at integration point
       *
       * \remark This is a static method as it is also called from
       * `scatra_timint_meshtying_strategy_s2i_elch.cpp` for the mortar implementation.
       *
       * @param[in] matelectrode   electrode material
       * @param[in] eslavetempnp   thermo state variables at slave-side nodes
       * @param[in] emastertempnp  thermo state variables at master-side nodes
       * @param[in] eslavephinp    scatra state variables at slave-side nodes
       * @param[in] emasterphinp   scatra state variables at master-side nodes
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] funct_slave                slave-side shape function values
       * @param[in] funct_master               master-side shape function values
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] timefacfac     time-integration factor times domain-integration factor
       * @param[in] timefacrhsfac  time-integration factor for right-hand side times
       *                           domain-integration factor
       * @param[in] detF           determinant of jacobian at current integration point
       * @param[out] k_ss          linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_sm          linearizations of slave-side residuals w.r.t. master-side dofs
       * @param[out] r_s           slave-side residual vector
       *
       * \tparam distype_master  This method is templated on the master-side discretization type.
       */
      template <Core::FE::CellType distype_master>
      static void evaluate_s2_i_coupling_at_integration_point(
          const Teuchos::RCP<const Mat::Electrode>& matelectrode,
          const Core::LinAlg::Matrix<nen_, 1>& eslavetempnp,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& emastertempnp,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& eslavephinp,
          const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>>&
              emasterphinp,
          double pseudo_contact_fac, const Core::LinAlg::Matrix<nen_, 1>& funct_slave,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& funct_master,
          const Discret::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          double timefacfac, double timefacrhsfac, double detF,
          Core::LinAlg::SerialDenseMatrix& k_ss, Core::LinAlg::SerialDenseMatrix& k_sm,
          Core::LinAlg::SerialDenseVector& r_s);

      /*!
       * \brief evaluate off-diagonal system matrix contributions associated with scatra-scatra
       *        interface coupling condition at integration point
       * @param[in] matelectrode  electrode material
       * @param[in] eslavetempnp  thermo state variables at slave-side nodes
       * @param[in] emastertempnp thermo state variables at master-side nodes
       * @param[in] eslavephinp   scatra state variables at slave-side nodes
       * @param[in] emasterphinp  scatra state variables at master-side nodes
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] funct_slave                slave-side shape function values
       * @param[in] funct_master               master-side shape function values
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] timefacfac      time-integration factor times domain-integration factor
       * @param[in] timefacwgt      time-integration factor times Gauss point weight
       * @param[in] detF            determinant of Jacobian from deformation at Gauss point
       * @param[in] differentiationtype        type of variable for linearization
       * @param[in] dsqrtdetg_dd               derivatives of the square root of the determinant of
       *                                       the metric tensor w.r.t. the displacement dofs
       * @param[in] shape_spatial_derivatives  spatial derivative of shape functions
       * @param[out] k_ss         linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_sm         linearizations of slave-side residuals w.r.t. master-side dofs
       */
      template <Core::FE::CellType distype_master>
      static void evaluate_s2_i_coupling_od_at_integration_point(
          const Teuchos::RCP<const Mat::Electrode>& matelectrode,
          const Core::LinAlg::Matrix<nen_, 1>& eslavetempnp,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& emastertempnp,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& eslavephinp,
          const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>>&
              emasterphinp,
          double pseudo_contact_fac, const Core::LinAlg::Matrix<nen_, 1>& funct_slave,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& funct_master,
          const Discret::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          double timefacfac, double timefacwgt, double detF,
          ScaTra::DifferentiationType differentiationtype,
          const Core::LinAlg::Matrix<nsd_, nen_>& dsqrtdetg_dd,
          const Core::LinAlg::Matrix<nsd_, nen_>& shape_spatial_derivatives,
          Core::LinAlg::SerialDenseMatrix& k_ss, Core::LinAlg::SerialDenseMatrix& k_sm);

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcSTIElectrode(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate main-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      void evaluate_s2_i_coupling(
          const Core::Elements::FaceElement* ele,          ///< current boundary element
          Teuchos::ParameterList& params,                  ///< parameter list
          Core::FE::Discretization& discretization,        ///< discretization
          Core::Elements::Element::LocationArray& la,      ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
          Core::LinAlg::SerialDenseMatrix& emastermatrix,  ///< element matrix for master side
          Core::LinAlg::SerialDenseVector& eslaveresidual  ///< element residual for slave side
          ) override;

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      void evaluate_s2_i_coupling_od(
          const Core::Elements::FaceElement* ele,         ///< current boundary element
          Teuchos::ParameterList& params,                 ///< parameter list
          Core::FE::Discretization& discretization,       ///< discretization
          Core::Elements::Element::LocationArray& la,     ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix,  ///< element matrix for slave side
          Core::LinAlg::SerialDenseMatrix& emastermatrix  ///< element matrix for master side
      );

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

      //! nodal electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
      std::vector<Core::LinAlg::Matrix<nen_, 1>> eelchnp_;
    };  // class ScaTraEleBoundaryCalcSTIElectrode
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
