/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes exhibiting surface layer
growth, e.g., lithium plating

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_GROWTH_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_GROWTH_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_boundary_calc_elch_electrode.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // class implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchElectrodeGrowth
        : public ScaTraEleBoundaryCalcElchElectrode<distype, probdim>
    {
      using my = Discret::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelch = Discret::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>;
      using myelectrode = Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);


     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchElectrodeGrowth(
          const int numdofpernode, const int numscal, const std::string& disname);

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

      //! evaluate minimum and maximum interfacial overpotential associated with scatra-scatra
      //! interface layer growth
      void evaluate_min_max_overpotential(
          const Core::Elements::FaceElement* ele,     //!< current boundary element
          Teuchos::ParameterList& params,             //!< parameter list
          Core::FE::Discretization& discretization,   //!< discretization
          Core::Elements::Element::LocationArray& la  //!< location array
      );

      /*!
       * \brief evaluate scatra-scatra interface coupling condition
       *
       * @param[in] ele              current boundary element
       * @param[in] params           parameter list
       * @param[in] discretization   discretization
       * @param[in] la               location array
       * @param[out] eslavematrix    element matrix for slave side
       * @param[out] emastermatrix   element matrix for master side
       * @param[out] eslaveresidual  element residual for slave side
       */
      void evaluate_s2_i_coupling(const Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& eslavematrix,
          Core::LinAlg::SerialDenseMatrix& emastermatrix,
          Core::LinAlg::SerialDenseVector& eslaveresidual) override;

      //! evaluate global growth-growth matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void evaluate_s2_i_coupling_growth_growth(
          const Core::Elements::FaceElement* ele,          ///< current boundary element
          Teuchos::ParameterList& params,                  ///< parameter list
          Core::FE::Discretization& discretization,        ///< discretization
          Core::Elements::Element::LocationArray& la,      ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
          Core::LinAlg::SerialDenseVector& eslaveresidual  ///< element residual for slave side
      );

      //! evaluate global growth-scatra matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void evaluate_s2_i_coupling_growth_scatra(
          const Core::Elements::FaceElement* ele,         ///< current boundary element
          Teuchos::ParameterList& params,                 ///< parameter list
          Core::FE::Discretization& discretization,       ///< discretization
          Core::Elements::Element::LocationArray& la,     ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix,  ///< element matrix for slave side
          Core::LinAlg::SerialDenseMatrix& emastermatrix  ///< element matrix for master side
      );

      //! evaluate global scatra-growth matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void evaluate_s2_i_coupling_scatra_growth(
          const Core::Elements::FaceElement* ele,        ///< current boundary element
          Teuchos::ParameterList& params,                ///< parameter list
          Core::FE::Discretization& discretization,      ///< discretization
          Core::Elements::Element::LocationArray& la,    ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
      );

      //! extract nodal state variables associated with boundary element
      void extract_node_values(const Core::FE::Discretization& discretization,  //!< discretization
          Core::Elements::Element::LocationArray& la                            //!< location array
          ) override;

      /*!
       * @brief Fill the rhs vector and the linearization matrix
       *
       * @param[in] numelectrons     number of electrons involved in charge transfer at
       *                             electrode-electrolyte interface
       * @param[in] timefacfac       time and domain integration factor for linearization terms
       * @param[in] timefacrhsfac    time and domain integration factor for RHS terms
       * @param[in] j                Butler-Volmer mass flux density
       * @param[in] dj_dc_slave      linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on slave-side
       * @param[in] dj_dc_master     linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on master-side
       * @param[in] dj_dpot_slave    linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on slave-side
       * @param[in] dj_dpot_master   linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on master-side
       * @param[out] eslavematrix    linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] emastermatrix   linearizations of slave-side residuals w.r.t. master-side dofs
       * @param[out] eslaveresidual  slave-side residual vector
       */
      void calculate_rhs_and_linearization(int numelectrons, double timefacfac,
          double timefacrhsfac, double j, double dj_dc_slave, double dj_dc_master,
          double dj_dpot_slave, double dj_dpot_master,
          Core::LinAlg::SerialDenseMatrix& eslavematrix,
          Core::LinAlg::SerialDenseMatrix& emastermatrix,
          Core::LinAlg::SerialDenseVector& eslaveresidual) const;

      //! nodal growth variables
      Core::LinAlg::Matrix<nen_, 1> egrowth_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeGrowth
  }     // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
