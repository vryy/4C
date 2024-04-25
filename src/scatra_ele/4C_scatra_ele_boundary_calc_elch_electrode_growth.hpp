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

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalcElchElectrodeGrowth
        : public ScaTraEleBoundaryCalcElchElectrode<distype, probdim>
    {
      using my = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelch = DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>;
      using myelectrode = DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>;
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
      int EvaluateAction(DRT::FaceElement* ele,             //!< boundary element
          Teuchos::ParameterList& params,                   //!< parameter list
          DRT::Discretization& discretization,              //!< discretization
          SCATRA::BoundaryAction action,                    //!< action
          DRT::Element::LocationArray& la,                  //!< location array
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
          ) override;

      //! evaluate minimum and maximum interfacial overpotential associated with scatra-scatra
      //! interface layer growth
      void EvaluateMinMaxOverpotential(const DRT::FaceElement* ele,  //!< current boundary element
          Teuchos::ParameterList& params,                            //!< parameter list
          DRT::Discretization& discretization,                       //!< discretization
          DRT::Element::LocationArray& la                            //!< location array
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
      void EvaluateS2ICoupling(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& eslavematrix,
          CORE::LINALG::SerialDenseMatrix& emastermatrix,
          CORE::LINALG::SerialDenseVector& eslaveresidual) override;

      //! evaluate global growth-growth matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void EvaluateS2ICouplingGrowthGrowth(
          const DRT::FaceElement* ele,                     ///< current boundary element
          Teuchos::ParameterList& params,                  ///< parameter list
          DRT::Discretization& discretization,             ///< discretization
          DRT::Element::LocationArray& la,                 ///< location array
          CORE::LINALG::SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
          CORE::LINALG::SerialDenseVector& eslaveresidual  ///< element residual for slave side
      );

      //! evaluate global growth-scatra matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void EvaluateS2ICouplingGrowthScatra(
          const DRT::FaceElement* ele,                    ///< current boundary element
          Teuchos::ParameterList& params,                 ///< parameter list
          DRT::Discretization& discretization,            ///< discretization
          DRT::Element::LocationArray& la,                ///< location array
          CORE::LINALG::SerialDenseMatrix& eslavematrix,  ///< element matrix for slave side
          CORE::LINALG::SerialDenseMatrix& emastermatrix  ///< element matrix for master side
      );

      //! evaluate global scatra-growth matrix block for scatra-scatra interface coupling involving
      //! interface layer growth
      void EvaluateS2ICouplingScatraGrowth(
          const DRT::FaceElement* ele,                   ///< current boundary element
          Teuchos::ParameterList& params,                ///< parameter list
          DRT::Discretization& discretization,           ///< discretization
          DRT::Element::LocationArray& la,               ///< location array
          CORE::LINALG::SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
      );

      //! extract nodal state variables associated with boundary element
      void ExtractNodeValues(const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la                                //!< location array
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
      void CalculateRHSAndLinearization(int numelectrons, double timefacfac, double timefacrhsfac,
          double j, double dj_dc_slave, double dj_dc_master, double dj_dpot_slave,
          double dj_dpot_master, CORE::LINALG::SerialDenseMatrix& eslavematrix,
          CORE::LINALG::SerialDenseMatrix& emastermatrix,
          CORE::LINALG::SerialDenseVector& eslaveresidual) const;

      //! nodal growth variables
      CORE::LINALG::Matrix<nen_, 1> egrowth_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeGrowth
  }     // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
