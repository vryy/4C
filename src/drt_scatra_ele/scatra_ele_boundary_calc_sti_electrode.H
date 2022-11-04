/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_H
#define SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_H

#include "scatra_ele_boundary_calc.H"

// forward declaration
namespace MAT
{
  class Electrode;
}

namespace DRT
{
  // forward declaration
  class Discretization;

  namespace ELEMENTS
  {
    class ScaTraEleBoundaryCalcElchElectrodeUtils;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcSTIElectrode : public ScaTraEleBoundaryCalc<distype>
    {
      using my = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>;
      using myelectrodeutils = DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcSTIElectrode<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcSTIElectrode() { return; };



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
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingAtIntegrationPoint(
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,
          const LINALG::Matrix<my::nen_, 1>& eslavetempnp,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              emastertempnp,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& eslavephinp,
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,
          const LINALG::Matrix<my::nen_, 1>& funct_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          const double timefacfac, const double timefacrhsfac, const double detF,
          Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm,
          Epetra_SerialDenseVector& r_s);

      /*!
       * \brief evaluate off-diagonal system matrix contributions associated with scatra-scatra
       *        interface coupling condition at integration point
       * @param[in] matelectrode  electrode material
       * @param[in] eslavetempnp  thermo state variables at slave-side nodes
       * @param[in] emastertempnp thermo state variables at master-side nodes
       * @param[in] eslavephinp   scatra state variables at slave-side nodes
       * @param[in] emasterphinp  scatra state variables at master-side nodes
       * @param[in] funct_slave                slave-side shape function values
       * @param[in] funct_master               master-side shape function values
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] timefacfac      time-integration factor times domain-integration factor
       * @param[in] timefacwgt      time-integration factor times Gauss point weight
       * @param[in] detF            determinant of Jacobian from deformation at Gauss point
       * @param[in] differentiationtype
       *                          type of variable for linearization
       * @param[i] shapederivatives
       *                          derivative of shape functions
       * @param[out] k_ss         linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_sm         linearizations of slave-side residuals w.r.t. master-side dofs
       */
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingODAtIntegrationPoint(
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,
          const LINALG::Matrix<my::nen_, 1>& eslavetempnp,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              emastertempnp,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& eslavephinp,
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,
          const LINALG::Matrix<my::nen_, 1>& funct_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          double timefacfac, double timefacwgt, double detF,
          SCATRA::DifferentiationType differentiationtype,
          const LINALG::Matrix<my::nsd_ + 1, my::nen_>& shapederivatives,
          Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm);

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcSTIElectrode(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate main-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      void EvaluateS2ICoupling(const DRT::FaceElement* ele,  ///< current boundary element
          Teuchos::ParameterList& params,                    ///< parameter list
          DRT::Discretization& discretization,               ///< discretization
          DRT::Element::LocationArray& la,                   ///< location array
          Epetra_SerialDenseMatrix& eslavematrix,            ///< element matrix for slave side
          Epetra_SerialDenseMatrix& emastermatrix,           ///< element matrix for master side
          Epetra_SerialDenseVector& eslaveresidual           ///< element residual for slave side
          ) override;

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      void EvaluateS2ICouplingOD(const DRT::FaceElement* ele,  ///< current boundary element
          Teuchos::ParameterList& params,                      ///< parameter list
          DRT::Discretization& discretization,                 ///< discretization
          DRT::Element::LocationArray& la,                     ///< location array
          Epetra_SerialDenseMatrix& eslavematrix,              ///< element matrix for slave side
          Epetra_SerialDenseMatrix& emastermatrix              ///< element matrix for master side
      );

      //! evaluate action
      int EvaluateAction(DRT::FaceElement* ele,      //!< boundary element
          Teuchos::ParameterList& params,            //!< parameter list
          DRT::Discretization& discretization,       //!< discretization
          SCATRA::BoundaryAction action,             //!< action
          DRT::Element::LocationArray& la,           //!< location array
          Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
          Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
          Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
          Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
          ) override;

      //! extract nodal state variables associated with boundary element
      void ExtractNodeValues(const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la                                //!< location array
          ) override;

      //! nodal electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
      std::vector<LINALG::Matrix<my::nen_, 1>> eelchnp_;
    };  // class ScaTraEleBoundaryCalcSTIElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
