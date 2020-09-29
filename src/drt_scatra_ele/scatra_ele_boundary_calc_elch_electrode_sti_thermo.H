/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_STI_THERMO_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_STI_THERMO_H

#include "../drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode.H"

namespace DRT
{
  // forward declaration
  class Discretization;

  namespace ELEMENTS
  {
    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcElchElectrodeSTIThermo
        : public ScaTraEleBoundaryCalcElchElectrode<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype> myelectrode;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcElchElectrodeSTIThermo* delete_me = NULL);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchElectrodeSTIThermo() { return; };

      //! singleton destruction
      virtual void Done();

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
       * @param[in] emasterphinp     scatra state variables at master-side nodes
       * @param[in] funct_slave      slave-side shape function values
       * @param[in] funct_master     master-side shape function values
       * @param[in] test_slave       slave-side test function values
       * @param[in] test_master      master-side test function values
       * @param[in] kineticmodel     kinetic model of scatra-scatra interface condition
       * @param[in] numelectrons     number of electrons involved in charge transfer
       * @param[in] stoichiometries  vector of stoichiometric coefficients for scatra-scatra
       *                             interface
       * @param[in] kr               charge transfer constant
       * @param[in] alphaa           symmetry coefficient of anodic intercalation reaction
       * @param[in] alphac           symmetry coefficient of cathodic intercalation reaction
       * @param[in] timefacfac       time-integration factor times domain-integration factor
       * @param[in] differentiationtype
       *                             type of variable for linearization
       * @param[out] k_ss            linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_ms            linearizations of master-side residuals w.r.t. slave-side dofs
       *
       * \tparam distype_master  This method is templated on the master-side discretization type.
       */
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingODAtIntegrationPoint(
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& eslavephinp,
          const LINALG::Matrix<my::nen_, 1>& eslavetempnp,
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,
          const LINALG::Matrix<my::nen_, 1>& funct_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,
          const LINALG::Matrix<my::nen_, 1>& test_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              test_master,
          int kineticmodel, int numelectrons, const std::vector<int>* stoichiometries, double kr,
          double alphaa, double alphac, double timefacfac, int differentiationtype,
          Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_ms);

     private:
      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchElectrodeSTIThermo(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      void EvaluateS2ICouplingOD(const DRT::FaceElement* ele,  ///< current boundary element
          Teuchos::ParameterList& params,                      ///< parameter list
          DRT::Discretization& discretization,                 ///< discretization
          DRT::Element::LocationArray& la,                     ///< location array
          Epetra_SerialDenseMatrix& eslavematrix               ///< element matrix for slave side
          ) override;

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

      //! evaluate factor F/RT
      double GetFRT() const override;

      //! nodal temperature variables associated with time t_{n+1} or t_{n+alpha_f}
      LINALG::Matrix<my::nen_, 1> etempnp_;
    };  // class ScaTraEleBoundaryCalcElchElectrodeSTIThermo
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
