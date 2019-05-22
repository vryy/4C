/*----------------------------------------------------------------------*/
/*!

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_H
#define SCATRA_ELE_BOUNDARY_CALC_STI_ELECTRODE_H

#include "../drt_scatra_ele/scatra_ele_boundary_calc.H"

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
    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcSTIElectrode : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcSTIElectrode<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcSTIElectrode* delete_me = NULL);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcSTIElectrode() { return; };

      //! singleton destruction
      virtual void Done();

      //! evaluate scatra-scatra interface coupling condition at integration point
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingAtIntegrationPoint(
          DRT::Condition& s2icondition,  //!< scatra-scatra interface coupling condition
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,  //!< electrode material
          const LINALG::Matrix<my::nen_, 1>&
              eslavetempnp,  //!< thermo state variables at slave-side nodes
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              eslavephinp,  //!< scatra state variables at slave-side nodes
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,  //!< scatra state variables at master-side nodes
          const LINALG::Matrix<my::nen_, 1>& funct_slave,  //!< slave-side shape function values
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,            //!< master-side shape function values
          const double timefacfac,     //!< time-integration factor times domain-integration factor
          const double timefacrhsfac,  //!< time-integration factor for right-hand side times
                                       //!< domain-integration factor
          Epetra_SerialDenseMatrix&
              k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
          Epetra_SerialDenseVector& r_s  //!< slave-side residual vector
      );

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition at integration point
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingODAtIntegrationPoint(
          DRT::Condition& s2icondition,  //!< scatra-scatra interface coupling condition
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,  //!< electrode material
          const LINALG::Matrix<my::nen_, 1>&
              eslavetempnp,  //!< thermo state variables at slave-side nodes
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              eslavephinp,  //!< scatra state variables at slave-side nodes
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,  //!< scatra state variables at master-side nodes
          const LINALG::Matrix<my::nen_, 1>& funct_slave,  //!< slave-side shape function values
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,         //!< master-side shape function values
          const double timefacfac,  //!< time-integration factor times domain-integration factor
          Epetra_SerialDenseMatrix&
              k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
          Epetra_SerialDenseMatrix&
              k_sm  //!< linearizations of slave-side residuals w.r.t. master-side dofs
      );

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
          Epetra_SerialDenseVector& eslaveresidual           ///< element residual for slave side
      );

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
      );

      //! extract nodal state variables associated with boundary element
      void ExtractNodeValues(const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la                                //!< location array
      );

      //! nodal electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
      std::vector<LINALG::Matrix<my::nen_, 1>> eelchnp_;
    };  // class ScaTraEleBoundaryCalcSTIElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
