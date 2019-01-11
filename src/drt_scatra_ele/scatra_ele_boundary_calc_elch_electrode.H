/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode.H

\brief evaluation of ScaTra boundary elements for isothermal electrodes

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_H

#include "scatra_ele_boundary_calc_elch.H"

// forward declaration
namespace MAT
{
  class Electrode;
}

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcElchElectrode : public ScaTraEleBoundaryCalcElch<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrode<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcElchElectrode* delete_me = NULL);

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchElectrode() { return; };

      //! singleton destruction
      virtual void Done();

      //! evaluate scatra-scatra interface coupling condition at integration point
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingAtIntegrationPoint(
          DRT::Condition& s2icondition,  //!< scatra-scatra interface coupling condition
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,  //!< electrode material
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              eslavephinp,  //!< state variables at slave-side nodes
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,                                //!< state variables at master-side nodes
          const LINALG::Matrix<my::nen_, 1>& funct_slave,  //!< slave-side shape function values
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,                               //!< master-side shape function values
          const LINALG::Matrix<my::nen_, 1>& test_slave,  //!< slave-side test function values
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              test_master,             //!< master-side test function values
          const double timefacfac,     //!< time-integration factor times domain-integration factor
          const double timefacrhsfac,  //!< time-integration factor for right-hand side times
                                       //!< domain-integration factor
          const double frt,            //!< factor F/(RT)
          Epetra_SerialDenseMatrix&
              k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
          Epetra_SerialDenseMatrix&
              k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
          Epetra_SerialDenseMatrix&
              k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
          Epetra_SerialDenseMatrix&
              k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
          Epetra_SerialDenseVector& r_s,  //!< slave-side residual vector
          Epetra_SerialDenseVector& r_m   //!< master-side residual vector
      );

     protected:
      //! protected constructor for singletons
      ScaTraEleBoundaryCalcElchElectrode(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate scatra-scatra interface coupling condition
      virtual void EvaluateS2ICoupling(const DRT::FaceElement* ele,  ///< current boundary element
          Teuchos::ParameterList& params,                            ///< parameter list
          DRT::Discretization& discretization,                       ///< discretization
          DRT::Element::LocationArray& la,                           ///< location array
          Epetra_SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
          Epetra_SerialDenseMatrix& emastermatrix,  ///< element matrix for master side
          Epetra_SerialDenseVector& eslaveresidual  ///< element residual for slave side
      );

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      virtual void EvaluateS2ICouplingOD(const DRT::FaceElement* ele,  ///< current boundary element
          Teuchos::ParameterList& params,                              ///< parameter list
          DRT::Discretization& discretization,                         ///< discretization
          DRT::Element::LocationArray& la,                             ///< location array
          Epetra_SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
      );

      //! extract valence of species k from element material
      virtual double GetValence(
          const Teuchos::RCP<const MAT::Material>& material,  //! element material
          const int k                                         //! species number
          ) const;

      //! evaluate factor F/RT
      virtual double GetFRT() const;
    };  // class ScaTraEleBoundaryCalcElchElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
