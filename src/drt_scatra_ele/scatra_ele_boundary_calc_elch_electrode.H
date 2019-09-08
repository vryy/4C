/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes

\level 2

\maintainer Christoph Schmidt
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

      /*!
       * \brief evaluate scatra-scatra interface coupling condition at integration point
       *
       * \remark This is a static method as it is also called from
       * `scatra_timint_meshtying_strategy_s2i_elch.cpp` for the mortar implementation.
       *
       * @param[in] matelectrode     electrode material
       * @param[in] eslavephinp      state variables at slave-side nodes
       * @param[in] emasterphinp     state variables at master-side nodes
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
       * @param[in] timefacrhsfac    time-integration factor for right-hand side times
       *                             domain-integration factor
       * @param[in] frt              factor F/(RT)
       * @param[out] k_ss            linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_sm            linearizations of slave-side residuals w.r.t. master-side dofs
       * @param[out] k_ms            linearizations of master-side residuals w.r.t. slave-side dofs
       * @param[out] k_mm            linearizations of master-side residuals w.r.t. master-side
       *                             dofs
       * @param[out] r_s             slave-side residual vector
       * @param[out] r_m             master-side residual vector
       *
       * \tparam distype_master  This method is templated on the master-side discretization type.
       */
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingAtIntegrationPoint(
          const Teuchos::RCP<const MAT::Electrode>& matelectrode,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& eslavephinp,
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
          const int kineticmodel, const int numelectrons, const std::vector<int>* stoichiometries,
          const double kr, const double alphaa, const double alphac, const double timefacfac,
          const double timefacrhsfac, const double frt, Epetra_SerialDenseMatrix& k_ss,
          Epetra_SerialDenseMatrix& k_sm, Epetra_SerialDenseMatrix& k_ms,
          Epetra_SerialDenseMatrix& k_mm, Epetra_SerialDenseVector& r_s,
          Epetra_SerialDenseVector& r_m);

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
