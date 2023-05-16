/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for isothermal electrodes

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_H

#include "scatra_ele_boundary_calc_elch.H"

// forward declarations
namespace MAT
{
  class Electrode;
}
namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleBoundaryCalcElchElectrodeUtils;
  }
}  // namespace DRT

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <DRT::Element::DiscretizationType distype,
        int probdim = DRT::UTILS::DisTypeToDim<distype>::dim + 1>
    class ScaTraEleBoundaryCalcElchElectrode : public ScaTraEleBoundaryCalcElch<distype, probdim>
    {
      using my = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>;
      using myelch = DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>;
      using myelectrodeutils = DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils;

     protected:
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! singleton access method
      static ScaTraEleBoundaryCalcElchElectrode<distype, probdim>* Instance(
          int numdofpernode, int numscal, const std::string& disname);


      /*!
       * \brief evaluate scatra-scatra interface coupling condition at integration point
       *
       * \remark This is a static method as it is also called from
       * `scatra_timint_meshtying_strategy_s2i_elch.cpp` for the mortar implementation.
       *
       * @param[in] matelectrode     electrode material
       * @param[in] eslavephinp      state variables at slave-side nodes
       * @param[in] emasterphinp     state variables at master-side nodes
       * @param[in] eslavetempnp     temperature at slave side of interface
       * @param[in] emastertempnp    temperature at master side of interface
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] funct_slave      slave-side shape function values
       * @param[in] funct_master     master-side shape function values
       * @param[in] test_slave                 slave-side test function values
       * @param[in] test_master                master-side test function values
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] timefacfac       time-integration factor times domain-integration factor
       * @param[in] timefacrhsfac    time-integration factor for right-hand side times
       *                             domain-integration factor
       * @param[in] detF             determinant of jacobian at current integration point
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
          const std::vector<LINALG::Matrix<nen_, 1>>& eslavephinp,
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,
          const LINALG::Matrix<nen_, 1>& eslavetempnp,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              emastertempnp,
          double pseudo_contact_fac, const LINALG::Matrix<nen_, 1>& funct_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,
          const LINALG::Matrix<nen_, 1>& test_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              test_master,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          double timefacfac, double timefacrhsfac, double detF, double frt,
          Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_sm,
          Epetra_SerialDenseMatrix& k_ms, Epetra_SerialDenseMatrix& k_mm,
          Epetra_SerialDenseVector& r_s, Epetra_SerialDenseVector& r_m);

      /*!
       * @brief evaluate capacitive part of the scatra-scatra interface coupling condition at
       * integration point
       *
       * @param[in] eslavephidtnp  time derivative of state variables at slave-side nodes at time
       *                           n+1
       * @param[in] emasterphidtnp time derivative of state variables at master-side nodes at time
       *                           n+1
       * @param[in] eslavephinp    state variables at slave-side nodes at time n+1
       * @param[in] emasterphinp   state variables at master-side nodes at time n+1
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] funct_slave    slave-side shape function values
       * @param[in] funct_master   master-side shape function values
       * @param[in] test_slave                 slave-side test function values
       * @param[in] test_master                master-side test function values
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] timederivfac   factor to account for time-derivative in linearization
       * @param[in] timefacfac     time-integration factor times domain-integration factor
       * @param[in] timefacrhsfac  time-integration factor for right-hand side times
       *                           domain-integration factor
       * @param[out] k_ss  linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_ms  linearizations of master-side residuals w.r.t. slave-side dofs
       * @param[out] r_s   slave-side residual vector
       * @param[out] r_m   master-side residual vector
       *
       * @tparam distype_master This method is templated on the master-side discretization type.
       */
      template <DRT::Element::DiscretizationType distype_master>
      static void EvaluateS2ICouplingCapacitanceAtIntegrationPoint(
          const std::vector<LINALG::Matrix<nen_, 1>>& eslavephidtnp,
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphidtnp,
          const std::vector<LINALG::Matrix<nen_, 1>>& eslavephinp,
          const std::vector<LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>>&
              emasterphinp,
          double pseudo_contact_fac, const LINALG::Matrix<nen_, 1>& funct_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,
          const LINALG::Matrix<nen_, 1>& test_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              test_master,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          double timederivfac, double timefacfac, double timefacrhsfac,
          Epetra_SerialDenseMatrix& k_ss, Epetra_SerialDenseMatrix& k_ms,
          Epetra_SerialDenseVector& r_s, Epetra_SerialDenseVector& r_m);

      /*!
       * \brief calculate out-parameters such as residual vectors and linearizations of residuals
       *
       * \remark This is a static method as it is called from
       * static method EvaluateS2ICouplingAtIntegrationPoint.
       *
       * @param[in] funct_slave      slave-side shape function values
       * @param[in] funct_master     master-side shape function values
       * @param[in] test_slave       slave-side test function values
       * @param[in] test_master      master-side test function values
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] numelectrons     number of electrons involved in charge transfer at
       *                             electrode-electrolyte interface
       * @param[in] nen_master       number of nodes of master-side mortar element
       * @param[in] timefacfac       time and domain integration factor for linearization terms
       * @param[in] timefacrhsfac    time and domain integration factor for RHS terms
       * @param[in] dj_dc_slave      linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on slave-side
       * @param[in] dj_dc_master     linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on master-side
       * @param[in] dj_dpot_slave    linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on slave-side
       * @param[in] dj_dpot_master   linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on master-side
       * @param[in] j                Butler-Volmer mass flux density
       * @param[out] k_ss            linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_sm            linearizations of slave-side residuals w.r.t. master-side dofs
       * @param[out] k_ms            linearizations of master-side residuals w.r.t. slave-side dofs
       * @param[out] k_mm            linearizations of master-side residuals w.r.t. master-side dofs
       * @param[out] r_s             slave-side residual vector
       * @param[out] r_m             master-side residual vector
       *
       * @tparam distype_master  This method is templated on the master-side discretization type.
       */
      template <DRT::Element::DiscretizationType distype_master>
      static void CalculateRHSandGlobalSystem(const LINALG::Matrix<nen_, 1>& funct_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              funct_master,
          const LINALG::Matrix<nen_, 1>& test_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              test_master,
          double pseudo_contact_fac, double numelectrons, int nen_master, double timefacfac,
          double timefacrhsfac, double dj_dc_slave, double dj_dc_master, double dj_dpot_slave,
          double dj_dpot_master, double j, Epetra_SerialDenseMatrix& k_ss,
          Epetra_SerialDenseMatrix& k_sm, Epetra_SerialDenseMatrix& k_ms,
          Epetra_SerialDenseMatrix& k_mm, Epetra_SerialDenseVector& r_s,
          Epetra_SerialDenseVector& r_m);

      /*!
       * @brief calculate residual vectors and linearizations of residuals due to capacitive flux at
       * the scatra-scatra interface
       *
       * @param[in] funct_slave     slave-side shape function values
       * @param[in] test_slave      slave-side test function values
       * @param[in] test_master     master-side test function values
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] numelectrons    number of electrons involved in charge transfer
       * @param[in] timefacfac      time-integration factor times domain-integration factor
       * @param[in] timefacrhsfac   time-integration factor for right-hand side times
       *                            domain-integration factor
       * @param[in] nen_master      number of nodes of master-side mortar element
       * @param[in] jC              capacitance mass flux density
       * @param[in] djC_dpot_slave  derivative of capacitance mass flux density w.r.t. electric
       *                            potential at the slave-side
       * @param[out] k_ss  linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_ms  linearizations of master-side residuals w.r.t. slave-side dofs
       * @param[out] r_s   slave-side residual vector
       * @param[out] r_m   master-side residual vector
       *
       * @tparam distype_master  This method is templated on the master-side discretization type.
       */
      template <DRT::Element::DiscretizationType distype_master>
      static void CalculateRHSandGlobalSystemCapacitiveFlux(
          const LINALG::Matrix<nen_, 1>& funct_slave, const LINALG::Matrix<nen_, 1>& test_slave,
          const LINALG::Matrix<
              DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement, 1>&
              test_master,
          double pseudo_contact_fac, int numelectrons, double timefacfac, double timefacrhsfac,
          int nen_master, double jC, double djC_dpot_slave, Epetra_SerialDenseMatrix& k_ss,
          Epetra_SerialDenseMatrix& k_ms, Epetra_SerialDenseVector& r_s,
          Epetra_SerialDenseVector& r_m);

     protected:
      //! protected constructor for singletons
      ScaTraEleBoundaryCalcElchElectrode(
          int numdofpernode, int numscal, const std::string& disname);

      void EvaluateS2ICoupling(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix,
          Epetra_SerialDenseVector& eslaveresidual) override;

      void EvaluateS2ICouplingCapacitance(const DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& eslavematrix,
          Epetra_SerialDenseMatrix& emastermatrix, Epetra_SerialDenseVector& eslaveresidual,
          Epetra_SerialDenseVector& emasterresidual) override;

      void EvaluateS2ICouplingOD(const DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& eslavematrix) override;

      void EvaluateS2ICouplingCapacitanceOD(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& eslavematrix, Epetra_SerialDenseMatrix& emastermatrix) override;

      double GetValence(const Teuchos::RCP<const MAT::Material>& material, int k) const override;

      void CalcS2ICouplingFlux(const DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseVector& scalars) override;

      //! evaluate factor F/RT
      virtual double GetFRT() const;
    };  // class ScaTraEleBoundaryCalcElchElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
