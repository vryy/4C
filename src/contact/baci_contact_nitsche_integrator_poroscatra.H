/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the poroscatra contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_INTEGRATOR_POROSCATRA_H
#define BACI_CONTACT_NITSCHE_INTEGRATOR_POROSCATRA_H

#include "baci_config.H"

#include "baci_contact_nitsche_integrator_poro.H"
#include "baci_utils_pairedvector.H"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  class IntegratorNitschePoroscatra : public IntegratorNitschePoro
  {
   public:
    using BaseType = IntegratorNitschePoro;

    /*!
     \brief Constructor  with shape function specification

     Constructs an instance of this class using a specific type of shape functions.<br>
     Note that this is \b not a collective call as overlaps are
     integrated in parallel by individual processes.<br>
     Note also that this constructor relies heavily on the
     CORE::DRT::UTILS::IntegrationPoints structs to get Gauss points
     and corresponding weights.

     */
    IntegratorNitschePoroscatra(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm);

   protected:
    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    void IntegrateGP_3D(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi) override;

    // /*!
    // \brief Evaluate cauchy stress component and its derivatives
    // */
    // template <int dim>
    // void SoEleCauchy(MORTAR::Element& moEle, double* boundary_gpcoord,
    //     std::vector<CORE::GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double
    //     gp_wgt, const CORE::LINALG::Matrix<dim, 1>& normal,
    //     std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
    //     const CORE::LINALG::Matrix<dim, 1>& direction,
    //     std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv, const double w,
    //     double& cauchy_nt, CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_d,
    //     CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_p,
    //     CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_s);

   protected:
    /*!
    \brief evaluate GPTS forces and linearization at this gp
    */
    template <int dim>
    void GPTSForces(MORTAR::Element& sele, MORTAR::Element& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxi,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double gap, const CORE::GEN::pairedvector<int, double>& dgapgp, const double* gpn,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi);

   private:
    /*!
     * @brief calculate the concentrations and derivatives at the current Gauss point
     *
     * @tparam dim  dimension of the problem
     * @param[in] ele   mortar contact element or integration cell mortar element
     * @param[in] shape_func   shape function evaluated at current Gauss point
     * @param[in] shape_deriv  shape function derivative at current Gauss point
     * @param[in] d_xi_dd     directional derivative of Gauss point coordinates
     * @param[out] gp_conc    concentration at current Gauss point
     * @param[out] d_conc_dc  derivative of concentration w.r.t. concentration
     * @param[out] d_conc_dd  directional derivative of concentration
     */
    template <int dim>
    void SetupGpConcentrations(MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape_func,
        const CORE::LINALG::SerialDenseMatrix& shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_xi_dd, double& gp_conc,
        CORE::GEN::pairedvector<int, double>& d_conc_dc,
        CORE::GEN::pairedvector<int, double>& d_conc_dd);

    /*!
     * @brief  integrate the ScaTra residual and linearizations
     *
     * @tparam dim  dimension of the problem
     * @param[in] fac  pre-factor to correct sign dependent on integration of master or slave side
     *                 terms
     * @param[in,out] ele      mortar contact element or integration cell mortar element
     * @param[in] shape_func   shape function evaluated at current Gauss point
     * @param[in] shape_deriv  shape function derivative at current Gauss point
     * @param[in] d_xi_dd      directional derivative of Gauss point coordinates
     * @param[in] jac          Jacobian determinant of integration cell
     * @param[in] d_jac_dd     directional derivative of cell Jacobian
     * @param[in] wgt          Gauss point weight
     * @param[in] test_val     quantity to be integrated
     * @param[in] d_test_val_dd  directional derivative of test_val
     * @param[in] d_test_val_ds  derivative of test_val w.r.t. scalar
     */
    template <int dim>
    void IntegrateScaTraTest(double fac, MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape_func,
        const CORE::LINALG::SerialDenseMatrix& shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_xi_dd, double jac,
        const CORE::GEN::pairedvector<int, double>& d_jac_dd, double wgt, double test_val,
        const CORE::GEN::pairedvector<int, double>& d_test_val_dd,
        const CORE::GEN::pairedvector<int, double>& d_test_val_ds,
        const CORE::GEN::pairedvector<int, double>& d_test_val_dp);

    /*!
     * @brief integrate the scatra-structure interaction interface condition
     *
     * @tparam dim  dimension of the problem
     * @param[in,out] slave_ele       slave-side mortar contact element
     * @param[in] slave_shape         slave-side shape function evaluated at current Gauss point
     * @param[in] slave_shape_deriv   slave-side shape function derivative at current Gauss point
     * @param[in] d_slave_xi_dd       slave-side directional derivative of Gauss point coordinates
     * @param[in,out] master_ele      master-side mortar contact element
     * @param[in] master_shape        master-side shape function evaluated at current Gauss point
     * @param[in] master_shape_deriv  master-side shape function derivative at current Gauss point
     * @param[in] d_master_xi_dd      master-side directional derivative of Gauss point coordinates
     * @param[in] jac                 Jacobian determinant of integration cell
     * @param[in] d_jac_dd            directional derivative of cell Jacobian
     * @param[in] wgt                 Gauss point weight
     */
    template <int dim>
    void IntegrateScatraInterface(MORTAR::Element& slave_ele,
        const CORE::LINALG::SerialDenseVector& slave_shape,
        const CORE::LINALG::SerialDenseMatrix& slave_shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_slave_xi_dd,
        MORTAR::Element& master_ele, const CORE::LINALG::SerialDenseVector& master_shape,
        const CORE::LINALG::SerialDenseMatrix& master_shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_master_xi_dd, double jac,
        const CORE::GEN::pairedvector<int, double>& d_jac_dd, double wgt, const double* gpn,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi,
        double contact_stress, const CORE::GEN::pairedvector<int, double>& contact_stress_deriv_d,
        const CORE::GEN::pairedvector<int, double>& contact_stress_deriv_p);
  };
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // BACI_CONTACT_NITSCHE_INTEGRATOR_POROSCATRA_H
