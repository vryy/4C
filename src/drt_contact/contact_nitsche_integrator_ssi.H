/*----------------------------------------------------------------------------*/
/*! \file
\brief A class to perform integration of nitsche related terms for the ssi contact case

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_INTEGRATOR_SSI_H
#define CONTACT_NITSCHE_INTEGRATOR_SSI_H

#include "contact_nitsche_integrator.H"

// forward declarations
namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterTimInt;
    class ScaTraEleParameterBoundary;
  }  // namespace ELEMENTS
}  // namespace DRT

namespace CONTACT
{
  class CoIntegratorNitscheSsi : public CoIntegratorNitsche
  {
   public:
    /*!
     * \brief Constructor with shape function specification
     *
     * Constructs an instance of this class using a specific type of shape functions.<br> Note that
     * this is \b not a collective call as overlaps are integrated in parallel by individual
     * processes.<br> Note also that this constructor relies heavily on the
     * DRT::UTILS::IntegrationPoints structs to get Gauss points and corresponding weights.
     *
     * @param[in] params   interface contact parameter list
     * @param[in] eletype  shape of integration cell for segment based integration or slave side
     *                     mortar contact element for element based integration
     * @param[in] comm     contact interface communicator
     */
    CoIntegratorNitscheSsi(Teuchos::ParameterList& params, DRT::Element::DiscretizationType eletype,
        const Epetra_Comm& comm);

   private:
    //! access method to scatra time integration factors
    const DRT::ELEMENTS::ScaTraEleParameterTimInt* GetScaTraEleParameterTimInt() const
    {
      return scatraparamstimint_;
    }

    //! access method to scatra-scatra coupling specific parameters
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* GetScaTraEleParameterBoundary() const
    {
      return scatraparamsboundary_;
    }

    void IntegrateGP_2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
        GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
        GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<GEN::pairedvector<int, double>>& derivsxi,
        std::vector<GEN::pairedvector<int, double>>& derivmxi) override;

    void IntegrateGP_3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
        GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
        GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<GEN::pairedvector<int, double>>& derivsxi,
        std::vector<GEN::pairedvector<int, double>>& derivmxi) override;

    /*!
     * @brief evaluate GPTS forces and linearization at this gp
     *
     * @tparam dim  dimension of the problem
     * @param[in,out] slave_ele       slave side mortar element
     * @param[in,out] master_ele      master side mortar element
     * @param[in] slave_shape         slave side shape function evaluated at current Gauss point
     * @param[in] slave_shape_deriv   slave side shape function derivative at current Gauss point
     * @param[in] d_slave_xi_dd       directional derivative of slave side Gauss point coordinates
     * @param[in] master_shape        master side shape function evaluated at current Gauss point
     * @param[in] master_shape_deriv  master side shape function derivative at current Gauss point
     * @param[in] d_master_xi_dd  directional derivative of master side Gauss point coordinates
     * @param[in] jac             Jacobian determinant of integration cell
     * @param[in] d_jac_dd        directional derivative of cell Jacobian
     * @param[in] gp_wgt          Gauss point weight
     * @param[in] gap             gap
     * @param[in] d_gap_dd        directional derivative of gap
     * @param[in] gp_normal       Gauss point normal
     * @param[in] d_gp_normal_dd  directional derivative of Gauss point normal
     * @param[in] slave_xi        slave side Gauss point coordinates
     * @param[in] master_xi       master side Gauss point coordinates
     */
    template <int dim>
    void GPTSForces(MORTAR::MortarElement& slave_ele, MORTAR::MortarElement& master_ele,
        const LINALG::SerialDenseVector& slave_shape,
        const LINALG::SerialDenseMatrix& slave_shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_slave_xi_dd,
        const LINALG::SerialDenseVector& master_shape,
        const LINALG::SerialDenseMatrix& master_shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_master_xi_dd, double jac,
        const GEN::pairedvector<int, double>& d_jac_dd, double gp_wgt, double gap,
        const GEN::pairedvector<int, double>& d_gap_dd, const double* gp_normal,
        const std::vector<GEN::pairedvector<int, double>>& d_gp_normal_dd, double* slave_xi,
        double* master_xi);

    /*!
     * @brief  Evaluate cauchy stress component and its derivatives
     *
     * @tparam dim  dimension of the problem
     * @param[in] mortar_ele     mortar element
     * @param[in] gp_coord       Gauss point coordinates
     * @param[in] d_gp_coord_dd  directional derivative of Gauss point coordinates
     * @param[in] gp_wgt         Gauss point weight
     * @param[in] gp_normal      Gauss point normal
     * @param[in] d_gp_normal_dd directional derivative of Gauss point normal
     * @param[in] test_dir       direction of evaluation (e.g. normal or tangential direction)
     * @param[in] d_test_dir_dd  directional derivative of direction of evaluation
     * @param[in] nitsche_wgt    Nitsche weight
     * @param[out] cauchy_nt_wgt   Cauchy stress tensor contracted with normal vector n and
     *                             direction vector t multiplied by nitsche_wgt \f[ nitsche_wgt
     *                             \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot \boldsymbol{t} \f]
     * @param[out] d_cauchy_nt_dd  directional derivative of cauchy_nt \f[ \frac{ \mathrm{d}
     *                             \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                             \boldsymbol{t}}{\mathrm{d} \boldsymbol{d}} \f]
     * @param[out] d_cauchy_nt_ds  derivative of cauchy_nt w.r.t. scalar s \f[ \frac{ \mathrm{d}
     *                             \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                             \boldsymbol{t}}{\mathrm{d} s} \f]
     */
    template <int dim>
    void SoEleCauchy(MORTAR::MortarElement& mortar_ele, double* gp_coord,
        const std::vector<GEN::pairedvector<int, double>>& d_gp_coord_dd, double gp_wgt,
        const LINALG::Matrix<dim, 1>& gp_normal,
        const std::vector<GEN::pairedvector<int, double>>& d_gp_normal_dd,
        const LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<GEN::pairedvector<int, double>>& d_test_dir_dd, double nitsche_wgt,
        double& cauchy_nt_wgt, GEN::pairedvector<int, double>& d_cauchy_nt_dd,
        GEN::pairedvector<int, double>& d_cauchy_nt_ds);

    /*!
     * @brief  integrate the structure residual and linearizations
     *
     * @tparam dim  dimension of the problem
     * @param[in] fac  pre-factor to correct sign dependent on integration of master or slave side
     *                 terms
     * @param[in,out] ele      mortar contact element or integration cell mortar element
     * @param[in] shape        shape function evaluated at current Gauss point
     * @param[in] shape_deriv  shape function derivative at current Gauss point
     * @param[in] d_xi_dd      directional derivative of Gauss point coordinates
     * @param[in] jac          Jacobian determinant of integration cell
     * @param[in] d_jac_dd     directional derivative of cell Jacobian
     * @param[in] wgt          Gauss point weight
     * @param[in] test_val     quantity to be integrated
     * @param[in] d_test_val_dd  directional derivative of quantity to be integrated
     * @param[in] d_test_val_ds  derivative of quantity to be integrated w.r.t. scalar s
     * @param[in] normal         normal
     * @param[in] d_normal_dd    directional derivative of normal
     */
    template <int dim>
    void IntegrateTest(double fac, MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_xi_dd, double jac,
        const GEN::pairedvector<int, double>& d_jac_dd, double wgt, double test_val,
        const GEN::pairedvector<int, double>& d_test_val_dd,
        const GEN::pairedvector<int, double>& d_test_val_ds, const LINALG::Matrix<dim, 1>& normal,
        const std::vector<GEN::pairedvector<int, double>>& d_normal_dd);

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
    void IntegrateScaTraTest(double fac, MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape_func, const LINALG::SerialDenseMatrix& shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_xi_dd, double jac,
        const GEN::pairedvector<int, double>& d_jac_dd, double wgt, double test_val,
        const GEN::pairedvector<int, double>& d_test_val_dd,
        const GEN::pairedvector<int, double>& d_test_val_ds);

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
    void IntegrateSSIInterfaceCondition(MORTAR::MortarElement& slave_ele,
        const LINALG::SerialDenseVector& slave_shape,
        const LINALG::SerialDenseMatrix& slave_shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_slave_xi_dd,
        MORTAR::MortarElement& master_ele, const LINALG::SerialDenseVector& master_shape,
        const LINALG::SerialDenseMatrix& master_shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_master_xi_dd, double jac,
        const GEN::pairedvector<int, double>& d_jac_dd, double wgt);

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
    void SetupGpConcentrations(MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape_func, const LINALG::SerialDenseMatrix& shape_deriv,
        const std::vector<GEN::pairedvector<int, double>>& d_xi_dd, double& gp_conc,
        GEN::pairedvector<int, double>& d_conc_dc, GEN::pairedvector<int, double>& d_conc_dd);

    //! scatra time integration factors
    const DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparamstimint_;
    //! scatra-scatra coupling specific parameters
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* scatraparamsboundary_;
  };
}  // namespace CONTACT
#endif  // CONTACT_NITSCHE_INTEGRATOR_SSI_H
