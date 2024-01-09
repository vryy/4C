/*----------------------------------------------------------------------------*/
/*! \file
\brief A class to perform integration of nitsche related terms for the ssi contact case

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_INTEGRATOR_SSI_H
#define BACI_CONTACT_NITSCHE_INTEGRATOR_SSI_H

#include "baci_config.H"

#include "baci_contact_nitsche_integrator.H"

BACI_NAMESPACE_OPEN

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
  /*!
   * @brief This class performs Gauss integration and the assembly to element matrices and vectors
   * that are relevant to the Nitsche contact formulation for scatra-structure interaction problems.
   *
   * @note Relevant methods are already templated w.r.t. the problem dimension. Currently only
   * 'dim=3' is used and tested but it should be quite easy to extend this if necessary.
   */
  class IntegratorNitscheSsi : public IntegratorNitsche
  {
   public:
    /*!
     * \brief Constructor with shape function specification
     *
     * Constructs an instance of this class using a specific type of shape functions.<br> Note that
     * this is \b not a collective call as overlaps are integrated in parallel by individual
     * processes.<br> Note also that this constructor relies heavily on the
     * CORE::FE::IntegrationPoints structs to get Gauss points and corresponding weights.
     *
     * @param[in] params   interface contact parameter list
     * @param[in] eletype  shape of integration cell for segment based integration or slave side
     *                     mortar contact element for element based integration
     * @param[in] comm     contact interface communicator
     */
    IntegratorNitscheSsi(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm);

   protected:
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
     * @param[out] d_sigma_nt_ds  derivative of cauchy_nt w.r.t. scalar s \f[ \frac{ \mathrm{d}
     *                             \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                             \boldsymbol{t}}{\mathrm{d} s} \f]
     */
    template <int dim>
    void SoEleCauchyStruct(MORTAR::Element& mortar_ele, double* gp_coord,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_gp_coord_dd, double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& gp_normal,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_gp_normal_dd,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_test_dir_dd, double nitsche_wgt,
        double& cauchy_nt_wgt, CORE::GEN::pairedvector<int, double>& d_cauchy_nt_dd,
        CORE::LINALG::SerialDenseMatrix* d_sigma_nt_ds);

   private:
    void IntegrateGP_2D(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi) override;

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
    void GPTSForces(MORTAR::Element& slave_ele, MORTAR::Element& master_ele,
        const CORE::LINALG::SerialDenseVector& slave_shape,
        const CORE::LINALG::SerialDenseMatrix& slave_shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_slave_xi_dd,
        const CORE::LINALG::SerialDenseVector& master_shape,
        const CORE::LINALG::SerialDenseMatrix& master_shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_master_xi_dd, double jac,
        const CORE::GEN::pairedvector<int, double>& d_jac_dd, double gp_wgt, double gap,
        const CORE::GEN::pairedvector<int, double>& d_gap_dd, const double* gp_normal,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_gp_normal_dd, double* slave_xi,
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
    void SoEleCauchy(MORTAR::Element& mortar_ele, double* gp_coord,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_gp_coord_dd, double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& gp_normal,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_gp_normal_dd,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_test_dir_dd, double nitsche_wgt,
        double& cauchy_nt_wgt, CORE::GEN::pairedvector<int, double>& d_cauchy_nt_dd,
        CORE::GEN::pairedvector<int, double>& d_cauchy_nt_ds);

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
    void IntegrateTest(double fac, MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape,
        const CORE::LINALG::SerialDenseMatrix& shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_xi_dd, double jac,
        const CORE::GEN::pairedvector<int, double>& d_jac_dd, double wgt, double test_val,
        const CORE::GEN::pairedvector<int, double>& d_test_val_dd,
        const CORE::GEN::pairedvector<int, double>& d_test_val_ds,
        const CORE::LINALG::Matrix<dim, 1>& normal,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_normal_dd);

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
        const CORE::GEN::pairedvector<int, double>& d_test_val_ds);

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
     * @param[in] cauchy_nn_average_pen_gap         normal contact stress
     * @param[in] d_cauchy_nn_weighted_average_dd   directional derivatives of normal contact stress
     * w.r.t displacement
     * @param[in] d_cauchy_nn_weighted_average_dc   directional derivatives of normal contact stress
     * w.r.t scalar/concentration
     * @param[in] jac                 Jacobian determinant of integration cell
     * @param[in] d_jac_dd            directional derivative of cell Jacobian
     * @param[in] wgt                 Gauss point weight
     */
    template <int dim>
    void IntegrateSSIInterfaceCondition(MORTAR::Element& slave_ele,
        const CORE::LINALG::SerialDenseVector& slave_shape,
        const CORE::LINALG::SerialDenseMatrix& slave_shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_slave_xi_dd,
        MORTAR::Element& master_ele, const CORE::LINALG::SerialDenseVector& master_shape,
        const CORE::LINALG::SerialDenseMatrix& master_shape_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& d_master_xi_dd,
        const double cauchy_nn_average_pen_gap,
        const CORE::GEN::pairedvector<int, double>& d_cauchy_nn_weighted_average_dd,
        const CORE::GEN::pairedvector<int, double>& d_cauchy_nn_weighted_average_dc, double jac,
        const CORE::GEN::pairedvector<int, double>& d_jac_dd, double wgt);

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

    //! scatra time integration factors
    const DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparamstimint_;
    //! scatra-scatra coupling specific parameters
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* scatraparamsboundary_;
  };
}  // namespace CONTACT
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_INTEGRATOR_SSI_H
