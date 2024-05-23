/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_INTEGRATOR_TSI_HPP
#define FOUR_C_CONTACT_NITSCHE_INTEGRATOR_TSI_HPP

#include "4C_config.hpp"

#include "4C_contact_nitsche_integrator.hpp"
#include "4C_utils_pairedvector.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
}

namespace CONTACT
{
  class IntegratorNitscheTsi : public CONTACT::IntegratorNitsche
  {
   public:
    /*!
     \brief Constructor  with shape function specification

     Constructs an instance of this class using a specific type of shape functions.<br>
     Note that this is \b not a collective call as overlaps are
     integrated in parallel by individual processes.<br>
     Note also that this constructor relies heavily on the
     CORE::FE::IntegrationPoints structs to get Gauss points
     and corresponding weights.

     */
    IntegratorNitscheTsi(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
        : IntegratorNitsche(params, eletype, comm),
          theta_thermo_(params.get<double>("NITSCHE_THETA_TSI")),
          nit_thr_(CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheThermoMethod>(
              params, "NITSCHE_METHOD_TSI")),
          pp_thermo_(params.get<double>("PENALTYPARAM_THERMO")),
          temp_ref_(params.get<double>("TEMP_REF")),
          temp_damage_(params.get<double>("TEMP_DAMAGE")),
          gamma_slave_(params.get<double>("HEATTRANSSLAVE")),
          gamma_master_(params.get<double>("HEATTRANSMASTER"))
    {
    }

   protected:
    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    void IntegrateGP_2D(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::Pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::Pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::Pairedvector<int, double>>& derivmxi) override;

    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    void IntegrateGP_3D(MORTAR::Element& sele, MORTAR::Element& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::Pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::Pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::Pairedvector<int, double>>& derivmxi) override;

   private:
    /*!
    \brief evaluate GPTS forces and linearization at this gp
    */
    template <int dim>
    void GPTSForces(MORTAR::Element& sele, MORTAR::Element& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dsxi,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dmxi, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double gap, const CORE::GEN::Pairedvector<int, double>& dgapgp, const double* gpn,
        std::vector<CORE::GEN::Pairedvector<int, double>>& deriv_contact_normal, double* sxi,
        double* mxi);

    template <int dim>
    void BuildAdjointTestTsi(MORTAR::Element& moEle, const double fac,
        const CORE::LINALG::SerialDenseMatrix& d2sntDdDT,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void IntegrateTest(const double fac, MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dxi, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const CORE::GEN::Pairedvector<int, double>& test_deriv_d,
        const CORE::GEN::Pairedvector<int, double>& test_deriv_T,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void integrate_adjoint_test(const double fac, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test, const CORE::GEN::Pairedvector<int, double>& deriv_test_d,
        const CORE::GEN::Pairedvector<int, double>& deriv_test_T, MORTAR::Element& moEle,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void integrate_thermal_test(const double fac, MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dxi, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const CORE::GEN::Pairedvector<int, double>& test_deriv_d,
        const CORE::GEN::Pairedvector<int, double>& test_deriv_T);

    template <int dim>
    void integrate_thermal_adjoint_test(const double fac, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test, const CORE::GEN::Pairedvector<int, double>& deriv_test_d,
        const CORE::GEN::Pairedvector<int, double>& deriv_test_T, MORTAR::Element& moEle,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);



    template <CORE::FE::CellType parentdistype, int dim>
    void inline SoEleGP(MORTAR::Element& sele, const double wgt, const double* gpcoord,
        CORE::LINALG::Matrix<dim, 1>& pxsi, CORE::LINALG::Matrix<dim, dim>& derivtrafo);

    template <int dim>
    void SoEleCauchy(MORTAR::Element& moEle, double* boundary_gpcoord,
        std::vector<CORE::GEN::Pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normal_deriv,
        const CORE::LINALG::Matrix<dim, 1>& direction,
        std::vector<CORE::GEN::Pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, CORE::GEN::Pairedvector<int, double>& deriv_sigma_nt_d,
        CORE::GEN::Pairedvector<int, double>& deriv_sigma_nt_T,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void SoEleCauchyHeatflux(MORTAR::Element& moEle, double* boundary_gpcoord,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& boundary_gpcoord_lin,
        const double gp_wgt, const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normal_deriv, const double w,
        double& heatflux, CORE::GEN::Pairedvector<int, double>& dq_dd,
        CORE::GEN::Pairedvector<int, double>& dq_dT, CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void build_adjoint_test_thermo(MORTAR::Element& moEle, const double fac,
        const CORE::LINALG::SerialDenseMatrix& dq_dT_ele,
        const CORE::LINALG::SerialDenseMatrix& d2q_dT_dd,
        const CORE::LINALG::SerialDenseMatrix& d2q_dT_dn,
        const CORE::LINALG::SerialDenseMatrix& d2q_dT_dpxi,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normal_deriv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& boundary_gpcoord_lin,
        CORE::LINALG::Matrix<dim, dim>& derivtravo_slave,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void SetupGpTemp(MORTAR::Element& moEle, const CORE::LINALG::SerialDenseVector& val,
        const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dxi, double& temp,
        CORE::GEN::Pairedvector<int, double>& d_temp_dT,
        CORE::GEN::Pairedvector<int, double>& d_temp_dd);

   private:
    double theta_thermo_;
    INPAR::CONTACT::NitscheThermoMethod nit_thr_;
    double pp_thermo_;
    double temp_ref_;
    double temp_damage_;
    double gamma_slave_;
    double gamma_master_;
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
