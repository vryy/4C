/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_INTEGRATOR_TSI_H
#define BACI_CONTACT_NITSCHE_INTEGRATOR_TSI_H

#include "baci_contact_nitsche_integrator.H"
#include "baci_utils_pairedvector.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
}

namespace CONTACT
{
  class CoIntegratorNitscheTsi : public CONTACT::CoIntegratorNitsche
  {
   public:
    /*!
     \brief Constructor  with shape function specification

     Constructs an instance of this class using a specific type of shape functions.<br>
     Note that this is \b not a collective call as overlaps are
     integrated in parallel by individual processes.<br>
     Note also that this constructor relies heavily on the
     CORE::DRT::UTILS::IntegrationPoints structs to get Gauss points
     and corresponding weights.

     */
    CoIntegratorNitscheTsi(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
        : CoIntegratorNitsche(params, eletype, comm),
          theta_thermo_(params.get<double>("NITSCHE_THETA_TSI")),
          nit_thr_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheThermoMethod>(
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
    void IntegrateGP_2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
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
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    void IntegrateGP_3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi) override;

   private:
    /*!
    \brief evaluate GPTS forces and linearization at this gp
    */
    template <int dim>
    void GPTSForces(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxi,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double gap, const CORE::GEN::pairedvector<int, double>& dgapgp, const double* gpn,
        std::vector<CORE::GEN::pairedvector<int, double>>& deriv_contact_normal, double* sxi,
        double* mxi);

    template <int dim>
    void BuildAdjointTestTsi(MORTAR::MortarElement& moEle, const double fac,
        const CORE::LINALG::SerialDenseMatrix& d2sntDdDT,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void IntegrateTest(const double fac, MORTAR::MortarElement& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const CORE::GEN::pairedvector<int, double>& test_deriv_d,
        const CORE::GEN::pairedvector<int, double>& test_deriv_T,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void IntegrateAdjointTest(const double fac, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test, const CORE::GEN::pairedvector<int, double>& deriv_test_d,
        const CORE::GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void IntegrateThermalTest(const double fac, MORTAR::MortarElement& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const CORE::GEN::pairedvector<int, double>& test_deriv_d,
        const CORE::GEN::pairedvector<int, double>& test_deriv_T);

    template <int dim>
    void IntegrateThermalAdjointTest(const double fac, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test, const CORE::GEN::pairedvector<int, double>& deriv_test_d,
        const CORE::GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);



    template <CORE::FE::CellType parentdistype, int dim>
    void inline SoEleGP(MORTAR::MortarElement& sele, const double wgt, const double* gpcoord,
        CORE::LINALG::Matrix<dim, 1>& pxsi, CORE::LINALG::Matrix<dim, dim>& derivtrafo);

    template <int dim>
    void SoEleCauchy(MORTAR::MortarElement& moEle, double* boundary_gpcoord,
        std::vector<CORE::GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
        const CORE::LINALG::Matrix<dim, 1>& direction,
        std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_d,
        CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_T,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void SoEleCauchyHeatflux(MORTAR::MortarElement& moEle, double* boundary_gpcoord,
        const std::vector<CORE::GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
        const double gp_wgt, const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv, const double w,
        double& heatflux, CORE::GEN::pairedvector<int, double>& dq_dd,
        CORE::GEN::pairedvector<int, double>& dq_dT, CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void BuildAdjointTestThermo(MORTAR::MortarElement& moEle, const double fac,
        const CORE::LINALG::SerialDenseMatrix& dq_dT_ele,
        const CORE::LINALG::SerialDenseMatrix& d2q_dT_dd,
        const CORE::LINALG::SerialDenseMatrix& d2q_dT_dn,
        const CORE::LINALG::SerialDenseMatrix& d2q_dT_dpxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
        CORE::LINALG::Matrix<dim, dim>& derivtravo_slave,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void SetupGpTemp(MORTAR::MortarElement& moEle, const CORE::LINALG::SerialDenseVector& val,
        const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, double& temp,
        CORE::GEN::pairedvector<int, double>& d_temp_dT,
        CORE::GEN::pairedvector<int, double>& d_temp_dd);

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
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_INTEGRATOR_TSI_H
