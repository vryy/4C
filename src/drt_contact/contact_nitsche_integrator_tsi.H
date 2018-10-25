/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_integrator_tsi.H

\brief A class to perform integrations of nitsche related terms

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_INTEGRATOR_NITSCHE_TSI_H
#define CONTACT_INTEGRATOR_NITSCHE_TSI_H

#include "contact_nitsche_integrator.H"
#include "../headers/pairedvector.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>

// forward declarations
namespace LINALG
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
     DRT::UTILS::IntegrationPoints structs to get Gauss points
     and corresponding weights.

     */
    CoIntegratorNitscheTsi(Teuchos::ParameterList& params, DRT::Element::DiscretizationType eletype,
        const Epetra_Comm& comm)
        : CoIntegratorNitsche(params, eletype, comm),
          theta_thermo_(params.get<double>("NITSCHE_THETA_TSI")),
          nit_wgt_thermo_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING_TSI")),
          nit_thr_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheThermoMethod>(
              params, "NITSCHE_METHOD_TSI")),
          pp_thermo_(params.get<double>("PENALTYPARAM_THERMO")),
          temp_ref_(params.get<double>("TEMP_REF")),
          temp_damage_(params.get<double>("TEMP_DAMAGE")),
          gamma_slave_(params.get<double>("HEATTRANSSLAVE")),
          gamma_master_(params.get<double>("HEATTRANSMASTER"))
    {
    }

    /*!
     \brief Destructor

     */
    virtual ~CoIntegratorNitscheTsi() {}

   protected:
    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    virtual void IntegrateGP_2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
        GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
        GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<GEN::pairedvector<int, double>>& derivsxi,
        std::vector<GEN::pairedvector<int, double>>& derivmxi);

    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    virtual void IntegrateGP_3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseVector& mval, LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmderiv,
        GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
        GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<GEN::pairedvector<int, double>>& derivsxi,
        std::vector<GEN::pairedvector<int, double>>& derivmxi);

   private:
    /*!
    \brief evaluate GPTS forces and linearization at this gp
    */
    template <int dim>
    void GPTS_forces(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::SerialDenseVector& sval, const LINALG::SerialDenseMatrix& sderiv,
        const std::vector<GEN::pairedvector<int, double>>& dsxi,
        const LINALG::SerialDenseVector& mval, const LINALG::SerialDenseMatrix& mderiv,
        const std::vector<GEN::pairedvector<int, double>>& dmxi, const double jac,
        const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double gap,
        const GEN::pairedvector<int, double>& dgapgp, double* gpn,
        std::vector<GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi);

    template <int dim>
    void BuildAdjointTestTsi(MORTAR::MortarElement& moEle, const double fac,
        const Epetra_SerialDenseMatrix& d2sntDdDT,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void IntegrateTest(const double fac, MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& deriv,
        const std::vector<GEN::pairedvector<int, double>>& dxi, const double jac,
        const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const GEN::pairedvector<int, double>& test_deriv_d,
        const GEN::pairedvector<int, double>& test_deriv_T, const LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<GEN::pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void IntegrateAdjointTest(const double fac, const double jac,
        const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double test,
        const GEN::pairedvector<int, double>& deriv_test_d,
        const GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void IntegrateThermalTest(const double fac, MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& deriv,
        const std::vector<GEN::pairedvector<int, double>>& dxi, const double jac,
        const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const GEN::pairedvector<int, double>& test_deriv_d,
        const GEN::pairedvector<int, double>& test_deriv_T);

    template <int dim>
    void IntegrateThermalAdjointTest(const double fac, const double jac,
        const GEN::pairedvector<int, double>& jacintcellmap, const double wgt, const double test,
        const GEN::pairedvector<int, double>& deriv_test_d,
        const GEN::pairedvector<int, double>& deriv_test_T, MORTAR::MortarElement& moEle,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T);



    template <DRT::Element::DiscretizationType parentdistype, int dim>
    void inline SoEleGP(MORTAR::MortarElement& sele, const double wgt, const double* gpcoord,
        LINALG::Matrix<dim, 1>& pxsi, LINALG::Matrix<dim, dim>& derivtrafo);

    template <int dim>
    void SoEleCauchy(MORTAR::MortarElement& moEle, double* boundary_gpcoord,
        std::vector<GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const LINALG::Matrix<dim, 1>& normal,
        std::vector<GEN::pairedvector<int, double>>& normal_deriv,
        const LINALG::Matrix<dim, 1>& direction,
        std::vector<GEN::pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, GEN::pairedvector<int, double>& deriv_sigma_nt_d,
        GEN::pairedvector<int, double>& deriv_sigma_nt_T, LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void SoEleCauchyHeatflux(MORTAR::MortarElement& moEle, double* boundary_gpcoord,
        const std::vector<GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
        const double gp_wgt, const LINALG::Matrix<dim, 1>& normal,
        std::vector<GEN::pairedvector<int, double>>& normal_deriv, const double w, double& heatflux,
        GEN::pairedvector<int, double>& dq_dd, GEN::pairedvector<int, double>& dq_dT,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void BuildAdjointTestThermo(MORTAR::MortarElement& moEle, const double fac,
        const Epetra_SerialDenseMatrix& dq_dT_ele, const Epetra_SerialDenseMatrix& d2q_dT_dd,
        const Epetra_SerialDenseMatrix& d2q_dT_dn, const Epetra_SerialDenseMatrix& d2q_dT_dpxi,
        std::vector<GEN::pairedvector<int, double>>& normal_deriv,
        const std::vector<GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
        LINALG::Matrix<dim, dim>& derivtravo_slave, LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_d,
        GEN::pairedvector<int, LINALG::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void SetupGpTemp(MORTAR::MortarElement& moEle, const LINALG::SerialDenseVector& val,
        const LINALG::SerialDenseMatrix& deriv,
        const std::vector<GEN::pairedvector<int, double>> dxi, double& temp,
        GEN::pairedvector<int, double>& d_temp_dT, GEN::pairedvector<int, double>& d_temp_dd);

   private:
    double theta_thermo_;
    INPAR::CONTACT::NitscheWeighting nit_wgt_thermo_;
    INPAR::CONTACT::NitscheThermoMethod nit_thr_;
    double pp_thermo_;
    double temp_ref_;
    double temp_damage_;
    double gamma_slave_;
    double gamma_master_;
  };
}  // namespace CONTACT
#endif  // CONTACT_INTEGRATOR_NITSCHE_TSI_H
