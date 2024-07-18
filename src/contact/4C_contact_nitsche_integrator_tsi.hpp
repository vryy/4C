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

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
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
     Core::FE::IntegrationPoints structs to get Gauss points
     and corresponding weights.

     */
    IntegratorNitscheTsi(
        Teuchos::ParameterList& params, Core::FE::CellType eletype, const Epetra_Comm& comm)
        : IntegratorNitsche(params, eletype, comm),
          theta_thermo_(params.get<double>("NITSCHE_THETA_TSI")),
          nit_thr_(Core::UTILS::IntegralValue<Inpar::CONTACT::NitscheThermoMethod>(
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
    void integrate_gp_2d(Mortar::Element& sele, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
        Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
        Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
        std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi) override;

    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    void integrate_gp_3d(Mortar::Element& sele, Mortar::Element& mele,
        Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
        Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
        Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
        double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
        Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
        std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi) override;

   private:
    /*!
    \brief evaluate GPTS forces and linearization at this gp
    */
    template <int dim>
    void gpts_forces(Mortar::Element& sele, Mortar::Element& mele,
        const Core::LinAlg::SerialDenseVector& sval, const Core::LinAlg::SerialDenseMatrix& sderiv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& dsxi,
        const Core::LinAlg::SerialDenseVector& mval, const Core::LinAlg::SerialDenseMatrix& mderiv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& dmxi, const double jac,
        const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double gap, const Core::Gen::Pairedvector<int, double>& dgapgp, const double* gpn,
        std::vector<Core::Gen::Pairedvector<int, double>>& deriv_contact_normal, double* sxi,
        double* mxi);

    template <int dim>
    void build_adjoint_test_tsi(Mortar::Element& moEle, const double fac,
        const Core::LinAlg::SerialDenseMatrix& d2sntDdDT,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void integrate_test(const double fac, Mortar::Element& ele,
        const Core::LinAlg::SerialDenseVector& shape, const Core::LinAlg::SerialDenseMatrix& deriv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& dxi, const double jac,
        const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const Core::Gen::Pairedvector<int, double>& test_deriv_d,
        const Core::Gen::Pairedvector<int, double>& test_deriv_T,
        const Core::LinAlg::Matrix<dim, 1>& test_dir,
        const std::vector<Core::Gen::Pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void integrate_adjoint_test(const double fac, const double jac,
        const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test, const Core::Gen::Pairedvector<int, double>& deriv_test_d,
        const Core::Gen::Pairedvector<int, double>& deriv_test_T, Mortar::Element& moEle,
        Core::LinAlg::SerialDenseVector& adjoint_test,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void integrate_thermal_test(const double fac, Mortar::Element& ele,
        const Core::LinAlg::SerialDenseVector& shape, const Core::LinAlg::SerialDenseMatrix& deriv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& dxi, const double jac,
        const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const Core::Gen::Pairedvector<int, double>& test_deriv_d,
        const Core::Gen::Pairedvector<int, double>& test_deriv_T);

    template <int dim>
    void integrate_thermal_adjoint_test(const double fac, const double jac,
        const Core::Gen::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test, const Core::Gen::Pairedvector<int, double>& deriv_test_d,
        const Core::Gen::Pairedvector<int, double>& deriv_test_T, Mortar::Element& moEle,
        Core::LinAlg::SerialDenseVector& adjoint_test,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T);



    template <Core::FE::CellType parentdistype, int dim>
    void inline so_ele_gp(Mortar::Element& sele, const double wgt, const double* gpcoord,
        Core::LinAlg::Matrix<dim, 1>& pxsi, Core::LinAlg::Matrix<dim, dim>& derivtrafo);

    template <int dim>
    void so_ele_cauchy(Mortar::Element& moEle, double* boundary_gpcoord,
        std::vector<Core::Gen::Pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const Core::LinAlg::Matrix<dim, 1>& normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& normal_deriv,
        const Core::LinAlg::Matrix<dim, 1>& direction,
        std::vector<Core::Gen::Pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, Core::Gen::Pairedvector<int, double>& deriv_sigma_nt_d,
        Core::Gen::Pairedvector<int, double>& deriv_sigma_nt_T,
        Core::LinAlg::SerialDenseVector& adjoint_test,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void so_ele_cauchy_heatflux(Mortar::Element& moEle, double* boundary_gpcoord,
        const std::vector<Core::Gen::Pairedvector<int, double>>& boundary_gpcoord_lin,
        const double gp_wgt, const Core::LinAlg::Matrix<dim, 1>& normal,
        std::vector<Core::Gen::Pairedvector<int, double>>& normal_deriv, const double w,
        double& heatflux, Core::Gen::Pairedvector<int, double>& dq_dd,
        Core::Gen::Pairedvector<int, double>& dq_dT, Core::LinAlg::SerialDenseVector& adjoint_test,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void build_adjoint_test_thermo(Mortar::Element& moEle, const double fac,
        const Core::LinAlg::SerialDenseMatrix& dq_dT_ele,
        const Core::LinAlg::SerialDenseMatrix& d2q_dT_dd,
        const Core::LinAlg::SerialDenseMatrix& d2q_dT_dn,
        const Core::LinAlg::SerialDenseMatrix& d2q_dT_dpxi,
        std::vector<Core::Gen::Pairedvector<int, double>>& normal_deriv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& boundary_gpcoord_lin,
        Core::LinAlg::Matrix<dim, dim>& derivtravo_slave,
        Core::LinAlg::SerialDenseVector& adjoint_test,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_d,
        Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseVector>& deriv_adjoint_test_T);

    template <int dim>
    void setup_gp_temp(Mortar::Element& moEle, const Core::LinAlg::SerialDenseVector& val,
        const Core::LinAlg::SerialDenseMatrix& deriv,
        const std::vector<Core::Gen::Pairedvector<int, double>>& dxi, double& temp,
        Core::Gen::Pairedvector<int, double>& d_temp_dT,
        Core::Gen::Pairedvector<int, double>& d_temp_dd);

   private:
    double theta_thermo_;
    Inpar::CONTACT::NitscheThermoMethod nit_thr_;
    double pp_thermo_;
    double temp_ref_;
    double temp_damage_;
    double gamma_slave_;
    double gamma_master_;
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
