/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the poro contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_INTEGRATOR_PORO_HPP
#define FOUR_C_CONTACT_NITSCHE_INTEGRATOR_PORO_HPP

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
  class IntegratorNitschePoro : public IntegratorNitsche
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
    IntegratorNitschePoro(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm);

   protected:
    /*!
     \brief Perform integration at GP
            This is where the distinction between methods should be,
            i.e. mortar, augmented, gpts,...
     */
    void integrate_gp_2_d(MORTAR::Element& sele, MORTAR::Element& mele,
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
    void integrate_gp_3_d(MORTAR::Element& sele, MORTAR::Element& mele,
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
    \brief Evaluate cauchy stress component and its derivatives
    */
    template <int dim>
    void so_ele_cauchy(MORTAR::Element& moEle, double* boundary_gpcoord,
        std::vector<CORE::GEN::Pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::Pairedvector<int, double>>& normal_deriv,
        const CORE::LINALG::Matrix<dim, 1>& direction,
        std::vector<CORE::GEN::Pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, CORE::GEN::Pairedvector<int, double>& deriv_sigma_nt,
        CORE::GEN::Pairedvector<int, double>& deriv_sigma_nt_p);

   private:
    /*!
    \brief evaluate GPTS forces and linearization at this gp
    */
    template <int dim>
    void gpts_forces(MORTAR::Element& sele, MORTAR::Element& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dsxi,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dmxi, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double gap, const CORE::GEN::Pairedvector<int, double>& dgapgp, const double* gpn,
        std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi);

   protected:
    template <int dim>
    void integrate_test(const double fac, MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dxi, const double jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const CORE::GEN::Pairedvector<int, double>& test_deriv_d,
        const CORE::GEN::Pairedvector<int, double>& test_deriv_p,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void integrate_poro_no_out_flow(const double fac, MORTAR::Element& ele, double* xi,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const double jac, const CORE::GEN::Pairedvector<int, double>& jacintcellmap,
        const double wgt, const CORE::LINALG::Matrix<dim, 1>& normal,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& normal_deriv,
        MORTAR::Element& otherele, const CORE::LINALG::SerialDenseVector& othershape);

    bool get_poro_pressure(MORTAR::Element& ele, const CORE::LINALG::SerialDenseVector& shape,
        MORTAR::Element& otherele, const CORE::LINALG::SerialDenseVector& othershape,
        double& poropressure);

    void get_poro_quantitiesat_gp(MORTAR::Element& ele, double* xi, double& spresgp,  //(in)
        double& sJ, std::map<int, double>& sJLin, double& sporosity, double& sdphi_dp,
        double& sdphi_dJ);  // out

   private:
    bool no_penetration_;  // no outflow in contact zone ...
    double dv_dd_;         // 1/(theta*dt) for OST
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
