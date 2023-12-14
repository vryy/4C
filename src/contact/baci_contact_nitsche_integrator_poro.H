/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the poro contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_INTEGRATOR_PORO_H
#define BACI_CONTACT_NITSCHE_INTEGRATOR_PORO_H

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
  class CoIntegratorNitschePoro : public CoIntegratorNitsche
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
    CoIntegratorNitschePoro(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm);

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

    /*!
    \brief Evaluate cauchy stress component and its derivatives
    */
    template <int dim>
    void SoEleCauchy(MORTAR::MortarElement& moEle, double* boundary_gpcoord,
        std::vector<CORE::GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
        const CORE::LINALG::Matrix<dim, 1>& direction,
        std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, CORE::GEN::pairedvector<int, double>& deriv_sigma_nt,
        CORE::GEN::pairedvector<int, double>& deriv_sigma_nt_p);

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
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double* sxi, double* mxi);

   protected:
    template <int dim>
    void IntegrateTest(const double fac, MORTAR::MortarElement& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, const double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const CORE::GEN::pairedvector<int, double>& test_deriv_d,
        const CORE::GEN::pairedvector<int, double>& test_deriv_p,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void IntegratePoroNoOutFlow(const double fac, MORTAR::MortarElement& ele, double* xi,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const double jac, const CORE::GEN::pairedvector<int, double>& jacintcellmap,
        const double wgt, const CORE::LINALG::Matrix<dim, 1>& normal,
        const std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
        MORTAR::MortarElement& otherele, const CORE::LINALG::SerialDenseVector& othershape);

    bool GetPoroPressure(MORTAR::MortarElement& ele, const CORE::LINALG::SerialDenseVector& shape,
        MORTAR::MortarElement& otherele, const CORE::LINALG::SerialDenseVector& othershape,
        double& poropressure);

    void GetPoroQuantitiesatGP(MORTAR::MortarElement& ele, double* xi, double& spresgp,  //(in)
        double& sJ, std::map<int, double>& sJLin, double& sporosity, double& sdphi_dp,
        double& sdphi_dJ);  // out

   private:
    bool no_penetration_;  // no outflow in contact zone ...
    double dv_dd_;         // 1/(theta*dt) for OST
  };
}  // namespace CONTACT
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_INTEGRATOR_PORO_H
