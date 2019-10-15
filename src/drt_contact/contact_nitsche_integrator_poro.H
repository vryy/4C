/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms for the poro contact case

\level 3

\maintainer Christoph Schmidt

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_INTEGRATOR_NITSCHE_PORO_H
#define CONTACT_INTEGRATOR_NITSCHE_PORO_H

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
  class CoIntegratorNitschePoro : public CoIntegratorNitsche
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
    CoIntegratorNitschePoro(Teuchos::ParameterList& params,
        DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm);

    /*!
     \brief Destructor

     */
    virtual ~CoIntegratorNitschePoro() {}

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

    /*!
    \brief Evaluate cauchy stress component and its derivatives
    */
    template <int dim>
    void SoEleCauchy(MORTAR::MortarElement& moEle, double* boundary_gpcoord,
        std::vector<GEN::pairedvector<int, double>> boundary_gpcoord_lin, const double gp_wgt,
        const LINALG::Matrix<dim, 1>& normal,
        std::vector<GEN::pairedvector<int, double>>& normal_deriv,
        const LINALG::Matrix<dim, 1>& direction,
        std::vector<GEN::pairedvector<int, double>>& direction_deriv, const double w,
        double& cauchy_nt, GEN::pairedvector<int, double>& deriv_sigma_nt,
        GEN::pairedvector<int, double>& deriv_sigma_nt_p);

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

   protected:
    template <int dim>
    void IntegrateTest(const double fac, MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& deriv,
        const std::vector<GEN::pairedvector<int, double>>& dxi, const double jac,
        const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const double test_val, const GEN::pairedvector<int, double>& test_deriv_d,
        const GEN::pairedvector<int, double>& test_deriv_p, const LINALG::Matrix<dim, 1>& normal,
        const std::vector<GEN::pairedvector<int, double>>& normal_deriv);

    template <int dim>
    void IntegratePoroNoOutFlow(const double fac, MORTAR::MortarElement& ele, double* xi,
        const LINALG::SerialDenseVector& shape, const LINALG::SerialDenseMatrix& deriv,
        const double jac, const GEN::pairedvector<int, double>& jacintcellmap, const double wgt,
        const LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<GEN::pairedvector<int, double>>& test_dir_deriv,
        MORTAR::MortarElement& otherele, const LINALG::SerialDenseVector& othershape);

    bool GetPoroPressure(MORTAR::MortarElement& ele, const LINALG::SerialDenseVector& shape,
        MORTAR::MortarElement& otherele, const LINALG::SerialDenseVector& othershape,
        double& poropressure);

    void GetPoroQuantitiesatGP(MORTAR::MortarElement& ele, double* xi, double& savpresgp,  //(in)
        double& sJ, std::map<int, double>& sJLin, double& sporosity, double& sdphi_dp,
        double& sdphi_dJ);  // out

   private:
    bool no_penetration_;  // no outflow in contact zone ...
    double dv_dd_;         // 1/(theta*dt) for OST
  };
}  // namespace CONTACT
#endif  // CONTACT_INTEGRATOR_NITSCHE_PORO_H
