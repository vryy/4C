/*---------------------------------------------------------------------*/
/*!
\brief A class to perform integrations of ehl related terms

\level 3

\maintainer Mostafa Faraji

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_INTEGRATOR_EHL_H
#define CONTACT_INTEGRATOR_EHL_H

#include "contact_integrator.H"
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
  class CoIntegratorEhl : public CONTACT::CoIntegrator
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
    CoIntegratorEhl(Teuchos::ParameterList& params, DRT::Element::DiscretizationType eletype,
        const Epetra_Comm& comm)
        : CoIntegrator(params, eletype, comm)
    {
    }

    /*!
     \brief Destructor

     */
    virtual ~CoIntegratorEhl() {}

   protected:
    /*!
     \brief Perform integration at GP
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
    // integrate surface gradient
    void GP_WeightedSurfGradAndDeriv(MORTAR::MortarElement& sele, const double* xi,
        const std::vector<GEN::pairedvector<int, double>>& dsxigp,
        const LINALG::SerialDenseVector& lmval, const LINALG::SerialDenseMatrix& lmderiv,
        const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap,
        const LINALG::SerialDenseVector& sval, const LINALG::SerialDenseMatrix& sderiv,
        const LINALG::SerialDenseMatrix& sderiv2, const double& wgt, const double& jac,
        const GEN::pairedvector<int, double>& jacintcellmap);

    // integrate relative and average tangential velocity
    void GP_WeightedAvRelVel(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::SerialDenseVector& sval, const LINALG::SerialDenseVector& lmval,
        const LINALG::SerialDenseVector& mval, const LINALG::SerialDenseMatrix& sderiv,
        const LINALG::SerialDenseMatrix& mderiv, const LINALG::SerialDenseMatrix& lmderiv,
        const GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, const double& wgt,
        const double& jac, const GEN::pairedvector<int, double>& derivjac, const double* normal,
        const std::vector<GEN::pairedvector<int, double>>& dnmap_unit, const double& gap,
        const GEN::pairedvector<int, double>& deriv_gap, const double* sxi, const double* mxi,
        const std::vector<GEN::pairedvector<int, double>>& derivsxi,
        const std::vector<GEN::pairedvector<int, double>>& derivmxi);
  };

}  // namespace CONTACT
#endif  // CONTACT_INTEGRATOR_EHL_H
