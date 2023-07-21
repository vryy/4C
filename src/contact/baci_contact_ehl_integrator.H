/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of ehl related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_EHL_INTEGRATOR_H
#define BACI_CONTACT_EHL_INTEGRATOR_H

#include "baci_contact_integrator.H"
#include "baci_utils_pairedvector.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>

// forward declarations
namespace CORE::LINALG
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
     CORE::DRT::UTILS::IntegrationPoints structs to get Gauss points
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
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
        CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi);

    /*!
     \brief Perform integration at GP
     */
    virtual void IntegrateGP_3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        CORE::LINALG::SerialDenseVector& sval, CORE::LINALG::SerialDenseVector& lmval,
        CORE::LINALG::SerialDenseVector& mval, CORE::LINALG::SerialDenseMatrix& sderiv,
        CORE::LINALG::SerialDenseMatrix& mderiv, CORE::LINALG::SerialDenseMatrix& lmderiv,
        CORE::GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, double& wgt, double& jac,
        CORE::GEN::pairedvector<int, double>& derivjac, double* normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit, double& gap,
        CORE::GEN::pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi);

   private:
    // integrate surface gradient
    void GP_WeightedSurfGradAndDeriv(MORTAR::MortarElement& sele, const double* xi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::LINALG::SerialDenseVector& lmval,
        const CORE::LINALG::SerialDenseMatrix& lmderiv,
        const CORE::GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const CORE::LINALG::SerialDenseMatrix& sderiv2, const double& wgt, const double& jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap);

    // integrate relative and average tangential velocity
    void GP_WeightedAvRelVel(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseVector& lmval,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const CORE::LINALG::SerialDenseMatrix& mderiv,
        const CORE::LINALG::SerialDenseMatrix& lmderiv,
        const CORE::GEN::pairedvector<int, Epetra_SerialDenseMatrix>& dualmap, const double& wgt,
        const double& jac, const CORE::GEN::pairedvector<int, double>& derivjac,
        const double* normal, const std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit,
        const double& gap, const CORE::GEN::pairedvector<int, double>& deriv_gap, const double* sxi,
        const double* mxi, const std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi);
  };

}  // namespace CONTACT
#endif  // CONTACT_EHL_INTEGRATOR_H
