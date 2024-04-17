/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of ehl related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_EHL_INTEGRATOR_HPP
#define FOUR_C_CONTACT_EHL_INTEGRATOR_HPP

#include "baci_config.hpp"

#include "baci_contact_integrator.hpp"
#include "baci_utils_pairedvector.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace CONTACT
{
  class IntegratorEhl : public CONTACT::Integrator
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
    IntegratorEhl(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
        : Integrator(params, eletype, comm)
    {
    }


   protected:
    /*!
     \brief Perform integration at GP
     */
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

    /*!
     \brief Perform integration at GP
     */
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

   private:
    // integrate surface gradient
    void GP_WeightedSurfGradAndDeriv(MORTAR::Element& sele, const double* xi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxigp,
        const CORE::LINALG::SerialDenseVector& lmval,
        const CORE::LINALG::SerialDenseMatrix& lmderiv,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const CORE::LINALG::SerialDenseMatrix& sderiv2, const double& wgt, const double& jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap);

    // integrate relative and average tangential velocity
    void GP_WeightedAvRelVel(MORTAR::Element& sele, MORTAR::Element& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseVector& lmval,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const CORE::LINALG::SerialDenseMatrix& mderiv,
        const CORE::LINALG::SerialDenseMatrix& lmderiv,
        const CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        const double& wgt, const double& jac, const CORE::GEN::pairedvector<int, double>& derivjac,
        const double* normal, const std::vector<CORE::GEN::pairedvector<int, double>>& dnmap_unit,
        const double& gap, const CORE::GEN::pairedvector<int, double>& deriv_gap, const double* sxi,
        const double* mxi, const std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi);
  };

}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
