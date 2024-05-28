/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of ehl related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_EHL_INTEGRATOR_HPP
#define FOUR_C_CONTACT_EHL_INTEGRATOR_HPP

#include "4C_config.hpp"

#include "4C_contact_integrator.hpp"
#include "4C_utils_pairedvector.hpp"

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

   private:
    // integrate surface gradient
    void gp_weighted_surf_grad_and_deriv(MORTAR::Element& sele, const double* xi,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& dsxigp,
        const CORE::LINALG::SerialDenseVector& lmval,
        const CORE::LINALG::SerialDenseMatrix& lmderiv,
        const CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const CORE::LINALG::SerialDenseMatrix& sderiv2, const double& wgt, const double& jac,
        const CORE::GEN::Pairedvector<int, double>& jacintcellmap);

    // integrate relative and average tangential velocity
    void gp_weighted_av_rel_vel(MORTAR::Element& sele, MORTAR::Element& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseVector& lmval,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const CORE::LINALG::SerialDenseMatrix& mderiv,
        const CORE::LINALG::SerialDenseMatrix& lmderiv,
        const CORE::GEN::Pairedvector<int, CORE::LINALG::SerialDenseMatrix>& dualmap,
        const double& wgt, const double& jac, const CORE::GEN::Pairedvector<int, double>& derivjac,
        const double* normal, const std::vector<CORE::GEN::Pairedvector<int, double>>& dnmap_unit,
        const double& gap, const CORE::GEN::Pairedvector<int, double>& deriv_gap, const double* sxi,
        const double* mxi, const std::vector<CORE::GEN::Pairedvector<int, double>>& derivsxi,
        const std::vector<CORE::GEN::Pairedvector<int, double>>& derivmxi);
  };

}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
