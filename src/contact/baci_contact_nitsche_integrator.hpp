/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of nitsche related terms

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_INTEGRATOR_HPP
#define FOUR_C_CONTACT_NITSCHE_INTEGRATOR_HPP

#include "baci_config.hpp"

#include "baci_contact_integrator.hpp"
#include "baci_utils_pairedvector.hpp"

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
  class IntegratorNitsche : public CONTACT::Integrator
  {
   public:
    /*!
     * \brief Constructor with shape function specification
     *
     * Constructs an instance of this class using a specific type of shape functions.<br> Note that
     * this is \b not a collective call as overlaps are integrated in parallel by individual
     * processes.<br> Note also that this constructor relies heavily on the
     * CORE::FE::IntegrationPoints structs to get Gauss points and corresponding weights.
     *
     * @param [in] params  interface contact parameter list
     * @param [in] eletype shape of integration cell for segment based integration or slave side
     *                     mortar contact element for element based integration
     * @param [in] comm    contact interface communicator
     */
    IntegratorNitsche(
        Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm)
        : Integrator(params, eletype, comm),
          theta_(params.get<double>("NITSCHE_THETA")),
          theta_2_(params.get<double>("NITSCHE_THETA_2")),
          nit_wgt_(CORE::UTILS::IntegralValue<INPAR::CONTACT::NitscheWeighting>(
              params, "NITSCHE_WEIGHTING")),
          ppn_(imortar_.get<double>("PENALTYPARAM")),
          ppt_(imortar_.get<double>("PENALTYPARAMTAN")),
          frcoeff_(imortar_.get<double>("FRCOEFF", -1.)),
          frbound_(imortar_.get<double>("FRBOUND", -1.)),
          frtype_(CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(imortar_, "FRICTION")),
          dt_(imortar_.get<double>("TIMESTEP")),
          two_half_pass_(imortar_.get<bool>("Two_half_pass"))
    {
    }

   protected:
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

    /*!
     * @brief Evaluate cauchy stress component and its derivatives
     *
     * @tparam dim   dimension of the problem
     * @param[in] moEle  mortar element
     * @param[in] boundary_gpcoord      Gauss point coordinates
     * @param[in] boundary_gpcoord_lin  directional derivative of Gauss point coordinates
     * @param[in] gp_wgt        Gauss point weight
     * @param[in] normal        normal
     * @param[in] normal_deriv  directional derivative of normal
     * @param[in] direction        direction of evaluation (e.g. normal or tangential direction)
     * @param[in] direction_deriv  directional derivative of direction of evaluation
     * @param[in] w                Nitsche weight
     * @param[out] cauchy_nt        Cauchy stress tensor contracted with normal vector and direction
     *                              vector \f[ \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                              \boldsymbol{t} \f]
     * @param[out] deriv_sigma_nt   directional derivative of cauchy_nt \f[ \frac{ \mathrm{d}
     *                              \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                              \boldsymbol{t}}{\mathrm{d} \boldsymbol{d}} \f]
     * @param[out] adjoint_test     directional derivative of cauchy_nt \f[ \frac{ \mathrm{d}
     *                              \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                              \boldsymbol{t}}{\mathrm{d} \boldsymbol{d}} \f]
     * @param[out] deriv_adjoint_test   second directional derivative of cauchy_nt \f[ \frac{
     *                                  \mathrm{d}^2 \boldsymbol{\sigma} \cdot \boldsymbol{n} \cdot
     *                                  \boldsymbol{t}}{\mathrm{d} \boldsymbol{d}^2} \f]
     */
    template <int dim>
    void SoEleCauchy(MORTAR::Element& moEle, double* boundary_gpcoord,
        std::vector<CORE::GEN::pairedvector<int, double>> boundary_gpcoord_lin, double gp_wgt,
        const CORE::LINALG::Matrix<dim, 1>& normal,
        std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
        const CORE::LINALG::Matrix<dim, 1>& direction,
        std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv, double w,
        double& cauchy_nt, CORE::GEN::pairedvector<int, double>& deriv_sigma_nt,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test);

   private:
    /*!
     * @brief evaluate GPTS forces and linearization at this gp
     *
     * @tparam dim  dimension of the problem
     * @param[in] sele  slave side mortar element
     * @param[in] mele  master side mortar element
     * @param[in] sval  slave side shape function evaluated at current Gauss point
     * @param[in] sderiv  slave side shape function derivative at current Gauss point
     * @param[in] dsxi    directional derivative of slave side Gauss point coordinates
     * @param[in] mval    master side shape function evaluated at current Gauss point
     * @param[in] mderiv  master side shape function derivative at current Gauss point
     * @param[in] dmxi    directional derivative of master side Gauss point coordinates
     * @param[in] jac            Jacobian determinant of integration cell
     * @param[in] jacintcellmap  directional derivative of cell Jacobian
     * @param[in] wgt     Gauss point weight
     * @param[in] gap     gap
     * @param[in] dgapgp  directional derivative of gap
     * @param[in] gpn         Gauss point normal
     * @param[in] dnmap_unit  directional derivative of Gauss point normal
     * @param[in] sxi         slave side Gauss point coordinates
     * @param[in] mxi         master side Gauss point coordinates
     */
    template <int dim>
    void GPTSForces(MORTAR::Element& sele, MORTAR::Element& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dsxi,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dmxi, double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, double wgt, double gap,
        const CORE::GEN::pairedvector<int, double>& dgapgp, const double* gpn,
        std::vector<CORE::GEN::pairedvector<int, double>>& deriv_contact_normal, double* sxi,
        double* mxi);

   protected:
    template <int dim>
    void BuildAdjointTest(MORTAR::Element& moEle, double fac,
        const CORE::LINALG::SerialDenseMatrix& dsntdd,
        const CORE::LINALG::SerialDenseMatrix& d2sntdd2,
        const CORE::LINALG::SerialDenseMatrix& d2sntDdDn,
        const CORE::LINALG::SerialDenseMatrix& d2sntDdDt,
        const CORE::LINALG::SerialDenseMatrix& d2sntDdDpxi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& boundary_gpcoord_lin,
        CORE::LINALG::Matrix<dim, dim> derivtravo_slave,
        const std::vector<CORE::GEN::pairedvector<int, double>>& normal_deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& direction_deriv,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test);

    template <int dim>
    void IntegrateTest(double fac, MORTAR::Element& ele,
        const CORE::LINALG::SerialDenseVector& shape, const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, double wgt, double test_val,
        const CORE::GEN::pairedvector<int, double>& test_deriv,
        const CORE::LINALG::Matrix<dim, 1>& test_dir,
        const std::vector<CORE::GEN::pairedvector<int, double>>& test_dir_deriv);

    template <int dim>
    void IntegrateAdjointTest(double fac, double jac,
        const CORE::GEN::pairedvector<int, double>& jacintcellmap, double wgt, double test,
        const CORE::GEN::pairedvector<int, double>& deriv_test, MORTAR::Element& moEle,
        CORE::LINALG::SerialDenseVector& adjoint_test,
        CORE::GEN::pairedvector<int, CORE::LINALG::SerialDenseVector>& deriv_adjoint_test);

    //! nitsche theta
    double theta_;
    //! nitsche theta 2
    double theta_2_;
    //! type of nitsche weighting
    INPAR::CONTACT::NitscheWeighting nit_wgt_;
    //! nitsche penalty parameter in normal direction
    double ppn_;
    //! nitsche penalty parameter in tangential direction
    double ppt_;
    //! coulomb friction coefficient
    double frcoeff_;
    //! tresca friction bound
    double frbound_;
    //! type of contact friction law
    INPAR::CONTACT::FrictionType frtype_;
    //! time step size
    double dt_;
    //! flag indicating if unbiased two half pass algorithm is activated or not
    bool two_half_pass_;
  };

  namespace UTILS
  {
    //! map local surface coordinate to local parent coordinate
    template <int dim>
    void MapGPtoParent(MORTAR::Element& moEle, double* boundary_gpcoord, double wgt,
        CORE::LINALG::Matrix<dim, 1>& pxsi, CORE::LINALG::Matrix<dim, dim>& derivtravo_slave);

    //! templated version of MapGPtoParent
    template <CORE::FE::CellType parentdistype, int dim>
    void inline SoEleGP(MORTAR::Element& sele, double wgt, const double* gpcoord,
        CORE::LINALG::Matrix<dim, 1>& pxsi, CORE::LINALG::Matrix<dim, dim>& derivtrafo);

    //! actually not the velocity but the displacement increment
    template <int dim>
    void RelVelInvariant(MORTAR::Element& sele, const double* sxi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& derivsxi,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseMatrix& sderiv,
        MORTAR::Element& mele, const double* mxi,
        const std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi,
        const CORE::LINALG::SerialDenseVector& mval, const CORE::LINALG::SerialDenseMatrix& mderiv,
        const double& gap, const CORE::GEN::pairedvector<int, double>& deriv_gap,
        CORE::LINALG::Matrix<dim, 1>& relVel,
        std::vector<CORE::GEN::pairedvector<int, double>>& relVel_deriv, double fac = 1.0);

    template <int dim>
    void RelVel(MORTAR::Element& ele, const CORE::LINALG::SerialDenseVector& shape,
        const CORE::LINALG::SerialDenseMatrix& deriv,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dxi, double fac,
        CORE::LINALG::Matrix<dim, 1>& relVel,
        std::vector<CORE::GEN::pairedvector<int, double>>& relVel_deriv);

    template <int dim>
    void VectorScalarProduct(const CORE::LINALG::Matrix<dim, 1>& v1,
        const std::vector<CORE::GEN::pairedvector<int, double>>& v1d,
        const CORE::LINALG::Matrix<dim, 1>& v2,
        const std::vector<CORE::GEN::pairedvector<int, double>>& v2d, double& val,
        CORE::GEN::pairedvector<int, double>& val_deriv);

    void BuildTangentVectors3D(const double* np,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dn, double* t1p,
        std::vector<CORE::GEN::pairedvector<int, double>>& dt1, double* t2p,
        std::vector<CORE::GEN::pairedvector<int, double>>& dt2);

    template <int dim>
    void BuildTangentVectors(const double* np,
        const std::vector<CORE::GEN::pairedvector<int, double>>& dn, double* t1p,
        std::vector<CORE::GEN::pairedvector<int, double>>& dt1, double* t2p,
        std::vector<CORE::GEN::pairedvector<int, double>>& dt2);

    /*!
     * @brief determine weights (harmonic/slave/master) and scale penalty
     *
     * @param[in] sele     slave-side nitsche weight
     * @param[in] mele     master-side nitsche weight
     * @param[in] nit_wgt  weighting methood for nitsche contact
     * @param[in] dt       time step size
     * @param[out] ws      slave-side nitsche weight
     * @param[out] wm      master-side nitsche weight
     * @param[out] pen     scaled nitsche penalty parameter in normal direction
     * @param[out] pet     scaled nitsche penalty parameter in tangential direction
     */
    void NitscheWeightsAndScaling(MORTAR::Element& sele, MORTAR::Element& mele,
        INPAR::CONTACT::NitscheWeighting nit_wgt, double dt, double& ws, double& wm, double& pen,
        double& pet);
  }  // namespace UTILS
}  // namespace CONTACT
BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NITSCHE_INTEGRATOR_H
