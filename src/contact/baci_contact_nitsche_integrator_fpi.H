/*---------------------------------------------------------------------*/
/*! \file

\brief A class to perform integrations of nitsche related terms for the fpi contact case

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_NITSCHE_INTEGRATOR_FPI_H
#define BACI_CONTACT_NITSCHE_INTEGRATOR_FPI_H

#include "baci_contact_nitsche_integrator_poro.H"
#include "baci_linalg_fixedsizematrix.H"

namespace XFEM
{
  class XFluid_Contact_Comm;
}

namespace CONTACT
{
  class CoElement;

  class CoIntegratorNitscheFpi : public CoIntegratorNitschePoro
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
    CoIntegratorNitscheFpi(Teuchos::ParameterList& params, DRT::Element::DiscretizationType eletype,
        const Epetra_Comm& comm);
    //! @name Derived functions
    //! @{

    //! @name currently unsupported derived methods
    //! @{
    void IntegrateDerivSegment2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
        MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override
    {
      dserror("Segment based integration is currently unsupported!");
    }

    void IntegrateDerivEle2D(MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
        bool* boundary_ele, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override
    {
      dserror("Element based integration in 2D is currently unsupported!");
    }

    void IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override
    {
      dserror("The auxiliary plane 3-D coupling integration case is currently unsupported!");
    }
    //! @}

    /*!
     \brief First, reevaluate which gausspoints should be used
     Second, Build all integrals and linearizations without segmentation -- 3D
     (i.e. M, g, LinM, Ling and possibly D, LinD)
     */
    void IntegrateDerivEle3D(MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
        bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override;

    //! @}

   protected:
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
        std::vector<CORE::GEN::pairedvector<int, double>>& derivmxi) override
    {
      dserror("2d problems not available for CoIntegratorNitscheFsi, as CutFEM is only 3D!");
    }

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


    template <int dim>
    double GetNormalContactTransition(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const CORE::LINALG::SerialDenseVector& sval, const CORE::LINALG::SerialDenseVector& mval,
        const double* sxi, const CORE::LINALG::Matrix<dim, 1>& pxsi,
        const CORE::LINALG::Matrix<dim, 1>& normal, bool& FSI_integrated, bool& gp_on_this_proc);

    /// Update Element contact state -2...not_specified, -1...no_contact, 0...mixed, 1...contact
    void UpdateEleContactState(MORTAR::MortarElement& sele, int state);

    /// Element contact state -2...not_specified, -1...no_contact, 0...mixed, 1...contact
    int ele_contact_state_;

    /// Xfluid Contact Communicator
    Teuchos::RCP<XFEM::XFluid_Contact_Comm> xf_c_comm_;
  };
}  // namespace CONTACT
#endif  // CONTACT_NITSCHE_INTEGRATOR_FPI_H
