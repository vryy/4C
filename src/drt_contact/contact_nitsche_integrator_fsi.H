/*---------------------------------------------------------------------*/
/*! \file

\brief A class to perform integrations of nitsche related terms for the fsi contact case

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_INTEGRATOR_NITSCHE_FSI_H
#define CONTACT_INTEGRATOR_NITSCHE_FSI_H

#include "contact_nitsche_integrator.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace XFEM
{
  class XFluid_Contact_Comm;
}

namespace CONTACT
{
  class CoElement;

  class CoIntegratorNitscheFsi : public CoIntegratorNitsche
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
    CoIntegratorNitscheFsi(Teuchos::ParameterList& params, DRT::Element::DiscretizationType eletype,
        const Epetra_Comm& comm);

    /*!
     \brief Destructor

     */
    virtual ~CoIntegratorNitscheFsi() {}

    //! @name Derived functions
    //! @{

    //! @name currently unsupported derived methods
    //! @{
    virtual void IntegrateDerivSegment2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
        MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
    {
      dserror("Segment based integration is currently unsupported!");
    }

    virtual void IntegrateDerivEle2D(MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
    {
      dserror("Element based integration in 2D is currently unsupported!");
    }

    virtual void IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
        const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
    {
      dserror("The auxiliary plane 3-D coupling integration case is currently unsupported!");
    }
    //! @}

    /*!
     \brief First, reevaluate which gausspoints should be used
     Second, Build all integrals and linearizations without segmentation -- 3D
     (i.e. M, g, LinM, Ling and possibly D, LinD)
     */
    virtual void IntegrateDerivEle3D(MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj_,
        const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! @}

   protected:
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
        std::vector<GEN::pairedvector<int, double>>& derivmxi)
    {
      dserror("2d problems not available for CoIntegratorNitscheFsi, as CutFEM is only 3D!");
      return;
    }

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

    /// Update Element contact state -2...not_specified, -1...no_contact, 0...mixed, 1...contact
    void UpdateEleContactState(MORTAR::MortarElement& sele, int state);

    /// Element contact state -2...not_specified, -1...no_contact, 0...mixed, 1...contact
    int ele_contact_state_;

    /// Xfluid Contact Communicator
    Teuchos::RCP<XFEM::XFluid_Contact_Comm> xf_c_comm_;
  };

  namespace UTILS
  {
    /// Compute Cauchy stress component sigma_{n dir} at local coord xsi
    double SolidCauchyAtXi(CONTACT::CoElement* cele,  // the contact element
        const LINALG::Matrix<2, 1>& xsi,              // local coord on the ele element
        const LINALG::Matrix<3, 1>& n,                // normal n
        const LINALG::Matrix<3, 1>& dir               // second directional vector
    );

    /// templated version of MapGPtoParent
    template <DRT::Element::DiscretizationType parentdistype, int dim>
    void inline SoEleGP(MORTAR::MortarElement& sele, const double wgt, const double* gpcoord,
        LINALG::Matrix<dim, 1>& pxsi, LINALG::Matrix<dim, dim>& derivtrafo);
  }  // namespace UTILS
}  // namespace CONTACT
#endif  // CONTACT_INTEGRATOR_NITSCHE_FSI_H
