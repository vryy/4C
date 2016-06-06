/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_integrator.H

\brief A class to perform integrations of Mortar matrices on the overlap
       of two MortarElements in 1D and 2D (derived version for
       augmented contact)

\level 2

\maintainer Michael Hiermeier

\date Apr 28, 2014

*/
/*---------------------------------------------------------------------*/

#ifndef CONTACT_AUGMENTED_INTEGRATOR_H_
#define CONTACT_AUGMENTED_INTEGRATOR_H_

# include "../drt_contact/contact_integrator.H"
# include "../drt_fem_general/drt_utils_integration.H"

namespace CONTACT
{
  class AugmentedIntegrator : public CoIntegrator
  {
  public:
    //! Default constructor without active set information
    AugmentedIntegrator(Teuchos::ParameterList& params,
        DRT::Element::DiscretizationType eletype,
        const Epetra_Comm& comm);

    //! Destructor
    virtual ~AugmentedIntegrator() {}

    //! @name 2D augmented Lagrange integration methods
    //! @{
    //! [derived]
    virtual void IntegrateDerivSegment2D(
         MORTAR::MortarElement& sele,
         double& sxia,
         double& sxib,
         MORTAR::MortarElement& mele,
         double& mxia,
         double& mxib,
         const Epetra_Comm& comm,
         const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! [derived]
    virtual void IntegrateDerivEle2D(
        MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles,
        bool *boundary_ele,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! Build remaining integrals and linearizations on the 1D slave interface
    void IntegrateDerivSlEle2D(
         MORTAR::MortarElement& sele,
         const Epetra_Comm& comm,
         const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    void IntegrateDerivSlEle2D(
         MORTAR::MortarElement& sele,const Epetra_Comm& comm,
         const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! @}

    //! @name 3D augmented Lagrange integration methods
    //! @{
    //! [derived]
    virtual void IntegrateDerivCell3DAuxPlane(
         MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
         Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
         const Epetra_Comm& comm,
         const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! [derived]
    virtual void IntegrateDerivEle3D(
         MORTAR::MortarElement& sele,
         std::vector<MORTAR::MortarElement*> meles,
         bool *boundary_ele, bool *proj_,
         const Epetra_Comm& comm,
         const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! Build the remaining integrals and linearizations on the 2D slave interface
    void IntegrateDerivSlEle3D(MORTAR::MortarElement& sele,
        const Epetra_Comm& comm,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    void IntegrateDerivSlEle3D(MORTAR::MortarElement& sele,
        const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);
    //! @}

  private:
    //! evaluate averaged weighted gap at gp
    void inline GP_2D_wG(
         MORTAR::MortarElement& sele,
         MORTAR::MortarElement& mele,
         LINALG::SerialDenseVector& sval,
         LINALG::SerialDenseVector& mval,
         LINALG::SerialDenseVector& lmval,
         LINALG::SerialDenseMatrix& scoord,
         LINALG::SerialDenseMatrix& mcoord,
         double* awgap,double* gpn,
         double* lengthn,
         double& dxdsxi,double& dsxideta,
         double& wgt);

    //! evaluate the scaling factor kappa at gp
    void inline GP_kappa(
         MORTAR::MortarElement& sele,
         LINALG::SerialDenseVector& lmval,
         double& wgt, double& jac);

    //! evaluate the weighted element Area at gp
    void inline GP_AugA(
        int& it,
        MORTAR::MortarElement& sele,
        LINALG::SerialDenseVector& sval,
        LINALG::SerialDenseVector& lmval,
        double& wgt,
        double& jac);

    //! evaluate the normal and the derivative of the normal at gp
    void inline GP_Normal_DerivNormal(
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector sval,
        LINALG::SerialDenseMatrix sderiv,
        const std::vector<GEN::pairedvector<int,double> >& dsxigp,
        double* gpn,
        std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
        int& linsize);

    /*! \brief evaluate the variation of the weighted gap (inner integral)
     *
     * ---> SLAVE and MASTER side */
    void inline GP_VarWGap(
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval,
        LINALG::SerialDenseVector& mval,
        LINALG::SerialDenseVector& lmval,
        double* gpn,
        double& wgt, double& jac);

    /*! evaluate the lin. of the scaling factor kappa at gp (2-D)
     *  (segment based integration) */
    void inline GP_2D_kappa_Lin(
        int& iter,
        MORTAR::MortarElement& sele,
        LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseMatrix& lmderiv,
        double& dsxideta, double& dxdsxi,
        double& dxdsxidsxi,
        double& wgt,
        const GEN::pairedvector<int,double>& dsxigp,
        const GEN::pairedvector<int,double>& derivjac,
        const std::vector<GEN::pairedvector<int,double> >& ximaps);

    /*! \brief evaluate the lin. of the scaling factor kappa at gp (2-D)
     *  (Element based integration) */
    void inline GP_2D_kappa_Ele_Lin(
        int& iter,
        MORTAR::MortarElement& sele,
        LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseMatrix& lmderiv,
        double& dxdsxi,double& wgt,
        const GEN::pairedvector<int,double>& derivjac);

    /*!
    \brief evaluate the lin. of the scaling factor kappa at gp (3-D)

    */
    void inline GP_3D_kappa_Lin(
        int& iter,
        MORTAR::MortarElement& sele,
        LINALG::SerialDenseVector& lmval,
        LINALG::SerialDenseMatrix& lmderiv,
        double& wgt, double& jac,
        const std::vector<GEN::pairedvector<int,double> >& dsxigp,
        const GEN::pairedvector<int,double>& jacintcellmap);

    //! evaluate the linearization of weighted element area at gp
    void inline GP_AugA_Lin(
        int& it,
        MORTAR::MortarElement& sele,
        LINALG::SerialDenseVector& sval,
        LINALG::SerialDenseVector& lmval,
        double& wgt,
        double& jac,
        const GEN::pairedvector<int,double>& derivjac);

    /*! \brief evaluate linearization of the variation of the weighted gap (2-D)
     *
     *  ---> SLAVE and MASTER side */
    void inline GP_2D_VarWGap_Lin(
        int& iter,
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval,
        LINALG::SerialDenseVector& mval,
        LINALG::SerialDenseVector& lmval,
        double* gpn,
        LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv,
        LINALG::SerialDenseMatrix& lmderiv,
        double& dsxideta, double& dxdsxi,
        double& dxdsxidsxi,
        double& wgt,
        const GEN::pairedvector<int,double>& dsxigp,
        const GEN::pairedvector<int,double>& dmxigp,
        const GEN::pairedvector<int,double>& derivjac,
        const std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
        const std::vector<GEN::pairedvector<int,double> >& ximaps);

    /*! \brief evaluate linearization of the variation of the weighted gap (2-D)
     *
     *  ---> SLAVE and MASTER side
     *  (element based integration) */
    void inline GP_2D_VarWGap_Ele_Lin(
        int& iter,
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval,
        LINALG::SerialDenseVector& mval,
        LINALG::SerialDenseVector& lmval,
        double* gpn,
        LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv,
        LINALG::SerialDenseMatrix& lmderiv,
        double& dxdsxi, double& wgt,
        const GEN::pairedvector<int,double>& dmxigp,
        const GEN::pairedvector<int,double>&derivjac,
        const std::vector<GEN::pairedvector<int,double> >& dnmap_unit);

    /*! \brief evaluate linearization of the variation of the weighted gap (3-D)
     *
     * ---> SLAVE and MASTER side */
    void inline GP_3D_VarWGap_Lin(
        int& iter,
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        LINALG::SerialDenseVector& sval,
        LINALG::SerialDenseVector& mval,
        LINALG::SerialDenseVector& lmval,
        double* gpn,
        LINALG::SerialDenseMatrix& sderiv,
        LINALG::SerialDenseMatrix& mderiv,
        LINALG::SerialDenseMatrix& lmderiv,
        double& wgt, double& jac,
        const std::vector<GEN::pairedvector<int,double> >& dsxigp,
        const std::vector<GEN::pairedvector<int,double> >& dmxigp,
        const GEN::pairedvector<int,double>& jacintcellmap,
        const std::vector<GEN::pairedvector<int,double> >& dnmap_unit);
    //@}

    //! @name Return methods
    //! @{
    //! Return the solution type
    const enum INPAR::CONTACT::SolvingStrategy& SolType() const { return stype_; };

    //@}

  protected:
    // don't want = operator and cctor
    AugmentedIntegrator operator = (const AugmentedIntegrator& old);
    AugmentedIntegrator(const AugmentedIntegrator& old);

  }; // class AugmentedIntegrator
} // namespace CONTACT

#endif /* CONTACT_AUGMENTED_INTEGRATOR_H_ */
