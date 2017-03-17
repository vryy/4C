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

#include "../drt_contact/contact_integrator.H"
#include "../drt_contact_aug/contact_integrator_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace CONTACT {
namespace AUG {

  // forward declaration
  class IntegratorGeneric;

  /*--------------------------------------------------------------------------*/
  /** \brief Intermediate class necessary to handle a vector of master elements
   *  with arbitrary shapes
   *
   *  \author hiermeier \date 03/17 */
  class IntegrationWrapper : public CoIntegrator
  {
    friend class IntegratorGeneric;

  public:
    /// constructor
    IntegrationWrapper( Teuchos::ParameterList& params,
        DRT::Element::DiscretizationType slavetype,
        const Epetra_Comm& comm );

    /// destructor
    virtual ~IntegrationWrapper() {};

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

    //! @}

    //! Build remaining integrals and linearizations on the 1-D/2-D slave interface
    void IntegrateDerivSlaveElement(
         MORTAR::MortarElement& sele,
         const Epetra_Comm& comm,
         const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    void IntegrateDerivSlaveElement(
         MORTAR::MortarElement& sele,const Epetra_Comm& comm,
         const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

  private:
    Teuchos::RCP<IntegratorGeneric> integrator_;

    static INTEGRATOR::UniqueProjInfoPair projInfo_;
  };  // class IntegratorWrapper


  /*--------------------------------------------------------------------------*/
  class IntegratorGeneric
  {
  public:
    static Teuchos::RCP<CONTACT::AUG::IntegratorGeneric> Create(
        int probdim,
        DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype,
        const CONTACT::ParamsInterface & cparams,
        CONTACT::CoIntegrator* wrapper );

  private:
    static CONTACT::AUG::IntegratorGeneric* Create2D(
        DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype,
        const CONTACT::ParamsInterface & cparams,
        CONTACT::CoIntegrator* wrapper );

    template < DRT::Element::DiscretizationType slavetype >
    static CONTACT::AUG::IntegratorGeneric * Create2D(
        DRT::Element::DiscretizationType mastertype,
        const CONTACT::ParamsInterface & cparams,
        CONTACT::CoIntegrator* wrapper );

    static CONTACT::AUG::IntegratorGeneric* Create3D(
        DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype,
        const CONTACT::ParamsInterface & cparams,
        CONTACT::CoIntegrator* wrapper );

    template < DRT::Element::DiscretizationType slavetype >
    static CONTACT::AUG::IntegratorGeneric * Create3D(
        DRT::Element::DiscretizationType mastertype,
        const CONTACT::ParamsInterface & cparams,
        CONTACT::CoIntegrator* wrapper );

  public:
    IntegratorGeneric(
        const CONTACT::ParamsInterface& cparams,
        CONTACT::CoIntegrator& wrapper)
        : cparams_(cparams),
          wrapper_( static_cast<IntegrationWrapper&>( wrapper ) )
    {};

    virtual ~IntegratorGeneric() {};

    //! @name 2D augmented Lagrange integration methods
    //! @{

    virtual void IntegrateDerivSegment2D(
       MORTAR::MortarElement& sele,
       double& sxia,
       double& sxib,
       MORTAR::MortarElement& mele,
       double& mxia,
       double& mxib ) = 0;

    virtual void IntegrateDerivEle2D(
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo ) = 0;

    //! @}

    //! @name 3D augmented Lagrange integration methods
    //! @{
    //! [derived]
    virtual void IntegrateDerivCell3DAuxPlane(
         MORTAR::MortarElement& sele,
         MORTAR::MortarElement& mele,
         MORTAR::IntCell& cell,
         double* auxn ) = 0;

    //! [derived]
    virtual void IntegrateDerivEle3D(
         MORTAR::MortarElement& sele,
         MORTAR::MortarElement& mele,
         bool boundary_ele,
         const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo ) = 0;

    //! @}

    //! Build the remaining integrals and linearizations on the 1-D/2-D slave interface
    virtual void IntegrateDerivSlaveElement( MORTAR::MortarElement& sele ) = 0;

  protected:
    INPAR::MORTAR::ShapeFcn ShapeFcn()
    {
      return wrapper_.shapefcn_;
    }

  protected:
    const CONTACT::ParamsInterface & cparams_;

    CONTACT::AUG::IntegrationWrapper & wrapper_;

  };  //

  /*--------------------------------------------------------------------------*/
  template < unsigned probdim,
             DRT::Element::DiscretizationType slavetype,
             DRT::Element::DiscretizationType mastertype,
             unsigned slavedim = DRT::UTILS::DisTypeToDim<slavetype>::dim,
             unsigned slavenumnode = DRT::UTILS::DisTypeToNumNodePerEle<slavetype>::numNodePerElement,
             unsigned masterdim = DRT::UTILS::DisTypeToDim<mastertype>::dim,
             unsigned masternumnode = DRT::UTILS::DisTypeToNumNodePerEle<mastertype>::numNodePerElement >
  class Integrator : public IntegratorGeneric
  {
  public:
    //! Default constructor without active set information
    Integrator(
        const CONTACT::ParamsInterface & cparams,
        CONTACT::CoIntegrator& wrapper );

    //! Destructor
    virtual ~Integrator() {}

  protected:
    //! @name 2D augmented Lagrange integration methods
    //! @{

    virtual void IntegrateDerivSegment2D(
       MORTAR::MortarElement& sele,
       double& sxia,
       double& sxib,
       MORTAR::MortarElement& mele,
       double& mxia,
       double& mxib );

    virtual void IntegrateDerivEle2D(
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo );

    //! @}

    //! @name 3D augmented Lagrange integration methods
    //! @{
    //! [derived]
    virtual void IntegrateDerivCell3DAuxPlane(
         MORTAR::MortarElement& sele,
         MORTAR::MortarElement& mele,
         MORTAR::IntCell& cell,
         double* auxn );

    //! [derived]
    virtual void IntegrateDerivEle3D(
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        bool boundary_ele,
        const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo );

    //! @}

    //! Build the remaining integrals and linearizations on the 1-D/2-D slave interface
    virtual void IntegrateDerivSlaveElement( MORTAR::MortarElement& sele );

  private:
    //! evaluate averaged weighted gap at gp
    void GP_2D_wG(
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
    void GP_kappa(
         MORTAR::MortarElement& sele,
         const LINALG::Matrix<slavenumnode,1>& lmval,
         double wgt,
         double jac) const;

    //! evaluate the weighted element Area at gp
    void GP_AugA(
        int it,
        MORTAR::MortarElement& sele,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        double wgt,
        double jac) const;

    //! evaluate the normal and the derivative of the normal at gp
    void GP_Normal_DerivNormal(
        MORTAR::MortarElement& sele,
        const MORTAR::MortarElement& mele,
        const LINALG::Matrix<slavenumnode,1>& sval,
        const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
        const std::vector<GEN::pairedvector<int,double> >& dsxigp,
        double* gpn,
        std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
        int linsize);

    /*! \brief evaluate the variation of the weighted gap (inner integral)
     *
     * ---> SLAVE and MASTER side */
    void GP_VarWGap(
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        const LINALG::Matrix<slavenumnode,1>& sval,
        const LINALG::Matrix<masternumnode,1>& mval,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const double* gpn,
        double wgt,
        double jac) const;

    /*! evaluate the lin. of the scaling factor kappa at gp (2-D)
     *  (segment based integration) */
    void GP_2D_kappa_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
        double dsxideta,
        double dxdsxi,
        double dxdsxidsxi,
        double wgt,
        const GEN::pairedvector<int,double>& dsxigp,
        const GEN::pairedvector<int,double>& derivjac,
        const std::vector<GEN::pairedvector<int,double> >& ximaps);

    /*! \brief evaluate the lin. of the scaling factor kappa at gp (2-D)
     *  (Element based integration) */
    void GP_2D_kappa_Ele_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
        double dxdsxi,
        double wgt,
        const GEN::pairedvector<int,double>& derivjac);

    /*!
    \brief evaluate the lin. of the scaling factor kappa at gp (3-D)

    */
    void GP_3D_kappa_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
        double wgt,
        double jac,
        const std::vector<GEN::pairedvector<int,double> >& dsxigp,
        const GEN::pairedvector<int,double>& jacintcellmap);

    //! evaluate the linearization of weighted element area at gp
    void GP_AugA_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        double wgt,
        double jac,
        const GEN::pairedvector<int,double>& derivjac) const;

    /*! \brief evaluate linearization of the variation of the weighted gap (2-D)
     *
     *  ---> SLAVE and MASTER side */
    void GP_2D_VarWGap_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        const LINALG::Matrix<slavenumnode,1>& sval,
        const LINALG::Matrix<masternumnode,1>& mval,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const double* gpn,
        const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
        const LINALG::Matrix<masternumnode,masterdim>& mderiv,
        const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
        double dsxideta,
        double dxdsxi,
        double dxdsxidsxi,
        double wgt,
        const GEN::pairedvector<int,double>& dsxigp,
        const GEN::pairedvector<int,double>& dmxigp,
        const GEN::pairedvector<int,double>& derivjac,
        const std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
        const std::vector<GEN::pairedvector<int,double> >& ximaps) const;

    /*! \brief evaluate linearization of the variation of the weighted gap (2-D)
     *
     *  ---> SLAVE and MASTER side
     *  (element based integration) */
    void GP_2D_VarWGap_Ele_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        const LINALG::Matrix<slavenumnode,1>& sval,
        const LINALG::Matrix<masternumnode,1>& mval,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const double* gpn,
        const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
        const LINALG::Matrix<masternumnode,masterdim>& mderiv,
        const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
        double dxdsxi,
        double wgt,
        const GEN::pairedvector<int,double>& dmxigp,
        const GEN::pairedvector<int,double>&derivjac,
        const std::vector<GEN::pairedvector<int,double> >& dnmap_unit) const;

    /*! \brief evaluate linearization of the variation of the weighted gap (3-D)
     *
     * ---> SLAVE and MASTER side */
    void GP_3D_VarWGap_Lin(
        unsigned iter,
        MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele,
        const LINALG::Matrix<slavenumnode,1>& sval,
        const LINALG::Matrix<masternumnode,1>& mval,
        const LINALG::Matrix<slavenumnode,1>& lmval,
        const double* gpn,
        const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
        const LINALG::Matrix<masternumnode,masterdim>& mderiv,
        const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
        double wgt,
        double jac,
        const std::vector<GEN::pairedvector<int,double> >& dsxigp,
        const std::vector<GEN::pairedvector<int,double> >& dmxigp,
        const GEN::pairedvector<int,double>& jacintcellmap,
        const std::vector<GEN::pairedvector<int,double> >& dnmap_unit) const;
    //@}


  protected:
    // don't want = operator and cctor
    Integrator operator = (const Integrator& old);
    Integrator(const Integrator& old);

  private:
    // vectors for shape fct. evaluation
    static LINALG::Matrix<slavenumnode,1> sval_;
    static LINALG::Matrix<slavenumnode,slavedim> sderiv_;
    static LINALG::Matrix<masternumnode,1> mval_;
    static LINALG::Matrix<masternumnode,masterdim> mderiv_;
    static LINALG::Matrix<slavenumnode,1> lmval_;
    static LINALG::Matrix<slavenumnode,slavedim> lmderiv_;

    // create empty vectors for shape fct. evaluation
    static LINALG::Matrix<slavenumnode,3> ssecderiv_;

    // slave and master nodal coords for Jacobian / GP evaluation
    static LINALG::Matrix<3,slavenumnode> scoord_;
    static LINALG::Matrix<3,masternumnode> mcoord_;

    // directional derivatives of sxia, sxib, mxia, mxib
    static std::vector<GEN::pairedvector<int,double> > ximaps_;

    // deriv of x and y comp. of gpn [a.k.a. GaussPoint normal] (unit)
    static std::vector<GEN::pairedvector<int,double> > dnmap_unit_;

    // GP slave coordinate derivatives
    static std::vector<GEN::pairedvector<int,double> > dsxigp_;

    // GP master coordinate derivatives
    static std::vector<GEN::pairedvector<int,double> > dmxigp_;

    // Jacobian derivative
    static GEN::pairedvector<int,double> derivjac_;

    static GEN::pairedvector<int,LINALG::Matrix<3,1> > lingp_;

    // directional derivative of slave GP normal (non-unit)
    static GEN::pairedvector<int,double> dmap_nxsl_gp_;
    static GEN::pairedvector<int,double> dmap_nysl_gp_;
    static GEN::pairedvector<int,double> dmap_nzsl_gp_;

  }; // class AugmentedIntegrator

} // namespace AUG
} // namespace CONTACT

/*----------------------------------------------------------------------------*/
// static members in CONTACT::AUG::Integrator
template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<slavenumnode,1>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::sval_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<slavenumnode,slavedim>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::sderiv_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<masternumnode,1>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::mval_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<masternumnode,masterdim>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::mderiv_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<slavenumnode,1>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::lmval_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<slavenumnode,slavedim>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::lmderiv_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<slavenumnode,3>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::ssecderiv_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<3,slavenumnode>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::scoord_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > LINALG::Matrix<3,masternumnode>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::mcoord_;

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > std::vector<GEN::pairedvector<int,double> >
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::ximaps_( 0, GEN::pairedvector<int,double>(0) );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > std::vector<GEN::pairedvector<int,double> >
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::dnmap_unit_( 0, GEN::pairedvector<int,double>(0) );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > std::vector<GEN::pairedvector<int,double> >
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::dsxigp_( 0, GEN::pairedvector<int,double>(0) );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > std::vector<GEN::pairedvector<int,double> >
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::dmxigp_( 0, GEN::pairedvector<int,double>(0) );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > GEN::pairedvector<int,double>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::derivjac_( 0 );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > GEN::pairedvector<int,LINALG::Matrix<3,1> >
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::lingp_( 0 );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > GEN::pairedvector<int,double>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::dmap_nxsl_gp_( 0 );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > GEN::pairedvector<int,double>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::dmap_nysl_gp_( 0 );

template < unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    unsigned slavedim, unsigned slavenumnode,
    unsigned masterdim, unsigned masternumnode > GEN::pairedvector<int,double>
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::dmap_nzsl_gp_( 0 );


#endif /* CONTACT_AUGMENTED_INTEGRATOR_H_ */
