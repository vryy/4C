/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_integrator.cpp

\brief A class to perform integrations of Mortar matrices on the overlap
       of two MortarElements in 1D and 2D (derived version for
       augmented contact)

\level 2

\maintainer Michael Hiermeier

\date Apr 28, 2014

*/
/*---------------------------------------------------------------------*/
#include "contact_augmented_integrator.H"
#include "contact_integrator_utils.H"

#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_serialdensevector.H"
#include <Epetra_Map.h>

// define and initialize static member
CONTACT::INTEGRATOR::UniqueProjInfoPair CONTACT::AUG::IntegrationWrapper::projInfo_( 0 );

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegrationWrapper::IntegrationWrapper(
    Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype,
    const Epetra_Comm& comm)
    : CONTACT::CoIntegrator::CoIntegrator(params,eletype,comm)

{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele,
     MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell,
     double* auxn,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  dsassert( Dim()==3, "ERROR: 3D integration method called for non-3D problem" );

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  if (ShapeFcn() == INPAR::MORTAR::shape_dual || ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane supports no Dual shape functions for the "
        "augmented Lagrange solving strategy!");

  integrator_ = IntegratorGeneric::Create( Dim(), sele.Shape(), mele.Shape(), *cparams_ptr, this );
  integrator_->IntegrateDerivCell3DAuxPlane( sele, mele, *cell, auxn );

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivEle3D(
     MORTAR::MortarElement& sele,
     std::vector<MORTAR::MortarElement*> meles,
     bool *boundary_ele,
     bool *proj_,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // explicitly defined shape function type needed
  if ( ShapeFcn() == INPAR::MORTAR::shape_undefined )
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  dsassert( Dim()==3, "ERROR: 3D integration method called for non-3D problem" );

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes)
    dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

  // check input data
  for (int test=0;test<(int)meles.size();++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
      dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  }

  // contact with wear
  if(wearlaw_!= INPAR::WEAR::wear_none)
    dserror( "Wear is not supported!" );

  //Boundary Segmentation check -- HasProj()-check
  *boundary_ele = BoundarySegmCheck3D(sele,meles);

  *proj_ = INTEGRATOR::FindFeasibleMasterElement3D( sele, meles, boundary_ele, *this, projInfo_ );

  for ( INTEGRATOR::UniqueProjInfoPair::const_iterator cit = projInfo_.begin();
        cit != projInfo_.end(); ++cit )
  {
    MORTAR::MortarElement& mele = *cit->first;
    integrator_ = IntegratorGeneric::Create( Dim(), sele.Shape(), mele.Shape(),
        *cparams_ptr, this );
    integrator_->IntegrateDerivEle3D( sele, mele, *boundary_ele, cit->second );
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivSlaveElement(
    MORTAR::MortarElement& sele,
    const Epetra_Comm& comm,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr);
  if (cparams_ptr.is_null())
    dserror("Cast to CONTACT::ParamsInterface failed!");
  IntegrateDerivSlaveElement(sele,comm,cparams_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivSlaveElement(
     MORTAR::MortarElement& sele,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  integrator_ = IntegratorGeneric::Create( Dim(), sele.Shape(), sele.Shape(),
      *cparams_ptr, this );
  integrator_->IntegrateDerivSlaveElement( sele );

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivSegment2D(
     MORTAR::MortarElement& sele,
     double& sxia,
     double& sxib,
     MORTAR::MortarElement& mele,
     double& mxia,
     double& mxib,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // *********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape()==MORTAR::MortarElement::line3 || ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("ERROR: Petrov-Galerkin / quadratic FE interpolation not yet implemented.");

  //check for problem dimension
  dsassert( Dim()==2, "ERROR: 2D integration method called for non-2D problem" );

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

  integrator_ = IntegratorGeneric::Create( Dim(), sele.Shape(), mele.Shape(),
      *cparams_ptr, this );
  integrator_->IntegrateDerivSegment2D( sele, sxia, sxib, mele, mxia, mxib );

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::IntegrationWrapper::IntegrateDerivEle2D(
    MORTAR::MortarElement& sele,
    std::vector<MORTAR::MortarElement*> meles,
    bool *boundary_ele,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // *********************************************************************
  // Check integrator input for non-reasonable quantities
  // *********************************************************************
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

// explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // check input data
  for (int i=0;i<(int)meles.size();++i)
  {
    if ((!sele.IsSlave()) || (meles[i]->IsSlave()))
      dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  }

  // number of nodes (slave) and problem dimension
  const int nrow = sele.NumNode();

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode)
      dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  //Boundary Segmentation check -- HasProj()-check
  if(IntegrationType()==INPAR::MORTAR::inttype_elements_BS)
    *boundary_ele=BoundarySegmCheck2D(sele,meles);

  if (*boundary_ele==false || IntegrationType()==INPAR::MORTAR::inttype_elements)
  {
    INTEGRATOR::FindFeasibleMasterElement2D( sele, meles, *this, projInfo_ );

    for ( INTEGRATOR::UniqueProjInfoPair::const_iterator cit = projInfo_.begin();
          cit != projInfo_.end(); ++cit )
    {
      MORTAR::MortarElement& mele = *cit->first;
      integrator_ = IntegratorGeneric::Create( Dim(), sele.Shape(), mele.Shape(),
          *cparams_ptr, this );

      integrator_->IntegrateDerivEle2D( sele, mele, cit->second );
    }
  }//boundary_ele check

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AUG::IntegratorGeneric>
    CONTACT::AUG::IntegratorGeneric::Create(
    int probdim,
    DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper )
{
  switch ( probdim )
  {
    case 2:
      return Teuchos::rcp( Create2D( slavetype, mastertype, cparams, wrapper ) );
    case 3:
      return Teuchos::rcp( Create3D( slavetype, mastertype, cparams, wrapper ) );
    default:
      dserror( "Unsupported problem dimension %d", probdim );
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create2D(
    DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper )
{
  switch ( slavetype )
  {
    case DRT::Element::line2:
      return Create2D<DRT::Element::line2>( mastertype, cparams, wrapper );
    default:
      dserror( "Unsupported master element type %s",
          DRT::DistypeToString( mastertype ).c_str() );
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < DRT::Element::DiscretizationType slavetype >
CONTACT::AUG::IntegratorGeneric * CONTACT::AUG::IntegratorGeneric::Create2D(
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper )
{
  switch ( mastertype )
  {
    case DRT::Element::line2:
      return new CONTACT::AUG::Integrator<2,slavetype,DRT::Element::line2>( cparams, *wrapper );
    default:
      dserror( "Unsupported master element type %s",
          DRT::DistypeToString( mastertype ).c_str() );
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::IntegratorGeneric* CONTACT::AUG::IntegratorGeneric::Create3D(
    DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper )
{
  switch ( slavetype )
  {
    case DRT::Element::line2:
      return Create3D<DRT::Element::line2>( mastertype, cparams, wrapper );
    case DRT::Element::quad4:
      return Create3D<DRT::Element::quad4>( mastertype, cparams, wrapper );
    case DRT::Element::tri3:
      return Create3D<DRT::Element::tri3>( mastertype, cparams, wrapper );
    default:
      dserror( "Unsupported master element type %s",
          DRT::DistypeToString( mastertype ).c_str() );
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < DRT::Element::DiscretizationType slavetype >
CONTACT::AUG::IntegratorGeneric * CONTACT::AUG::IntegratorGeneric::Create3D(
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper )
{
  switch ( mastertype )
  {
    case DRT::Element::line2:
      return new CONTACT::AUG::Integrator<3,slavetype,DRT::Element::line2>( cparams, *wrapper );
    case DRT::Element::quad4:
      return new CONTACT::AUG::Integrator<3,slavetype,DRT::Element::quad4>( cparams, *wrapper );
    case DRT::Element::tri3:
      return new CONTACT::AUG::Integrator<3,slavetype,DRT::Element::tri3>( cparams, *wrapper );
    default:
      dserror( "Unsupported master element type %s",
          DRT::DistypeToString( mastertype ).c_str() );
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::Integrator(
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator& wrapper )
    : CONTACT::AUG::IntegratorGeneric( cparams, wrapper )
{
  // empty constructor body
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::IntegrateDerivSegment2D(
       MORTAR::MortarElement& sele,
       double& sxia,
       double& sxib,
       MORTAR::MortarElement& mele,
       double& mxia,
       double& mxib )
{
  // *********************************************************************
  // Prepare integration
  // *********************************************************************

  // Do the augmented Lagrange stuff...
  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes)
    dserror("ERROR: AugmentedIntegrator::IntegrateDerivSegment2D: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  bool linlm = false;
  for ( unsigned k=0; k<slavenumnode; ++k )
  {
   MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
   if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
   bound += mymrtrnode->IsOnBound();
  }

  if ( bound )
   dserror("ERROR: AugmentedIntegrator::IntegrateDerivSegment2D: "
     "Boundary modification is not considered for the augmented Lagrange formulation!");

  if ( slavetype == DRT::Element::line3 )
   dserror("ERROR: AugmentedIntegrator::IntegrateDerivSegment2D: "
         "Quadratic shape functions are not yet considered.");

  // get slave and master nodal coords for Jacobian / GP evaluation
  sele.GetNodalCoords(scoord_);
  mele.GetNodalCoords(mcoord_);

  // nodal lagrange multiplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  lagmult   = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
  sele.GetNodalLagMult(*lagmult);

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  int linsize = 0;
  for ( unsigned i=0; i<slavenumnode; ++i )
  {
   CoNode* cnode = static_cast<CoNode*> (mynodes[i]);
   linsize += cnode->GetLinsize();
  }

  // *********************************************************************
  // Find out about whether start / end of overlap are slave or master!
  // CAUTION: be careful with positive rotation direction ("Umlaufsinn")
  // sxia -> belongs to sele.Nodes()[0]
  // sxib -> belongs to sele.Nodes()[1]
  // mxia -> belongs to mele.Nodes()[0]
  // mxib -> belongs to mele.Nodes()[1]
  // but slave and master have different positive rotation directions,
  // counter-clockwise for slave side, clockwise for master side!
  // this means that mxia belongs to sxib and vice versa!
  // *********************************************************************

  bool startslave = false;
  bool endslave = false;

  if ( sxia!=-1.0 and mxib!=1.0 )
   dserror("ERROR: First outer node is neither slave nor master node");
  if ( sxib!=1.0 and mxia!=-1.0 )
   dserror("ERROR: Second outer node is neither slave nor master node");

  if (sxia==-1.0)
    startslave = true;
  else
    startslave = false;

  if (sxib==1.0)
    endslave   = true;
  else
    endslave   = false;

  // get directional derivatives of sxia, sxib, mxia, mxib
  INTEGRATOR::ResetPairedVector<4>( linsize+probdim*masternumnode, ximaps_ );
  this->wrapper_.DerivXiAB2D(sele,sxia,sxib,mele,mxia,mxib,ximaps_,startslave,endslave,linsize);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<this->wrapper_.nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {this->wrapper_.Coordinate(gp,0), 0.0};
    double wgt = this->wrapper_.Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] = {0.0, 0.0};

    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    MORTAR::MortarProjector::Impl(sele,mele)->ProjectGaussPoint2D(sele,sxi,mele,mxi);

    // check GP projection
    if ( (mxi[0]<mxia) or (mxi[0]>mxib) )
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      std::cout << __FILE__ << ", line: " << __LINE__ << std::endl;
      std::cout << "EXCEPTION-WARNING: AugmentedIntegrator::IntegrateAvgWgapSegment2D: "
         "Gauss point projection failed! mxi=" << mxi[0] << std::endl;
      throw 100;
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (linlm)
    {
      sele.EvaluateShapeLagMultLin(this->ShapeFcn(),sxi,lmval_,lmderiv_,slavenumnode);
    }
    else
    {
      sele.EvaluateShapeLagMult(this->ShapeFcn(),sxi,lmval_,lmderiv_,slavenumnode,true);
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape( sxi, sval_, sderiv_, slavenumnode, false );
    mele.EvaluateShape( mxi, mval_, mderiv_, masternumnode, false );

    // evaluate the two slave side Jacobians
    const double dxdsxi = sele.Jacobian(sxi);
    const double dsxideta = -0.5*sxia + 0.5*sxib;

    // evaluate linearizations *******************************************
    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape( sxi, ssecderiv_, slavenumnode );

    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    static_cast<CONTACT::CoElement&>(sele).DJacDXi( djacdxi, sxi, ssecderiv_ );
    const double dxdsxidsxi=djacdxi[0]; // only 2D here

    // evaluate the GP slave coordinate derivatives
    INTEGRATOR::ResetPairedVector<1>( linsize+probdim*masternumnode, dsxigp_ );

    for (CI p=ximaps_[0].begin();p!=ximaps_[0].end();++p)
      dsxigp_[0][p->first] += 0.5*(1-eta[0])*(p->second);
    for (CI p=ximaps_[1].begin();p!=ximaps_[1].end();++p)
      dsxigp_[0][p->first] += 0.5*(1+eta[0])*(p->second);

    // evaluate the GP master coordinate derivatives
    INTEGRATOR::ResetPairedVector<1>( linsize+probdim*masternumnode, dmxigp_ );

    this->wrapper_.DerivXiGP2D(sele,mele,sxi[0],mxi[0],dsxigp_[0],dmxigp_[0],linsize);

    // evaluate the Jacobian derivative
    INTEGRATOR::ResetPairedVector( slavenumnode*probdim, derivjac_ );
    sele.DerivJacobian(sxi,derivjac_);


    //**********************************************************************
    // auxiliary variables
    //**********************************************************************
    // normalized normal at gp
    double gpn[3]      = {0.0,0.0,0.0};
    // deriv of x and y comp. of gpn (unit)
    INTEGRATOR::ResetPairedVector<2>( linsize+probdim*masternumnode, dnmap_unit_ );

    const double jac = dxdsxi*dsxideta;
    // calculate the averaged normal + derivative at gp level
    GP_Normal_DerivNormal(sele,mele,sval_,sderiv_,dsxigp_,&gpn[0],dnmap_unit_,linsize);
    // integrate scaling factor kappa
    GP_kappa(sele,lmval_,wgt,jac);
    // integrate the inner integral for later usage (for all found slave nodes)
    GP_VarWGap(sele,mele,sval_,mval_,lmval_,&gpn[0],wgt,jac);

    //**********************************************************************
    // compute segment LINEARIZATION
    //**********************************************************************
    for ( unsigned iter=0; iter<slavenumnode; ++iter )
    {
     GP_2D_VarWGap_Lin(iter,sele,mele,sval_,mval_,lmval_,&gpn[0],sderiv_,
         mderiv_,lmderiv_,dsxideta,dxdsxi,dxdsxidsxi,wgt,dsxigp_[0],dmxigp_[0],
         derivjac_,dnmap_unit_,ximaps_);

     GP_2D_kappa_Lin(iter,sele,lmval_,lmderiv_,dsxideta,dxdsxi,dxdsxidsxi,
         wgt,dsxigp_[0],derivjac_,ximaps_);
    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::IntegrateDerivCell3DAuxPlane(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    MORTAR::IntCell& cell,
    double* auxn )
{
  // number of nodes (slave, master)
  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // get slave and master nodal coords for Jacobian / GP evaluation
  sele.GetNodalCoords(scoord_);
  mele.GetNodalCoords(mcoord_);

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  int linsize = 0;
  for ( unsigned i=0; i<slavenumnode; ++i )
  {
    CoNode* cnode = static_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<this->wrapper_.nGP();++gp)
  {
    // coordinates and weight
    const double eta[2] = {this->wrapper_.Coordinate(gp,0), this->wrapper_.Coordinate(gp,1)};
    const double wgt = this->wrapper_.Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell.LocalToGlobal(eta,globgp,0);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // project Gauss point onto slave element
    // project Gauss point onto master element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(globgp,auxn,sele,sxi,sprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(globgp,auxn,mele,mxi,mprojalpha);

    // check GP projection (SLAVE)
    const double tol = 0.01;
    switch ( slavetype )
    {
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
        {
          std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
          std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
          std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
          std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
        }
        break;
      }
      default:
      {
        if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
        {
          std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
          std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
          std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
          std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
        }
        break;
      }
    }

    // check GP projection (MASTER)
    switch ( mastertype )
    {
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
        {
          std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
          std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
          std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
          std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
        }
        break;
      }
      default:
      {
        if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
        {
          std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
          std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
          std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
          std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
        }
        break;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(this->ShapeFcn(),sxi,lmval_,lmderiv_,slavenumnode,true);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval_,sderiv_,slavenumnode,false);
    mele.EvaluateShape(mxi,mval_,mderiv_,masternumnode,false);

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi,ssecderiv_,slavenumnode);

    // evaluate linearizations *******************************************
    // evaluate the intcell Jacobian derivative
    INTEGRATOR::ResetPairedVector( (slavenumnode+masternumnode)*probdim, derivjac_ );
    cell.DerivJacobian(derivjac_);

    // evaluate global GP coordinate derivative
    static LINALG::Matrix<3,1> svalcell;
    static LINALG::Matrix<3,2> sderivcell;
    cell.EvaluateShape(eta,svalcell,sderivcell);

    INTEGRATOR::ResetPairedVector( (slavenumnode+masternumnode)*probdim, lingp_ );

    for ( int v=0; v<3; ++v )
      for ( int d=0; d<3; ++d )
        for ( CI p=(cell.GetDerivVertex(v))[d].begin(); p!=(cell.GetDerivVertex(v))[d].end(); ++p )
          lingp_[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    INTEGRATOR::ResetPairedVector<2>( (slavenumnode+masternumnode)*probdim, dsxigp_ );

    this->wrapper_.DerivXiGP3DAuxPlane(sele,sxi,cell.Auxn(),dsxigp_,sprojalpha,cell.GetDerivAuxn(),lingp_);

    // evalute the GP master coordinate derivatives
    INTEGRATOR::ResetPairedVector<2>( (slavenumnode+masternumnode)*probdim, dmxigp_ );

    this->wrapper_.DerivXiGP3DAuxPlane(mele,mxi,cell.Auxn(),dmxigp_,mprojalpha,cell.GetDerivAuxn(),lingp_);

    // evaluate the integration cell Jacobian
    const double jac = cell.Jacobian();

    //**********************************************************************
    // frequently reused quantities
    //**********************************************************************
    double gpn[3] = {0.0,0.0,0.0};
    // deriv of x,y and z comp. of gpn (unit)
    INTEGRATOR::ResetPairedVector<3>( linsize, dnmap_unit_ );

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    // calculate the averaged normal + derivative at gp level
    GP_Normal_DerivNormal(sele,mele,sval_,sderiv_,dsxigp_,gpn,dnmap_unit_,linsize);
    // integrate scaling factor kappa
    GP_kappa(sele,lmval_,wgt,jac);
    // integrate the inner integral for later usage (for all found slave nodes)
    GP_VarWGap(sele,mele,sval_,mval_,lmval_,gpn,wgt,jac);

    //********************************************************************
    // compute cell linearization
    //********************************************************************
    for ( unsigned iter=0; iter<slavenumnode; ++iter )
    {
      GP_3D_kappa_Lin(iter,sele,lmval_,lmderiv_,wgt,jac,dsxigp_,derivjac_);

      GP_3D_VarWGap_Lin(iter,sele,mele,sval_,mval_,lmval_,gpn,sderiv_,mderiv_,
          lmderiv_,wgt,jac,dsxigp_,dmxigp_,derivjac_,dnmap_unit_);
    }

  }
  //**********************************************************************
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::IntegrateDerivSlaveElement( MORTAR::MortarElement& sele )
{

  for ( int gp=0; gp<this->wrapper_.nGP(); ++gp )
  {
    const double eta[2] = {this->wrapper_.Coordinate(gp,0), this->wrapper_.Coordinate(gp,1)};
    const double wgt    = this->wrapper_.Weight(gp);
    double sxi[2] = {0.0, 0.0};

    // get Gauss point in slave element coordinates
    sxi[0] = eta[0];
    sxi[1] = eta[1];

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(this->ShapeFcn(),sxi,lmval_,lmderiv_,slavenumnode,true);
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval_,sderiv_,slavenumnode,false);

    // integrate the slave jacobian
    const double jac = sele.Jacobian(sxi);

    // evaluate the slave Jacobian derivative
    INTEGRATOR::ResetPairedVector( slavenumnode*probdim, derivjac_ );
    sele.DerivJacobian(sxi,derivjac_);

    // *** SLAVE NODES ****************************************************
    for ( unsigned it=0; it<slavenumnode; ++it )
    {
      GP_AugA(it,sele,lmval_,wgt,jac);
      /*-----------------------------------------------------------------------*
       |   compute LINEARIZATION                                               |
       *-----------------------------------------------------------------------*/
      GP_AugA_Lin(it,sele,lmval_,wgt,jac,derivjac_);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::IntegrateDerivEle2D(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo )
{
  // *********************************************************************
  // Define slave quantities
  // *********************************************************************

  //consider entire slave element --> parameter space [-1,1]
  double sxia=-1.0;
  double sxib=1.0;

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  int linsize = 0;
  DRT::Node** mynodes = sele.Nodes();
  for ( unsigned i=0; i<slavenumnode; ++i )
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // get the gausspoints of this slave / master element pair
  const unsigned num_gps = projInfo.gaussPoints_.size();

  for ( unsigned i=0; i<num_gps; ++i )
  {
    int gp = projInfo.gaussPoints_[i];

    // coordinates and weight
    const double eta = this->wrapper_.Coordinate(gp,0);
    const double wgt = this->wrapper_.Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    const double sxi[2] = {eta, 0.0};

    // evaluate the two slave side Jacobians
    const double dxdsxi = sele.Jacobian(sxi);
//      double dsxideta = -0.5*sxia + 0.5*sxib; // dummy for gap

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult( this->ShapeFcn(), sxi, lmval_, lmderiv_,slavenumnode,true);

    // evaluate trace space shape functions
    sele.EvaluateShape(sxi,sval_,sderiv_,slavenumnode,false);

    const double* mxi = projInfo.uniqueMxi_[i].A();
//    int ncol      =   meles[nummaster]->NumNode();
//    LINALG::SerialDenseVector mval(ncol);
//    LINALG::SerialDenseMatrix mderiv(ncol,1);
//
//    // get master nodal coords for Jacobian / GP evaluation
//    LINALG::SerialDenseMatrix mcoord(3,meles[nummaster]->NumNode());
    mele.GetNodalCoords(mcoord_);

    // evaluate trace space shape functions
    mele.EvaluateShape(mxi,mval_,mderiv_,masternumnode,false);

    // get directional derivatives of sxia, sxib, mxia, mxib
    // --> derivatives of mxia/mxib not required
    INTEGRATOR::ResetPairedVector<4>( linsize+probdim*masternumnode, ximaps_ );
    bool startslave = true;
    bool endslave   = true;
    double mxia = -0.1;  //--> arbitrary value
    double mxib =  0.1;  //--> arbitrary value
    this->wrapper_.DerivXiAB2D(sele,sxia,sxib,mele,mxia,mxib,ximaps_,startslave,endslave,linsize);

    // evaluate the GP slave coordinate derivatives --> no entries
    INTEGRATOR::ResetPairedVector<1>( linsize+probdim*masternumnode, dsxigp_ );

    for (CI p=ximaps_[0].begin();p!=ximaps_[0].end();++p)
      dsxigp_[0][p->first] = 0.0;

    // evaluate the GP master coordinate derivatives
    INTEGRATOR::ResetPairedVector<1>( linsize+probdim*masternumnode, dmxigp_ );

    this->wrapper_.DerivXiGP2D(sele,mele,sxi[0],mxi[0],dsxigp_[0],dmxigp_[0],linsize);

    // evaluate the Jacobian derivative
    INTEGRATOR::ResetPairedVector( slavenumnode*probdim, derivjac_ );
    sele.DerivJacobian(sxi,derivjac_); //direct derivative if xi^1_g does not change

    //**********************************************************************
    // frequently reused quantities
    //**********************************************************************
    // normalized normal at gp
    double gpn[3]      = {0.0,0.0,0.0};
    // deriv of x and y comp. of gpn (unit)
    INTEGRATOR::ResetPairedVector<2>( linsize+probdim*masternumnode, dnmap_unit_ );

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************

    // calculate the averaged normal + derivative at gp level
    GP_Normal_DerivNormal(sele,mele,sval_,sderiv_,dsxigp_,&gpn[0],dnmap_unit_,linsize);
    // integrate scaling factor kappa
    GP_kappa(sele,lmval_,wgt,dxdsxi);
    // integrate the inner integral for later usage (for all found slave nodes)
    GP_VarWGap(sele,mele,sval_,mval_,lmval_,&gpn[0],wgt,dxdsxi);

    //**********************************************************************
    // compute LINEARIZATION
    //**********************************************************************
    switch ( cparams_.GetActionType() )
    {
      case MORTAR::eval_force_stiff:
      {
        for ( unsigned iter=0; iter<slavenumnode; ++iter )
        {
          GP_2D_VarWGap_Ele_Lin(iter,sele,mele,sval_,mval_,lmval_,&gpn[0],sderiv_,
              mderiv_,lmderiv_,dxdsxi,wgt,dmxigp_[0],derivjac_,dnmap_unit_);

          GP_2D_kappa_Ele_Lin(iter,sele,lmval_,lmderiv_,dxdsxi,wgt,derivjac_);
        }
        break;
      }
      default:
        // do nothing
        break;
    }
  } // End Loop over all GP
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::IntegrateDerivEle3D(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    bool boundary_ele,
    const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo )
{

  // get slave and master nodal coords for Jacobian / GP evaluation
  sele.GetNodalCoords(scoord_);


  int linsize = 0;
  DRT::Node** mynodes = sele.Nodes();
  for (unsigned i=0;i<slavenumnode;++i)
  {
    CoNode* cnode = static_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // Start integration if fast integration should be used or if there is no boundary element
  // for the fast_BS integration
  if (boundary_ele==false || this->wrapper_.IntegrationType()==INPAR::MORTAR::inttype_elements)
  {
    // get the gausspoints of this slave / master element pair
    const unsigned num_gps = projInfo.gaussPoints_.size();

    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************

    for ( unsigned i=0; i<num_gps; ++i )
    {
      const int gp = projInfo.gaussPoints_[i];

      // coordinates and weight
      const double eta[2] = {this->wrapper_.Coordinate(gp,0), this->wrapper_.Coordinate(gp,1)};
      const double wgt = this->wrapper_.Weight(gp);

      // get Gauss point in slave element coordinates
      const double sxi[2] = { eta[0], eta[1] };

      // evaluate Lagrange multiplier shape functions (on slave element)
      sele.EvaluateShapeLagMult(this->ShapeFcn(),sxi,lmval_,lmderiv_,slavenumnode,true);

      // evaluate trace space shape functions (on both elements)
      sele.EvaluateShape(sxi,sval_,sderiv_,slavenumnode,false);

      // evaluate the two Jacobians (int. cell and slave element)
      double jacslave = sele.Jacobian(sxi);

      // evaluate linearizations *******************************************
      // evaluate the slave Jacobian derivative
      INTEGRATOR::ResetPairedVector( slavenumnode*probdim, derivjac_ );
      sele.DerivJacobian(sxi,derivjac_);

      const double uniqueProjalpha  = projInfo.uniqueProjAlpha_[i];
      const double* uniqueMxi = projInfo.uniqueMxi_[i].A();


      mele.GetNodalCoords(mcoord_);

      // get mval
      mele.EvaluateShape(uniqueMxi,mval_,mderiv_,masternumnode,false);

      // evaluate the GP slave coordinate derivatives
      INTEGRATOR::ResetPairedVector<2>( 0, dsxigp_ );
      INTEGRATOR::ResetPairedVector<2>( linsize+masternumnode*probdim, dmxigp_ );

      this->wrapper_.DerivXiGP3D(sele,mele,sxi,uniqueMxi,dsxigp_,dmxigp_,uniqueProjalpha);

      //**********************************************************************
      // frequently reused quantities
      //**********************************************************************
      double gpn[3]      = {0.0,0.0,0.0};
      // deriv of x,y and z comp. of gpn (unit)
      INTEGRATOR::ResetPairedVector<3>( linsize+probdim*masternumnode, dnmap_unit_ );

      //**********************************************************************
      // evaluate at GP and lin char. quantities
      //**********************************************************************
      // calculate the averaged normal + derivative at gp level
      GP_Normal_DerivNormal(sele,mele,sval_,sderiv_,dsxigp_,&gpn[0],dnmap_unit_,linsize);
      // integrate scaling factor kappa
      GP_kappa(sele,lmval_,wgt,jacslave);
      // integrate the inner integral for later usage (for all found slave nodes)
      GP_VarWGap(sele,mele,sval_,mval_,lmval_,&gpn[0],wgt,jacslave);

      switch ( cparams_.GetActionType() )
      {
        case MORTAR::eval_force_stiff:
        {
          //********************************************************************
          // compute ele linearization
          //********************************************************************
          // loop over all slave nodes
          for (unsigned iter=0; iter<slavenumnode; ++iter)
          {
            GP_3D_kappa_Lin(iter,sele,lmval_,lmderiv_,wgt,jacslave,dsxigp_,derivjac_);

            GP_3D_VarWGap_Lin(iter,sele,mele,sval_,mval_,lmval_,&gpn[0],sderiv_,mderiv_,
                lmderiv_,wgt,jacslave,dsxigp_,dmxigp_,derivjac_,dnmap_unit_);
          }
          break;
        }
        default:
          // do nothing
          break;
      }
    }//GP-loop
  }
}


template CONTACT::AUG::IntegratorGeneric *
CONTACT::AUG::IntegratorGeneric::Create2D<DRT::Element::line2>(
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper );

template CONTACT::AUG::IntegratorGeneric *
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::line2>(
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper );
template CONTACT::AUG::IntegratorGeneric *
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::quad4>(
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper );
template CONTACT::AUG::IntegratorGeneric *
CONTACT::AUG::IntegratorGeneric::Create3D<DRT::Element::tri3>(
    DRT::Element::DiscretizationType mastertype,
    const CONTACT::ParamsInterface & cparams,
    CONTACT::CoIntegrator* wrapper );

#include "contact_augmented_integrator_list.H"
