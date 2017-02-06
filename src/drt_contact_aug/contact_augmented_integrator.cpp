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
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_serialdensevector.H"
#include <Epetra_Map.h>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AugmentedIntegrator::AugmentedIntegrator(Teuchos::ParameterList& params,
    DRT::Element::DiscretizationType eletype,
    const Epetra_Comm& comm)
    : CONTACT::CoIntegrator::CoIntegrator(params,eletype,comm)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivCell3DAuxPlane(
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
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  if (ShapeFcn() == INPAR::MORTAR::shape_dual || ShapeFcn() == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane supports no Dual shape functions for the "
        "augmented Lagrange solving strategy!");

  if (DRT::INPUT::IntegralValue<int>(imortar_,"LM_NODAL_SCALE"))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane supports currently no LM nodal scaling for the "
        "augmented Lagrange solving strategy!");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  sele.GetNodalCoords(scoord);
  mele.GetNodalCoords(mcoord);

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  int linsize = 0;
  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp,0);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // project Gauss point onto slave element
    // project Gauss point onto master element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector::Impl(sele)->ProjectGaussPointAuxn3D(globgp,auxn,sele,sxi,sprojalpha);
    MORTAR::MortarProjector::Impl(mele)->ProjectGaussPointAuxn3D(globgp,auxn,mele,mxi,mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt==DRT::Element::quad4 || sdt==DRT::Element::quad8 || sdt==DRT::Element::quad9)
    {
      if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);

    // evaluate linearizations *******************************************
    // evaluate the intcell Jacobian derivative
    GEN::pairedvector<int,double> jacintcellmap((nrow+ncol)*ndof);
    cell->DerivJacobian(eta,jacintcellmap);

    // evaluate global GP coordinate derivative
    static LINALG::Matrix<3,1> svalcell;
    static LINALG::Matrix<3,2> sderivcell;
    cell->EvaluateShape(eta,svalcell,sderivcell);

    GEN::pairedvector<int,LINALG::Matrix<3,1> > lingp((nrow+ncol)*ndof);

    for (int v=0;v<3;++v)
      for (int d=0; d<3; ++d)
        for (CI p=(cell->GetDerivVertex(v))[d].begin();p!=(cell->GetDerivVertex(v))[d].end();++p)
          lingp[p->first](d) += svalcell(v) * (p->second);

    // evalute the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int,double> > dsxigp(2,(nrow+ncol)*ndof);
    DerivXiGP3DAuxPlane(sele,sxi,cell->Auxn(),dsxigp,sprojalpha,cell->GetDerivAuxn(),lingp);

    // evalute the GP master coordinate derivatives
    std::vector<GEN::pairedvector<int,double> > dmxigp(2,(nrow+ncol)*ndof);
    DerivXiGP3DAuxPlane(mele,mxi,cell->Auxn(),dmxigp,mprojalpha,cell->GetDerivAuxn(),lingp);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    //**********************************************************************
    // frequently reused quantities
    //**********************************************************************
    double gpn[3]      = {0.0,0.0,0.0};
    // deriv of x,y and z comp. of gpn (unit)
    std::vector<GEN::pairedvector<int,double> > dnmap_unit(3,linsize);

    //**********************************************************************
    // evaluate at GP and lin char. quantities
    //**********************************************************************
    // calculate the averaged normal + derivative at gp level
    GP_Normal_DerivNormal(sele,mele,sval,sderiv,dsxigp,gpn,dnmap_unit,linsize);
    // integrate scaling factor kappa
    GP_kappa(sele,lmval,wgt,jac);
    // integrate the inner integral for later usage (for all found slave nodes)
    GP_VarWGap(sele,mele,sval,mval,lmval,gpn,wgt,jac);

    //********************************************************************
    // compute cell linearization
    //********************************************************************
    for (int iter=0;iter<nrow;++iter)
    {
      GP_3D_kappa_Lin(iter,sele,lmval,lmderiv,wgt,jac,dsxigp,jacintcellmap);

      GP_3D_VarWGap_Lin(iter,sele,mele,sval,mval,lmval,gpn,sderiv,mderiv,
          lmderiv,wgt,jac,dsxigp,dmxigp,jacintcellmap,dnmap_unit);
    }

  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivEle3D(
     MORTAR::MortarElement& sele,
     std::vector<MORTAR::MortarElement*> meles,
     bool *boundary_ele,
     bool *proj_,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (SolType()!=INPAR::CONTACT::solution_augmented)
    dserror("AugmentedIntegrator::IntegrateDerivEle3D: Method call for wrong SolType.");

  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

  // check input data
  for (int test=0;test<(int)meles.size();++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
      dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  }

  // contact with wear
  bool wear = false;
  if(wearlaw_!= INPAR::WEAR::wear_none)
    wear = true;

  int msize   = meles.size();
  int nrow    = sele.NumNode();
  int ndof    = dynamic_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);
  LINALG::SerialDenseVector svalmod(nrow);
  LINALG::SerialDenseMatrix sderivmod(nrow,2,true);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  //********************************************************************
  //  Boundary_segmentation test -- HasProj() check
  //  if a slave-node has no projection onto each master element
  //  --> Boundary_ele==true
  //********************************************************************
  INPAR::MORTAR::IntType integrationtype =
    DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_,"INTTYPE");

  //************************************************************************
  //Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  *boundary_ele=BoundarySegmCheck3D(sele,meles);

  int linsize = 0;
  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  // Start integration if fast integration should be used or if there is no boundary element
  // for the fast_BS integration
  if (*boundary_ele==false || integrationtype==INPAR::MORTAR::inttype_elements)
  {
    //**********************************************************************
    // loop over all Gauss points for integration
    //**********************************************************************
    for (int gp=0;gp<nGP();++gp)
    {
      int iter_proj=0;
      // coordinates and weight
      double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
      double wgt = Weight(gp);

      // note that the third component of sxi is necessary!
      // (although it will always be 0.0 of course)
      double sxi[2] = {0.0, 0.0};
      double mxi[2] = {0.0, 0.0};
      double projalpha = 0.0;

      // get Gauss point in slave element coordinates
      sxi[0] = eta[0];
      sxi[1] = eta[1];

      bool is_on_mele=true;

      // evaluate Lagrange multiplier shape functions (on slave element)
      sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

      // evaluate trace space shape functions (on both elements)
      sele.EvaluateShape(sxi,sval,sderiv,nrow);

      // evaluate the two Jacobians (int. cell and slave element)
      double jacslave = sele.Jacobian(sxi);

      // evaluate linearizations *******************************************
      // evaluate the slave Jacobian derivative
      GEN::pairedvector<int,double> jacslavemap(nrow*ndof);
      sele.DerivJacobian(sxi,jacslavemap);

      //**********************************************************************
      // loop over all mele
      //**********************************************************************
      bool isprojected = false;
      int uniqueMaEle = -1;
      double uniqueProjalpha = 0.0;
      double uniqueMxi[2] = {0.0,0.0};
      // --> find master elements with an unique projection
      for(int nummaster=0;nummaster<msize;++nummaster)
      {
        DRT::Element::DiscretizationType dt = meles[nummaster]->Shape();
        // project Gauss point onto master element
        MORTAR::MortarProjector::Impl(sele,*meles[nummaster])->ProjectGaussPoint3D(sele,sxi,*meles[nummaster],mxi,projalpha);

        is_on_mele=true;

        // check GP projection
        const double tol = 0.00;
        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
          {
            is_on_mele=false;
          }
        }
        else
        {
          if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
          {
            is_on_mele=false;
          }
        }

        // gp is valid and the current master element is the first feasible one
        if (is_on_mele==true and !isprojected)
        {
          isprojected = true;
          uniqueMaEle = nummaster;
          uniqueProjalpha = projalpha;
          uniqueMxi[0] = mxi[0];
          uniqueMxi[1] = mxi[1];
        }
        // found a second master element with a feasible projection
        else if (is_on_mele and isprojected)
        {
          if (projalpha<uniqueProjalpha)
          {
            uniqueMaEle = nummaster;
            uniqueProjalpha = projalpha;
            uniqueMxi[0] = mxi[0];
            uniqueMxi[1] = mxi[1];
          }
        }
      }//mele loop

      // use the found master element
      if (uniqueMaEle!=-1)
      {
        int nmnode  = meles[uniqueMaEle]->NumNode();
        LINALG::SerialDenseVector mval(nmnode);
        LINALG::SerialDenseMatrix mderiv(nmnode,2,true);
        LINALG::SerialDenseMatrix mcoord(3,meles[uniqueMaEle]->NumNode());
        meles[uniqueMaEle]->GetNodalCoords(mcoord);

        // get them in the case of tsi
        if (wear or gpslip_)
        {
          scoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
          mcoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,meles[uniqueMaEle]->NumNode()));
          lagmult   = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
          sele.GetNodalCoordsOld(*scoordold);
          meles[uniqueMaEle]->GetNodalCoordsOld(*mcoordold);
          sele.GetNodalLagMult(*lagmult);
        }

        // for both-sided wear
        LINALG::SerialDenseVector lm2val(nmnode);
        LINALG::SerialDenseMatrix lm2deriv(nmnode,2,true);

        // evaluate Lagrange multiplier shape functions (on slave element)
        if (WearSide() != INPAR::WEAR::wear_slave)
          meles[uniqueMaEle]->EvaluateShapeLagMult(ShapeFcn(),uniqueMxi,lm2val,lm2deriv,nmnode);

        *proj_=true;
        iter_proj+=1;

        // get mval
        meles[uniqueMaEle]->EvaluateShape(uniqueMxi,mval,mderiv,nmnode);

        // evaluate the GP slave coordinate derivatives
        std::vector<GEN::pairedvector<int,double> > dsxigp(2,0);
        std::vector<GEN::pairedvector<int,double> > dmxigp(2,linsize+nmnode*ndof);
        DerivXiGP3D(sele,*meles[uniqueMaEle],sxi,uniqueMxi,dsxigp,dmxigp,uniqueProjalpha);

        //**********************************************************************
        // frequently reused quantities
        //**********************************************************************
        double gpn[3]      = {0.0,0.0,0.0};
        // deriv of x,y and z comp. of gpn (unit)
        std::vector<GEN::pairedvector<int,double> > dnmap_unit(3,((nmnode*ndof)+linsize));
        //**********************************************************************
        // evaluate at GP and lin char. quantities
        //**********************************************************************
        // calculate the averaged normal + derivative at gp level
        GP_Normal_DerivNormal(sele,*meles[uniqueMaEle],sval,sderiv,dsxigp,&gpn[0],dnmap_unit,linsize);
        // integrate scaling factor kappa
        GP_kappa(sele,lmval,wgt,jacslave);
        // integrate the inner integral for later usage (for all found slave nodes)
        GP_VarWGap(sele,*meles[uniqueMaEle],sval,mval,lmval,&gpn[0],wgt,jacslave);

        if (cparams_ptr->GetActionType()==MORTAR::eval_force_stiff)
        {
          //********************************************************************
          // compute ele linearization
          //********************************************************************
          // loop over all slave nodes
          for (int iter=0; iter<nrow; ++iter)
          {
            GP_3D_kappa_Lin(iter,sele,lmval,lmderiv,wgt,jacslave,dsxigp,jacslavemap);

            GP_3D_VarWGap_Lin(iter,sele,*meles[uniqueMaEle],sval,mval,lmval,&gpn[0],sderiv,mderiv,
                lmderiv,wgt,jacslave,dsxigp,dmxigp,jacslavemap,dnmap_unit);
          }
        }
      }
      // warning, if an element which is declared not to be on the boundary by the above test
      // has non-projectable Gauss points
      else if (uniqueMaEle == -1 && *boundary_ele==false)
        std::cout << "*** warning *** Non-boundary element has non-projectable Gauss point \n" ;
    }//GP-loop
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivSlEle3D(
    MORTAR::MortarElement& sele,
    const Epetra_Comm& comm,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr);
  if (cparams_ptr.is_null())
    dserror("Cast to CONTACT::ParamsInterface failed!");
  IntegrateDerivSlEle3D(sele,comm,cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivSlEle3D(
     MORTAR::MortarElement& sele,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  if (Dim()!=3) dserror("ERROR: IntegrateDerivSlEle3D has been called for a 2-D problem!");

  const int nrow = sele.NumNode();
  const int ndof = Dim();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  for (int gp=0;gp<nGP();++gp)
  {
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt    = Weight(gp);
    double sxi[2] = {0.0, 0.0};

    // get Gauss point in slave element coordinates
    sxi[0] = eta[0];
    sxi[1] = eta[1];

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);

    // integrate the slave jacobian
    double jac = sele.Jacobian(sxi);

    // evaluate the slave Jacobian derivative
    GEN::pairedvector<int,double> derivjac(nrow*ndof);
    sele.DerivJacobian(sxi,derivjac);

    // *** SLAVE NODES ****************************************************
    for (int it=0;it<nrow;++it)
    {
      GP_AugA(it,sele,sval,lmval,wgt,jac);
      /*-----------------------------------------------------------------------*
       |   compute LINEARIZATION                                               |
       *-----------------------------------------------------------------------*/
      GP_AugA_Lin(it,sele,sval,lmval,wgt,jac,derivjac);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave/master seg.   (2D)hiermeier 04/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivSegment2D(
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
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

  // *********************************************************************
  // Prepare integration
  // *********************************************************************

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // Do the augmented Lagrange stuff...
  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: AugmentedIntegrator::IntegrateDerivSegment2D: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  bool linlm = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  if (bound) dserror("ERROR: AugmentedIntegrator::IntegrateDerivSegment2D: "
      "Boundary modification is not considered for the augmented Lagrange formulation!");

  if (sele.Shape() == DRT::Element::line3)
    dserror("ERROR: AugmentedIntegrator::IntegrateDerivSegment2D: "
          "Quadratic shape functions are not yet considered.");


  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,1);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  sele.GetNodalCoords(scoord);
  mele.GetNodalCoords(mcoord);

  // nodal lagrange multiplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  if (SolType() == INPAR::CONTACT::solution_augmented)
  {
    lagmult   = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalLagMult(*lagmult);
  }
  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  int linsize = 0;
  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);
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

  if (sxia!=-1.0 && mxib!=1.0)
    dserror("ERROR: First outer node is neither slave nor master node");
  if (sxib!=1.0 && mxia!=-1.0)
      dserror("ERROR: Second outer node is neither slave nor master node");

  if (sxia==-1.0) startslave = true;
  else            startslave = false;
  if (sxib==1.0)  endslave   = true;
  else            endslave   = false;

  // get directional derivatives of sxia, sxib, mxia, mxib
  std::vector<GEN::pairedvector<int,double> > ximaps(4,linsize+ndof*ncol);
  DerivXiAB2D(sele,sxia,sxib,mele,mxia,mxib,ximaps,startslave,endslave,linsize);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] = {0.0, 0.0};

    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    MORTAR::MortarProjector::Impl(sele,mele)->ProjectGaussPoint2D(sele,sxi,mele,mxi);

    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
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
      sele.EvaluateShapeLagMultLin(ShapeFcn(),sxi,lmval,lmderiv,nrow);
    else
    {
      sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    // evaluate linearizations *******************************************
    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);

    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    dynamic_cast<CONTACT::CoElement&>(sele).DJacDXi(djacdxi,sxi,ssecderiv);
    double dxdsxidsxi=djacdxi[0]; // only 2D here

    // evaluate the GP slave coordinate derivatives
    std::vector<GEN::pairedvector<int,double> > dsxigp(1,linsize+ndof*ncol);
    for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
      dsxigp[0][p->first] += 0.5*(1-eta[0])*(p->second);
    for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
      dsxigp[0][p->first] += 0.5*(1+eta[0])*(p->second);

    // evaluate the GP master coordinate derivatives
    GEN::pairedvector<int,double> dmxigp(linsize+ndof*ncol);
    DerivXiGP2D(sele,mele,sxi[0],mxi[0],dsxigp[0],dmxigp,linsize);

    // evaluate the Jacobian derivative
    GEN::pairedvector<int,double> derivjac(nrow*ndof);
    sele.DerivJacobian(sxi,derivjac);


    //**********************************************************************
    // auxiliary variables
    //**********************************************************************
    double gpn[3]      = {0.0,0.0,0.0};                                          // normalized normal at gp
    std::vector<GEN::pairedvector<int,double> > dnmap_unit(2,linsize+ndof*ncol); // deriv of x and y comp. of gpn (unit)

    double jac = dxdsxi*dsxideta;
    // calculate the averaged normal + derivative at gp level
    GP_Normal_DerivNormal(sele,mele,sval,sderiv,dsxigp,&gpn[0],dnmap_unit,linsize);
    // integrate scaling factor kappa
    GP_kappa(sele,lmval,wgt,jac);
    // integrate the inner integral for later usage (for all found slave nodes)
    GP_VarWGap(sele,mele,sval,mval,lmval,&gpn[0],wgt,jac);

    //**********************************************************************
    // compute segment LINEARIZATION
    //**********************************************************************
    for (int iter=0;iter<nrow;++iter)
    {
      GP_2D_VarWGap_Lin(iter,sele,mele,sval,mval,lmval,&gpn[0],sderiv,mderiv,lmderiv,dsxideta,
          dxdsxi,dxdsxidsxi,wgt,dsxigp[0],dmxigp,derivjac,dnmap_unit,ximaps);

      GP_2D_kappa_Lin(iter,sele,lmval,lmderiv,dsxideta,dxdsxi,dxdsxidsxi,wgt,dsxigp[0],
          derivjac,ximaps);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivEle2D(
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

  // *********************************************************************
  // Define slave quantities
  // *********************************************************************

  //consider entire slave element --> parameter space [-1,1]
  double sxia=-1.0;
  double sxib=1.0;

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ndof = Dim();

  // create empty vectors for shape fct. evaluation
  static LINALG::SerialDenseVector sval(nrow);
  static LINALG::SerialDenseMatrix sderiv(nrow,1);
  static LINALG::SerialDenseVector lmval(nrow);
  static LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // get slave nodal coords for Jacobian / GP evaluation
  static LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is only necessary for quadratic dual shape functions in 2D
  //bool duallin = false; // --> coming soon
  GEN::pairedvector<int,Epetra_SerialDenseMatrix> dualmap(ndof*nrow,0,Epetra_SerialDenseMatrix(nrow,nrow));
  if ((ShapeFcn() == INPAR::MORTAR::shape_dual || ShapeFcn() ==INPAR::MORTAR::shape_petrovgalerkin)
      && (sele.Shape()==MORTAR::MortarElement::line3 || sele.MoData().DerivDualShape()!=Teuchos::null))
  {
    //duallin=true;
    sele.DerivShapeDual(dualmap);
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // get numerical integration type
  INPAR::MORTAR::IntType inttype =
    DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_,"INTTYPE");

  //************************************************************************
  //Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  if(inttype==INPAR::MORTAR::inttype_elements_BS)
    *boundary_ele=BoundarySegmCheck2D(sele,meles);

  int linsize = 0;
  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (mynodes[i]);
    linsize += cnode->GetLinsize();
  }

  if (*boundary_ele==false || inttype==INPAR::MORTAR::inttype_elements)
  {
    //*************************************************************************
    //                loop over all Gauss points for integration
    //*************************************************************************
    for (int gp=0;gp<nGP();++gp)
    {
      bool kink_projection=false;

      // coordinates and weight
      double eta[2] = {Coordinate(gp,0), 0.0};
      double wgt = Weight(gp);

      // coordinate transformation sxi->eta (slave MortarElement->Overlap)
      double sxi[2] = {0.0, 0.0};
      sxi[0]= eta[0];

      // evaluate the two slave side Jacobians
      double dxdsxi = sele.Jacobian(sxi);
//      double dsxideta = -0.5*sxia + 0.5*sxib; // dummy for gap

      // evaluate Lagrange multiplier shape functions (on slave element)
      sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

      // evaluate trace space shape functions
      sele.EvaluateShape(sxi,sval,sderiv,nrow);

      //****************************************************************************************************************
      //                loop over all Master Elements
      //****************************************************************************************************************
      for (int nummaster=0;nummaster<(int)meles.size();++nummaster)
      {
        // project Gauss point onto master element
        double mxi[2] = {0.0, 0.0};
        MORTAR::MortarProjector::Impl(sele,*meles[nummaster])->ProjectGaussPoint2D(sele,sxi,*meles[nummaster],mxi);

        // gp on mele?
        if ((mxi[0]>=-1.0) && (mxi[0]<=1.0) && (kink_projection==false))
        {
          kink_projection=true;

          int ncol      =   meles[nummaster]->NumNode();
          LINALG::SerialDenseVector mval(ncol);
          LINALG::SerialDenseMatrix mderiv(ncol,1);

          // get master nodal coords for Jacobian / GP evaluation
          LINALG::SerialDenseMatrix mcoord(3,meles[nummaster]->NumNode());
          meles[nummaster]->GetNodalCoords(mcoord);

          // evaluate trace space shape functions
          meles[nummaster]->EvaluateShape(mxi,mval,mderiv,ncol);

          // get directional derivatives of sxia, sxib, mxia, mxib --> derivatives of mxia/mxib not required
          std::vector<GEN::pairedvector<int,double> > ximaps(4,linsize+ndof*ncol);
          bool startslave = true;
          bool endslave   = true;
          double mxia = -0.1;  //--> arbitrary value
          double mxib =  0.1;  //--> arbitrary value
          DerivXiAB2D(sele,sxia,sxib,*meles[nummaster],mxia,mxib,ximaps,startslave,endslave,linsize);

          // evaluate the GP slave coordinate derivatives --> no entries
          std::vector<GEN::pairedvector<int,double> >dsxigp(1,linsize+ndof*ncol);
          for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
            dsxigp[0][p->first] = 0.0;

          // evaluate the GP master coordinate derivatives
          GEN::pairedvector<int,double> dmxigp(linsize+ndof*ncol);
          DerivXiGP2D(sele,*meles[nummaster],sxi[0],mxi[0],dsxigp[0],dmxigp,linsize);

          // evaluate the Jacobian derivative
          GEN::pairedvector<int,double> derivjac(nrow*ndof);
          sele.DerivJacobian(sxi,derivjac); //direct derivative if xi^1_g does not change

          //**********************************************************************
          // frequently reused quantities
          //**********************************************************************
          double gpn[3]      = {0.0,0.0,0.0};  // normalized normal at gp
          std::vector<GEN::pairedvector<int,double> > dnmap_unit(2,(linsize+ndof*ncol)); // deriv of x and y comp. of gpn (unit)

          //**********************************************************************
          // evaluate at GP and lin char. quantities
          //**********************************************************************

          // calculate the averaged normal + derivative at gp level
          GP_Normal_DerivNormal(sele,*meles[nummaster],sval,sderiv,dsxigp,&gpn[0],dnmap_unit,linsize);
          // integrate scaling factor kappa
          GP_kappa(sele,lmval,wgt,dxdsxi);
          // integrate the inner integral for later usage (for all found slave nodes)
          GP_VarWGap(sele,*meles[nummaster],sval,mval,lmval,&gpn[0],wgt,dxdsxi);

          //**********************************************************************
          // compute LINEARIZATION
          //**********************************************************************
          for (int iter=0;iter<nrow;++iter)
          {
            GP_2D_VarWGap_Ele_Lin(iter,sele,*meles[nummaster],sval,mval,lmval,&gpn[0],sderiv,mderiv,lmderiv,
                dxdsxi,wgt,dmxigp,derivjac,dnmap_unit);

            GP_2D_kappa_Ele_Lin(iter,sele,lmval,lmderiv,dxdsxi,wgt,derivjac);
          }
        }
      }//End Loop over all Master Elements
    } // End Loop over all GP
  }//boundary_ele check

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivSlEle2D(
     MORTAR::MortarElement& sele,
     const Epetra_Comm& comm,
     const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcp_dynamic_cast<CONTACT::ParamsInterface>(mparams_ptr);
  if (cparams_ptr.is_null())
    dserror("Cast to CONTACT::ParamsInterface failed!");
  IntegrateDerivSlEle2D(sele,comm,cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedIntegrator::IntegrateDerivSlEle2D(
     MORTAR::MortarElement& sele,
     const Epetra_Comm& comm,
     const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  if (cparams_ptr.is_null())
    dserror("ERROR: The contact parameter interface pointer is undefined!");

  if (Dim()!=2) dserror("ERROR: IntegrateDerivSlEle2D has been called for a 3-D problem!");

  // Input was already checked in the previous call of IntegrateAvgWgapSegment2D,
  // so we can skip it here.

  // *********************************************************************
  // Prepare integration
  // *********************************************************************
  // number of nodes (slave, master)
  const int nrow = sele.NumNode();
  const int ndof = Dim();

//  // Split the nodes of the current slave element in an active and an inactive set.
//  std::vector<int> rowNumActive;
//  std::vector<int> rowNumInactive;
//
//  for (int i=0;i<nrow;++i)
//    if (augActiveSlaveNodes_->LID(mynodes[i]->Id())!=-1)
//      rowNumActive.push_back(i);
//
//  int nrowActive   = (int) rowNumActive.size();
//  int nrowInactive = (int) rowNumInactive.size();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double sxi[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);

    sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions
    sele.EvaluateShape(sxi,sval,sderiv,nrow);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);

    // evaluate linearizations *******************************************
    // evaluate the Jacobian derivative
    GEN::pairedvector<int,double> derivjac(nrow*ndof);
    sele.DerivJacobian(sxi,derivjac);

    // *** SLAVE NODES ****************************************************
    for (int it=0;it<nrow;++it)
    {
      GP_AugA(it,sele,sval,lmval,wgt,dxdsxi);
      /*-----------------------------------------------------------------------*
       |   compute LINEARIZATION                                               |
       *-----------------------------------------------------------------------*/
      GP_AugA_Lin(it,sele,sval,lmval,wgt,dxdsxi,derivjac);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_kappa(
     MORTAR::MortarElement& sele,
     LINALG::SerialDenseVector& lmval,
     double& wgt, double& jac)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!");

  // number of nodes (slave)
  int nrow = sele.NumNode();

  // add to node
  for (int j=0;j<nrow;++j)
  {
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[j]);

    double val = 0.0;
    val = lmval[j]*jac*wgt;

    // add current Gauss point's contribution kappaseg
    cnode->AddKappaValue(val);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_2D_kappa_Lin(
    int& iter,
    MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& lmderiv,
    double& dsxideta, double& dxdsxi,
    double& dxdsxidsxi,
    double& wgt,
    const GEN::pairedvector<int,double>& dsxigp,
    const GEN::pairedvector<int,double>& derivjac,
    const std::vector<GEN::pairedvector<int,double> >& ximaps)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[iter]);
  std::map<int,double>& kappaLinMap = cnode->CoData().GetKappaLin();

  double fac = 0.0;

  // (0) Lin(LmSlave) - slave GP coordinates
  fac  = lmderiv(iter,0)*dxdsxi;
  // (1) Lin(dxdsxi) - slave GP coordinates
  fac += lmval[iter]*dxdsxidsxi;
  fac *= wgt*dsxideta;
  for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
    kappaLinMap[p->first] += fac*(p->second);

  // (2) Lin(dsxideta) - segment end coordinates
  fac = wgt*lmval[iter]*dxdsxi;
  for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
    kappaLinMap[p->first] -= 0.5*fac*(p->second);
  for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
    kappaLinMap[p->first] += 0.5*fac*(p->second);

  // (3) Lin(dxdsxi) - slave GP Jacobian
  fac = wgt*lmval[iter]*dsxideta;
  for (CI p=derivjac.begin();p!=derivjac.end();++p)
    kappaLinMap[p->first] += fac*(p->second);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_2D_kappa_Ele_Lin(
    int& iter,
    MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& lmderiv,
    double& dxdsxi,double& wgt,
    const GEN::pairedvector<int,double>& derivjac)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[iter]);
  std::map<int,double>& kappaLinMap = cnode->CoData().GetKappaLin();

  // (0) Lin(LmSlave) - slave GP coordinates --> 0
  // (1) Lin(dxdsxi) - slave GP coordinates --> 0
  // (2) Lin(dsxideta) - segment end coordinates --> 0
  // (3) Lin(dxdsxi) - slave GP Jacobian
  double fac = wgt*lmval[iter];
  for (CI p=derivjac.begin();p!=derivjac.end();++p)
    kappaLinMap[p->first] += fac*(p->second);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_3D_kappa_Lin(
    int& iter,
    MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseMatrix& lmderiv,
    double& wgt, double& jac,
    const std::vector<GEN::pairedvector<int,double> >& dsxigp,
    const GEN::pairedvector<int,double>& jacintcellmap)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[iter]);
  std::map<int,double>& kappaLinMap = cnode->CoData().GetKappaLin();
  double fac = 0.0;

  // (1) Lin(Phi) - dual shape functions
  // this vanishes here since there are no deformation-dependent dual functions

  // (2) Lin(LmSlave) - slave GP coordinates
  for (int i=0;i<(int) dsxigp.size();++i)
  {
    fac = wgt*lmderiv(iter,i)*jac;
    for (CI p=dsxigp[i].begin();p!=dsxigp[i].end();++p)
      kappaLinMap[p->first] += fac*(p->second);
  }

  // (3) Lin(dsxideta) - intcell GP Jacobian
  fac = wgt*lmval[iter];
  for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
    kappaLinMap[p->first] += fac*(p->second);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_Normal_DerivNormal(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseMatrix& sderiv,
    const std::vector<GEN::pairedvector<int,double> >& dsxigp,
    double* gpn,
    std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
    int& linsize)
{
  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: GP_2D_Normal_DerivNormal: Null Pointer!");

  // number of slave nodes
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  for (int i=0;i<nrow;++i)
  {
    MORTAR::MortarNode* mymrtnode = dynamic_cast<MORTAR::MortarNode*> (snodes[i]);
    gpn[0]+=sval[i]*mymrtnode->MoData().n()[0];
    gpn[1]+=sval[i]*mymrtnode->MoData().n()[1];
    gpn[2]+=sval[i]*mymrtnode->MoData().n()[2];
  }

  // normalize interpolated GP normal back to length 1.0 !!!
  double lengthn = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
  if (lengthn<1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i=0;i<3;++i)
    gpn[i]/=lengthn;

  // ******************************
  // Linearization of the gp-normal
  // ******************************
  // *** 2-D case ***********************************************
  if (Dim()==2)
  {
    // build directional derivative of slave GP normal (non-unit)
    GEN::pairedvector<int,double> dmap_nxsl_gp(ncol*ndof+linsize);
    GEN::pairedvector<int,double> dmap_nysl_gp(ncol*ndof+linsize);

    for (int i=0;i<nrow;++i)
    {
      CoNode* snode = dynamic_cast<CoNode*> (snodes[i]);

      GEN::pairedvector<int,double>& dmap_nxsl_i = dynamic_cast<CoNode*>(snodes[i])->CoData().GetDerivN()[0];
      GEN::pairedvector<int,double>& dmap_nysl_i = dynamic_cast<CoNode*>(snodes[i])->CoData().GetDerivN()[1];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*snode->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*snode->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
      }
    }

    // build directional derivative of slave GP normal (unit)
    double ll = lengthn*lengthn;
    double sxsx = gpn[0]*gpn[0]*ll; // gpn is the unit normal --> multiplication with ll
    double sxsy = gpn[0]*gpn[1]*ll; // to get the non-unit normal
    double sysy = gpn[1]*gpn[1]*ll;

    for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    {
      dnmap_unit[0][p->first] += 1/lengthn*(p->second);
      dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsx*(p->second);
      dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
    }

    for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
    {
      dnmap_unit[1][p->first] += 1/lengthn*(p->second);
      dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sysy*(p->second);
      dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
    }
  }
  // *** 3-D case ***********************************************
  else if (Dim()==3)
  {
    // build directional derivative of slave GP normal (non-unit)
    GEN::pairedvector<int,double> dmap_nxsl_gp(linsize);
    GEN::pairedvector<int,double> dmap_nysl_gp(linsize);
    GEN::pairedvector<int,double> dmap_nzsl_gp(linsize);

    for (int i=0;i<nrow;++i)
    {
      CoNode* cnode = dynamic_cast<CoNode*> (snodes[i]);

      GEN::pairedvector<int,double>& dmap_nxsl_i = cnode->CoData().GetDerivN()[0];
      GEN::pairedvector<int,double>& dmap_nysl_i = cnode->CoData().GetDerivN()[1];
      GEN::pairedvector<int,double>& dmap_nzsl_i = cnode->CoData().GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*cnode->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*cnode->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*cnode->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*cnode->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*cnode->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*cnode->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }
    }

    double ll = lengthn*lengthn;
    double sxsx = gpn[0]*gpn[0]*ll;
    double sxsy = gpn[0]*gpn[1]*ll;
    double sxsz = gpn[0]*gpn[2]*ll;
    double sysy = gpn[1]*gpn[1]*ll;
    double sysz = gpn[1]*gpn[2]*ll;
    double szsz = gpn[2]*gpn[2]*ll;

    for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    {
      dnmap_unit[0][p->first] += 1/lengthn*(p->second);
      dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsx*(p->second);
      dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
      dnmap_unit[2][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsz*(p->second);
    }

    for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
    {
      dnmap_unit[1][p->first] += 1/lengthn*(p->second);
      dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sysy*(p->second);
      dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
      dnmap_unit[2][p->first] -= 1/(lengthn*lengthn*lengthn)*sysz*(p->second);
    }

    for (CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
    {
      dnmap_unit[2][p->first] += 1/lengthn*(p->second);
      dnmap_unit[2][p->first] -= 1/(lengthn*lengthn*lengthn)*szsz*(p->second);
      dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsz*(p->second);
      dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sysz*(p->second);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_VarWGap(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseVector& lmval,
    double* gpn,
    double& wgt,double& jac)
{
  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  if(!snodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  if(!mnodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();

  // loop over all possible active slave nodes
  for (int i=0;i<nrow;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*>(snodes[i]);

    // *** slave node contributions ***
    // loop over all slave nodes
    for (int k=0;k<nrow;++k)
    {
      CoNode* snode = dynamic_cast<CoNode*>(snodes[k]);
      for (int kdof=0;kdof<Dim();++kdof)
      {
        int sGid   = snode->Id();
        int sDofID = snode->Dofs()[kdof];
        double val = lmval[i] * sval[k] * gpn[kdof] * wgt * jac;
        cnode->AddVarWGapSl(sDofID,sGid,val);
      }
    }

    // *** master node contributions ***
    // loop over all master nodes
    for (int j=0;j<ncol;++j)
    {
      CoNode* mnode = dynamic_cast<CoNode*>(mnodes[j]);
      for (int jdof=0;jdof<Dim();++jdof)
      {
        int mGid   = mnode->Id();
        int mDofID = mnode->Dofs()[jdof];
        double val = lmval[i] * mval[j] * gpn[jdof] * wgt * jac;
        cnode->AddVarWGapMa(mDofID,mGid,val);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_2D_VarWGap_Lin(
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
    const std::vector<GEN::pairedvector<int,double> >& ximaps)
{
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CoNode* cnode = dynamic_cast<CoNode*>(snodes[iter]);
  if (!cnode) dserror("ERROR: GP_2D_VarWGap_Lin: Null pointer!");

  // *** integrate lin varWGapSl ****************************************
  for (int k=0;k<nrow;++k)
  {
    CoNode* snode = dynamic_cast<CoNode*>(snodes[k]);

    double val[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates
    val[0] = wgt*lmderiv(iter,0)*sval[k]*dsxideta*dxdsxi;
    // (1) Lin(n-Slave) - normal direction
    val[1] = wgt*lmval[iter]*sval[k]*dsxideta*dxdsxi;
    // (2) Lin(NSlave) - slave GP coordinates
    val[2] = wgt*lmval[iter]*sderiv(k,0)*dsxideta*dxdsxi;
    // (3) Lin(dsxideta) - segment end coordinates
    val[3] = wgt*lmval[iter]*sval[k]*dxdsxi;
    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[4] = wgt*lmval[iter]*sval[k]*dsxideta;
    // (5) Lin(dxdsxi) - slave GP coordinates
    val[5] = wgt*lmval[iter]*sval[k]*dsxideta*dxdsxidsxi;

    for (int kdof=0;kdof<snode->NumDof();++kdof)
    {
      int sDofId = snode->Dofs()[kdof];
      std::map<int,double>& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl()[sDofId];

      // (0,2,5) - slave GP coordinates
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        varWGapLinSlMap[p->first] += (val[0]+val[2]+val[5])*p->second*gpn[kdof];
      // (1) - normal direction
      for (CI p=dnmap_unit[kdof].begin();p!=dnmap_unit[kdof].end();++p)
        varWGapLinSlMap[p->first] += val[1]*p->second;
      // (3) - segment end coordinates
      for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
        varWGapLinSlMap[p->first] -= 0.5*val[3]*p->second*gpn[kdof];
      for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
        varWGapLinSlMap[p->first] += 0.5*val[3]*p->second*gpn[kdof];
      // (4) - Slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinSlMap[p->first] += val[4]*p->second*gpn[kdof];
    }
  }
  // *** integrate lin varWGapMa ****************************************
  for (int l=0;l<ncol;++l)
  {
    CoNode* mnode = dynamic_cast<CoNode*>(mnodes[l]);

    double val[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates
    val[0] = wgt*lmderiv(iter,0)*mval[l]*dsxideta*dxdsxi;
    // (1) Lin(n-Slave) - normal direction
    val[1] = wgt*lmval[iter]*mval[l]*dsxideta*dxdsxi;
    // (2) Lin(NMaster) - master GP coordinates
    val[2] = wgt*lmval[iter]*mderiv(l,0)*dsxideta*dxdsxi;
    // (3) Lin(dsxideta) - segment end coordinates
    val[3] = wgt*lmval[iter]*mval[l]*dxdsxi;
    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[4] = wgt*lmval[iter]*mval[l]*dsxideta;
    // (5) Lin(dxdsxi) - slave GP coordinates
    val[5] = wgt*lmval[iter]*mval[l]*dsxideta*dxdsxidsxi;

    for (int ldof=0;ldof<mnode->NumDof();++ldof)
    {
      int mDofId = mnode->Dofs()[ldof];
      std::map<int,double>& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa()[mDofId];

      // (0,5) - slave GP coordinates
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        varWGapLinMaMap[p->first] += (val[0]+val[5])*p->second*gpn[ldof];
      // (1) - normal direction
      for (CI p=dnmap_unit[ldof].begin();p!=dnmap_unit[ldof].end();++p)
        varWGapLinMaMap[p->first] += val[1]*p->second;
      // (2) - master GP coordinates
      for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
        varWGapLinMaMap[p->first] += val[2]*p->second*gpn[ldof];
      // (3) - segement end coordinates
      for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
        varWGapLinMaMap[p->first] -= 0.5*val[3]*p->second*gpn[ldof];
      for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
        varWGapLinMaMap[p->first] += 0.5*val[3]*p->second*gpn[ldof];
      // (4) - slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinMaMap[p->first] += val[4]*p->second*gpn[ldof];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_2D_VarWGap_Ele_Lin(
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
    const std::vector<GEN::pairedvector<int,double> >& dnmap_unit)
{
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CoNode* cnode = dynamic_cast<CoNode*>(snodes[iter]);
  if (!cnode) dserror("ERROR: GP_2D_VarWGap_Ele_Lin: Null pointer!");

  // *** integrate lin varWGapSl ****************************************
  for (int k=0;k<nrow;++k)
  {
    CoNode* snode = dynamic_cast<CoNode*>(snodes[k]);

    double val[2] = {0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates --> 0

    // (1) Lin(n-Slave) - normal direction
    val[0] = wgt*lmval[iter]*sval[k]*dxdsxi;
    // (2) Lin(NSlave) - slave GP coordinates --> 0

    // (3) Lin(dsxideta) - segment end coordinates --> 0

    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[1] = wgt*lmval[iter]*sval[k];
    // (5) Lin(dxdsxi) - slave GP coordinates --> 0

    for (int kdof=0;kdof<snode->NumDof();++kdof)
    {
      int sDofId = snode->Dofs()[kdof];
      std::map<int,double>& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl()[sDofId];

      // (1) - normal direction
      for (CI p=dnmap_unit[kdof].begin();p!=dnmap_unit[kdof].end();++p)
        varWGapLinSlMap[p->first] += val[0]*p->second;
      // (4) - Slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinSlMap[p->first] += val[1]*p->second*gpn[kdof];
    }
  }
  // *** integrate lin varWGapMa ****************************************
  for (int l=0;l<ncol;++l)
  {
    CoNode* mnode = dynamic_cast<CoNode*>(mnodes[l]);

    double val[3] = {0.0, 0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates --> 0

    // (1) Lin(n-Slave) - normal direction
    val[0] = wgt*lmval[iter]*mval[l]*dxdsxi;
    // (2) Lin(NMaster) - master GP coordinates
    val[1] = wgt*lmval[iter]*mderiv(l,0)*dxdsxi;
    // (3) Lin(dsxideta) - segment end coordinates --> 0

    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[2] = wgt*lmval[iter]*mval[l];
    // (5) Lin(dxdsxi) - slave GP coordinates --> 0

    for (int ldof=0;ldof<mnode->NumDof();++ldof)
    {
      int mDofId = mnode->Dofs()[ldof];
      std::map<int,double>& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa()[mDofId];

      // (1) - normal direction
      for (CI p=dnmap_unit[ldof].begin();p!=dnmap_unit[ldof].end();++p)
        varWGapLinMaMap[p->first] += val[0]*p->second;
      // (2) - master GP coordinates
      for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
        varWGapLinMaMap[p->first] += val[1]*p->second*gpn[ldof];
      // (4) - slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinMaMap[p->first] += val[2]*p->second*gpn[ldof];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_3D_VarWGap_Lin(
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
    const std::vector<GEN::pairedvector<int,double> >& dnmap_unit)
{
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CoNode* cnode = dynamic_cast<CoNode*>(snodes[iter]);
  if (!cnode) dserror("ERROR: GP_2D_VarWGap_Lin: Null pointer!");

  // *** integrate lin varWGapSl ****************************************
  for (int k=0;k<nrow;++k)
  {
    CoNode* snode = dynamic_cast<CoNode*>(snodes[k]);
    double val[4] = {0.0,0.0,0.0,0.0};

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions
    // (2) Lin(LMShape) & Lin(NSlave) - 1st slave GP coordinate
    val[0]  = wgt*lmderiv(iter,0)*sval[k]*jac;
    val[0] += wgt*lmval[iter]*sderiv(k,0)*jac;
    // (3) Lin(LMShape) & Lin(NSlave) - 2nd slave GP coordinate
    val[1]  = wgt*lmderiv(iter,1)*sval[k]*jac;
    val[1] += wgt*lmval[iter]*sderiv(k,1)*jac;
    // (4) Lin(dsxideta) - intcell GP Jacobian
    val[2]  = wgt*lmval[iter]*sval[k];
    // (5) Lin(n-Slave) - normal direction
    val[3]  = wgt*lmval[iter]*sval[k]*jac;

    for (int kdof=0;kdof<snode->NumDof();++kdof)
    {
      int sDofId = snode->Dofs()[kdof];
      std::map<int,double>& varWGapLinSlMap = cnode->CoData().GetVarWGapLinSl()[sDofId];

      for (int i=0; i<(int) dsxigp.size();++i)
        for (CI p=dsxigp[i].begin(); p!=dsxigp[i].end(); ++p)
          varWGapLinSlMap[p->first] += val[i]*(p->second)*gpn[kdof];

      for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
        varWGapLinSlMap[p->first] += val[2]*(p->second)*gpn[kdof];

      for (CI p=dnmap_unit[kdof].begin();p!=dnmap_unit[kdof].end();++p)
        varWGapLinSlMap[p->first] += val[3]*(p->second);
    }
  }
  // *** integrate lin varWGapMa ****************************************
  for (int l=0;l<ncol;++l)
  {
    CoNode* mnode = dynamic_cast<CoNode*>(mnodes[l]);

    double val[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions
    // (2) Lin(LMShape) - 1st slave GP coordinate
    val[0] = wgt*lmderiv(iter,0)*mval[l]*jac;
    // (3) Lin(LMShape) - 2nd slave GP coordinate
    val[1] = wgt*lmderiv(iter,1)*mval[l]*jac;
    // (4) Lin(NMaster) - 1st master GP coordinate
    val[2] = wgt*lmval[iter]*mderiv(l,0)*jac;
    // (5) Lin(NMaster) - 1st master GP coordinate
    val[3] += wgt*lmval[iter]*mderiv(l,1)*jac;
    // (6) Lin(dsxideta) - intcell GP Jacobian
    val[4]  = wgt*lmval[iter]*mval[l];
    // (7) Lin(n-Slave) - normal direction
    val[5]  = wgt*lmval[iter]*mval[l]*jac;

    for (int ldof=0;ldof<mnode->NumDof();++ldof)
    {
      int mDofId = mnode->Dofs()[ldof];
      std::map<int,double>& varWGapLinMaMap = cnode->CoData().GetVarWGapLinMa()[mDofId];

      for (int i=0;i<(int) dsxigp.size();++i)
        for (CI p=dsxigp[i].begin(); p!=dsxigp[i].end(); ++p)
          varWGapLinMaMap[p->first] += val[i]*(p->second)*gpn[ldof];

      for (int i=0;i<(int) dmxigp.size();++i)
        for (CI p=dmxigp[i].begin(); p!=dmxigp[i].end(); ++p)
          varWGapLinMaMap[p->first] += val[i+2]*(p->second)*gpn[ldof];

      for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
        varWGapLinMaMap[p->first] += val[4]*(p->second)*gpn[ldof];

      for (CI p=dnmap_unit[ldof].begin();p!=dnmap_unit[ldof].end();++p)
        varWGapLinMaMap[p->first] += val[5]*(p->second);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_AugA(
    int& it,
    MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& lmval,
    double& wgt,
    double& jac)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!");

  CoNode* cnode = dynamic_cast<CoNode*>(snodes[it]);
  double& augA = cnode->CoData().GetAugA();

  augA += lmval[it]*jac*wgt;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void inline CONTACT::AugmentedIntegrator::GP_AugA_Lin(
    int& it,
    MORTAR::MortarElement& sele,
    LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& lmval,
    double& wgt,
    double& jac,
    const GEN::pairedvector<int,double>& derivjac)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!");

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(snodes[it]);
  std::map<int,double>& augALinMap = cnode->CoData().GetAugALin();

  // Lin(dxdsxi) - slave GP Jacobian
  double val = wgt*lmval[it];

  // (2) - slave GP Jacobian
  for (CI p=derivjac.begin();p!=derivjac.end();++p)
    augALinMap[p->first] += val*p->second;


  return;
}
