/*!----------------------------------------------------------------------
\file contact_coupling2d.cpp

\maintainer Philipp Farah, Alexander Seitz

\brief Classes for mortar contact coupling in 2D.

*-----------------------------------------------------------------------*/

#include "contact_coupling2d.H"
#include "contact_integrator.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "contact_node.H"
#include "contact_interpolator.H"

#include "../drt_contact_aug/contact_augmented_integrator.H"

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_projector.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
CONTACT::CoCoupling2d::CoCoupling2d(DRT::Discretization& idiscret,
                                    int dim, bool quad,
                                    Teuchos::ParameterList& params,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::MortarElement& mele) :
MORTAR::Coupling2d(idiscret,dim,quad,params,sele,mele),
stype_(DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params,"STRATEGY"))
{
  // empty constructor

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate slave / master overlap (public)                 popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoCoupling2d::IntegrateOverlap()
{
  // explicitly defined shape function type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror(
        "ERROR: IntegrateOverlap called without specific shape function defined!");

  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj_ entries integrate the Mortar  */
  /* matrices D and M and the weighted gap function g~ on the overlap   */
  /* of the current sl / ma pair.                                       */
  /**********************************************************************/

  // no integration if no overlap
  if (!overlap_)
    return false;

  // set segmentation status of all slave nodes
  // (hassegment_ of a slave node is true if ANY segment/cell
  // is integrated that contributes to this slave node)
  int nnodes = SlaveElement().NumNode();
  DRT::Node** mynodes = SlaveElement().Nodes();
  if (!mynodes)
    dserror("ERROR: Null pointer!");
  for (int k = 0; k < nnodes; ++k)
  {
    MORTAR::MortarNode* mycnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mycnode)
      dserror("ERROR: Null pointer!");
    mycnode->HasSegment() = true;
  }

  //local working copies of input variables
  double sxia = xiproj_[0];
  double sxib = xiproj_[1];
  double mxia = xiproj_[2];
  double mxib = xiproj_[3];

  // create a CONTACT integrator instance with correct NumGP and Dim
  Teuchos::RCP<CONTACT::CoIntegrator> integrator = Teuchos::null;
  if (stype_==INPAR::CONTACT::solution_augmented)
    integrator = Teuchos::rcp(new CONTACT::AugmentedIntegrator(imortar_,SlaveElement().Shape(),Comm(),Teuchos::null));
  else
    integrator=Teuchos::rcp(new CONTACT::CoIntegrator(imortar_,SlaveElement().Shape(),Comm()));

  // *******************************************************************
  // different options for mortar integration
  // *******************************************************************
  // (1) no quadratic element(s) involved -> linear LM interpolation
  // (2) quadratic element(s) involved -> quadratic LM interpolation
  // (3) quadratic element(s) involved -> linear LM interpolation
  // (4) quadratic element(s) involved -> piecew. linear LM interpolation
  // *******************************************************************
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  // *******************************************************************
  // cases (1), (2) and (3)
  // *******************************************************************
  if (!Quad() || (Quad() && lmtype == INPAR::MORTAR::lagmult_quad)
      || (Quad() && lmtype == INPAR::MORTAR::lagmult_lin))
  {
    // ***********************************************************
    //                   Integrate stuff !!!                    //
    // ***********************************************************
    if (stype_==INPAR::CONTACT::solution_augmented)
      Teuchos::rcp_dynamic_cast<CONTACT::AugmentedIntegrator>(integrator)->IntegrateDerivSegment2D(SlaveElement(),sxia,sxib,MasterElement(),mxia,mxib,Comm());
    else
      integrator->IntegrateDerivSegment2D(SlaveElement(),sxia,sxib,MasterElement(),mxia,mxib,Comm());
    // ***********************************************************
    //                   END INTEGRATION !!!                    //
    // ***********************************************************
  }

  // *******************************************************************
  // case (4)
  // *******************************************************************
  else if (Quad() && lmtype == INPAR::MORTAR::lagmult_pwlin)
  {
    dserror("ERROR: Piecewise linear LM not (yet?) implemented in 2D");
  }

  // *******************************************************************
  // undefined case
  // *******************************************************************
  else if (Quad() && lmtype == INPAR::MORTAR::lagmult_undefined)
  {
    dserror("Lagrange multiplier interpolation for quadratic elements undefined\n"
            "If you are using 2nd order mortar elements, you need to specify LM_QUAD in MORTAR COUPLING section");
  }

  // *******************************************************************
  // other cases
  // *******************************************************************
  else
  {
    dserror("ERROR: IntegrateOverlap: Invalid case for 2D mortar coupling LM interpolation");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
CONTACT::CoCoupling2dManager::CoCoupling2dManager(DRT::Discretization& idiscret,
    int dim, bool quad, Teuchos::ParameterList& params,
    MORTAR::MortarElement* sele, std::vector<MORTAR::MortarElement*> mele) :
    MORTAR::Coupling2dManager(
        idiscret,
        dim,
        quad,
        params,
        sele,
        mele),
    stype_(DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params,"STRATEGY"))
{
  // evaluate coupling
  EvaluateCoupling();

  return;
}


/*----------------------------------------------------------------------*
 |  get communicator  (public)                               farah 01/13|
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::CoCoupling2dManager::Comm() const
{
  return idiscret_.Comm();
}

/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                  farah 10/14|
 *----------------------------------------------------------------------*/
bool CONTACT::CoCoupling2dManager::EvaluateCoupling()
{
  if(MasterElements().size() == 0)
    return false;

  // decide which type of coupling should be evaluated
  INPAR::MORTAR::AlgorithmType algo =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if(algo==INPAR::MORTAR::algorithm_mortar || algo==INPAR::MORTAR::algorithm_gpts)
    IntegrateCoupling();

  //*********************************
  // Error
  //*********************************
  else
    dserror("ERROR: chose contact algorithm not supported!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar coupling pairs                           Popp 03/09 |
 *----------------------------------------------------------------------*/
void CONTACT::CoCoupling2dManager::IntegrateCoupling()
{
  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (IntType() == INPAR::MORTAR::inttype_segments)
  {
      // loop over all master elements associated with this slave element
    for (int m = 0; m < (int) MasterElements().size(); ++m)
    {
      // create Coupling2d object and push back
      Coupling().push_back(
          Teuchos::rcp(
              new CoCoupling2d(idiscret_, dim_, quad_, imortar_,
                  SlaveElement(), MasterElement(m))));

      // project the element pair
      Coupling()[m]->Project();

      // check for element overlap
      Coupling()[m]->DetectOverlap();
    }

    // calculate consistent dual shape functions for this element
    ConsistDualShape();

    // do mortar integration
    for (int m = 0; m < (int) MasterElements().size(); ++m)
      Coupling()[m]->IntegrateOverlap();

    // free memory of consistent dual shape function coefficient matrix
    SlaveElement().MoData().ResetDualShape();
    SlaveElement().MoData().ResetDerivDualShape();
  }
  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (IntType() == INPAR::MORTAR::inttype_elements
      || IntType() == INPAR::MORTAR::inttype_elements_BS)
  {
    if ((int) MasterElements().size() == 0)
      return;

    // create an integrator instance with correct NumGP and Dim
    Teuchos::RCP<CONTACT::CoIntegrator> integrator = Teuchos::null;
    if (stype_==INPAR::CONTACT::solution_augmented)
      integrator = Teuchos::rcp(new CONTACT::AugmentedIntegrator(imortar_,SlaveElement().Shape(),
          Comm(),Teuchos::null));
    else
      integrator=Teuchos::rcp(new CONTACT::CoIntegrator(imortar_,SlaveElement().Shape(),Comm()));

    // *******************************************************************
    // different options for mortar integration
    // *******************************************************************
    // (1) no quadratic element(s) involved -> linear LM interpolation
    // (2) quadratic element(s) involved -> quadratic LM interpolation
    // (3) quadratic element(s) involved -> linear LM interpolation
    // (4) quadratic element(s) involved -> piecew. linear LM interpolation
    // *******************************************************************
    INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

    // *******************************************************************
    // cases (1), (2) and (3)
    // *******************************************************************
    if (!Quad() || (Quad() && lmtype == INPAR::MORTAR::lagmult_quad)
        || (Quad() && lmtype == INPAR::MORTAR::lagmult_lin))
    {
      // Test whether projection from slave to master surface is feasible --> important for dual LM Fnc.
      // Contact_interface.cpp --> AssembleG
      for (int m = 0; m < (int) MasterElements().size(); ++m)
      {
        // create CoCoupling2d object and push back
        Coupling().push_back(
            Teuchos::rcp(
                new CoCoupling2d(idiscret_, dim_, quad_, imortar_,
                    SlaveElement(), MasterElement(m))));

        // project the element pair
        Coupling()[m]->Project();
      }

      //Bool for identification of boundary elements
      bool boundary_ele = false;

      // ***********************************************************
      //                  START INTEGRATION !!!                   //
      // ***********************************************************
      if (stype_==INPAR::CONTACT::solution_augmented)
        Teuchos::rcp_dynamic_cast<CONTACT::AugmentedIntegrator>(integrator)->IntegrateDerivEle2D(SlaveElement(),
            MasterElements(),&boundary_ele);
      else
        integrator->IntegrateDerivEle2D(SlaveElement(), MasterElements(),&boundary_ele);
      // ***********************************************************
      //                   END INTEGRATION !!!                    //
      // ***********************************************************

      if (IntType() == INPAR::MORTAR::inttype_elements_BS
          and boundary_ele == true)
      {
        // switch, if consistent boundary modification chosen
        if (DRT::INPUT::IntegralValue<int>(imortar_, "LM_DUAL_CONSISTENT")
            == true && ShapeFcn() != INPAR::MORTAR::shape_standard // so for petrov-Galerkin and dual
            )
        {
          // loop over all master elements associated with this slave element
          for (int m = 0; m < (int) MasterElements().size(); ++m)
          {
            // create Coupling2d object and push back
            Coupling().push_back(Teuchos::rcp(new CoCoupling2d(idiscret_, dim_, quad_, imortar_,
                SlaveElement(), MasterElement(m))));

            // project the element pair
            Coupling()[m]->Project();

            // check for element overlap
            Coupling()[m]->DetectOverlap();
          }

          // calculate consistent dual shape functions for this element
          ConsistDualShape();

          // do mortar integration
          for (int m = 0; m < (int) MasterElements().size(); ++m)
            Coupling()[m]->IntegrateOverlap();

          // free memory of consistent dual shape function coefficient matrix
          SlaveElement().MoData().ResetDualShape();
          SlaveElement().MoData().ResetDerivDualShape();
        }

        // segment-based integration for boundary elements
        else
        {
          // loop over all master elements associated with this slave element
          for (int m = 0; m < (int) MasterElements().size(); ++m)
          {
            // create CoCoupling2d object and push back
            Coupling().push_back(
                Teuchos::rcp(
                    new CoCoupling2d(idiscret_, dim_, quad_, imortar_,
                        SlaveElement(), MasterElement(m))));

            // project the element pair
            Coupling()[m]->Project();

            // check for element overlap
            Coupling()[m]->DetectOverlap();

            // integrate the element overlap
            Coupling()[m]->IntegrateOverlap();
          }
        }
      }
      else
      {
        // nothing to do...
      }
    }
    // *******************************************************************
    // case (4)
    // *******************************************************************
    else if (Quad() && lmtype == INPAR::MORTAR::lagmult_pwlin)
    {
      dserror("ERROR: Piecewise linear LM not (yet?) implemented in 2D");
    }

    // *******************************************************************
    // undefined case
    // *******************************************************************
    else if (Quad() && lmtype == INPAR::MORTAR::lagmult_undefined)
    {
      dserror("Lagrange multiplier interpolation for quadratic elements undefined\n"
              "If you are using 2nd order mortar elements, you need to specify LM_QUAD in MORTAR COUPLING section");
    }

    // *******************************************************************
    // other cases
    // *******************************************************************
    else
    {
      dserror(
          "ERROR: Integrate: Invalid case for 2D mortar coupling LM interpolation");
    }
  }
  //**********************************************************************
  // INVALID TYPE OF NUMERICAL INTEGRATION
  //**********************************************************************
  else
  {
    dserror("ERROR: Invalid type of numerical integration");
  }
}


/*----------------------------------------------------------------------*
 |  Calculate dual shape functions                           seitz 07/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoCoupling2dManager::ConsistDualShape()
{

  // For standard shape functions no modification is necessary
  // A switch earlier in the process improves computational efficiency
  if (ShapeFcn() == INPAR::MORTAR::shape_standard)
    return;

  INPAR::MORTAR::ConsistentDualType consistent=DRT::INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(imortar_,"LM_DUAL_CONSISTENT");
  if (consistent==INPAR::MORTAR::consistent_none)
    return;

  // Consistent modification only for linear LM interpolation or NURBS
  if (Quad()==true && consistent!=INPAR::MORTAR::consistent_none && !SlaveElement().IsNurbs())
    dserror("Consistent dual shape functions in boundary elements only for linear LM interpolation");

  // do nothing if there are no coupling pairs
  if (Coupling().size()==0)
    return;

  const int nnodes = SlaveElement().NumNode();
  const int ndof = 2;

  int linsize = 0;
  for (int i=0;i<nnodes;++i)
  {
    CoNode* cnode = dynamic_cast<CoNode*> (SlaveElement().Nodes()[i]);
    linsize += cnode->GetLinsize();
  }

  int mnodes = 0;
  for (int m=0;m<(int)Coupling().size();++m)
    mnodes += MasterElements()[m]->NumNode();

  // detect entire overlap
  double ximin = 1.0;
  double ximax = -1.0;
  GEN::pairedvector<int, double> dximin(linsize + ndof*mnodes);
  GEN::pairedvector<int,double> dximax(linsize + ndof*mnodes);

  // loop over all master elements associated with this slave element
  for (int m=0;m<(int)Coupling().size();++m)
  {
    double sxia = Coupling()[m]->XiProj()[0];
    double sxib = Coupling()[m]->XiProj()[1];
    double mxia = Coupling()[m]->XiProj()[2];
    double mxib = Coupling()[m]->XiProj()[3];

    // no overlap for this slave-master pair --> continue with next pair
    if (sxia==0.0 && sxib==0.0) continue;

    // for contact we need the derivatives as well
    bool startslave = false;
    bool endslave = false;

    if (sxia==-1.0) startslave = true;
    else startslave = false;
    if (sxib==1.0) endslave = true;
    else endslave = false;

    // create an integrator for this segment
    CONTACT::CoIntegrator integrator(imortar_,SlaveElement().Shape(),Comm());

    std::vector<GEN::pairedvector<int,double> > ximaps(4,linsize + ndof*mnodes);
    // get directional derivatives of sxia, sxib, mxia, mxib
    integrator.DerivXiAB2D(SlaveElement(),sxia,sxib,MasterElement(m),mxia,mxib,ximaps,startslave,endslave,linsize);

    // get element contact integration area
    // and for contact derivatives of beginning and end
    if ((sxia!=0.0 || sxib!=0.0) && (sxia>=-1.0 && sxia<=1.0) && (sxib>=-1.0 && sxib<=1.0))
    {
      if (sxia<ximin && sxia>=-1. && sxia<=1.)
      {
        ximin=sxia;
        dximin=ximaps[0];
      }
      if (sxib>ximax && sxib>=-1. && sxib<=1.)
      {
        ximax=sxib;
        dximax=ximaps[1];
      }
    }
  }

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  // no overlap: the applied dual shape functions don't matter, as the integration domain is void
  if ((ximax==-1.0 && ximin==1.0) || (ximax-ximin<4.*MORTARINTLIM))
    return;

  // fully projecting element: no modification necessary
  if (ximin==-1.0 && ximax==1.0)
    return;

  // calculate consistent dual schape functions (see e.g. Cichosz et.al.:
  // Consistent treatment of boundaries with mortar contact formulations, CMAME 2010

  // Dual shape functions coefficient matrix
  LINALG::SerialDenseMatrix ae(nnodes,nnodes,true);

  // store derivae into element
  SlaveElement().MoData().DerivDualShape() =
      Teuchos::rcp(new GEN::pairedvector<int,Epetra_SerialDenseMatrix>(linsize + 2*ndof*mnodes,0,Epetra_SerialDenseMatrix(nnodes,nnodes)));
  GEN::pairedvector<int,Epetra_SerialDenseMatrix>& derivae=*(SlaveElement().MoData().DerivDualShape());

  // compute entries to bi-ortho matrices me/de with Gauss quadrature
  MORTAR::ElementIntegrator integrator(SlaveElement().Shape());

  // prepare for calculation of dual shape functions
  LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
  LINALG::SerialDenseMatrix de(nnodes,nnodes,true);
  // two-dim arrays of maps for linearization of me/de
  std::vector<std::vector<GEN::pairedvector<int,double> > > derivme(nnodes,std::vector<GEN::pairedvector<int,double> >(nnodes,linsize + 2*ndof*mnodes));
  std::vector<std::vector<GEN::pairedvector<int,double> > > derivde(nnodes,std::vector<GEN::pairedvector<int,double> >(nnodes,linsize + 2*ndof*mnodes));

  LINALG::SerialDenseVector sval(nnodes);
  LINALG::SerialDenseMatrix sderiv(nnodes,1,true);
  LINALG::SerialDenseMatrix ssecderiv(nnodes,1);

  for (int gp=0;gp<integrator.nGP();++gp)
  {
    // coordinates and weight
    double eta[2] =
    { integrator.Coordinate(gp,0), 0.0};
    double wgt = integrator.Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] =
    { 0.0, 0.0};
    sxi[0] = 0.5*(1.0-eta[0])*ximin + 0.5*(1.0+eta[0])*ximax;

    // evaluate trace space shape functions
    SlaveElement().EvaluateShape(sxi,sval,sderiv,nnodes);
    SlaveElement().Evaluate2ndDerivShape(sxi,ssecderiv,nnodes);

    // evaluate the two slave side Jacobians
    double dxdsxi = SlaveElement().Jacobian(sxi);
    double dsxideta = -0.5*ximin + 0.5*ximax;

    // evaluate linearizations *******************************************
    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] =
    { 0.0, 0.0};
    dynamic_cast<CONTACT::CoElement&> (SlaveElement()).DJacDXi(djacdxi,sxi,ssecderiv);
    double dxdsxidsxi=djacdxi[0]; // only 2D here

    // evalute the GP slave coordinate derivatives
    GEN::pairedvector<int,double> dsxigp(linsize + ndof*mnodes);
    for (_CI p=dximin.begin();p!=dximin.end();++p)
      dsxigp[p->first] += 0.5*(1-eta[0])*(p->second);
    for (_CI p=dximax.begin();p!=dximax.end();++p)
      dsxigp[p->first] += 0.5*(1+eta[0])*(p->second);

    // evaluate the Jacobian derivative
    GEN::pairedvector<int,double> derivjac(SlaveElement().NumNode()*Dim());
    SlaveElement().DerivJacobian(sxi,derivjac);

    // integrate dual shape matrices de, me and their linearizations
    for (int j=0; j<nnodes; ++j)
    {
      double fac;
      // de and linearization
      de(j,j) += wgt * sval[j] * dxdsxi*dsxideta;

      // (1) linearization of slave gp coordinates in ansatz function j for derivative of de
      fac=wgt * sderiv(j,0) * dxdsxi*dsxideta;
      for (_CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
        derivde[j][j][p->first] += fac*(p->second);

      // (2) linearization dsxideta - segment end coordinates
      fac=0.5*wgt*sval[j]*dxdsxi;
      for (_CI p=dximin.begin();p!=dximin.end();++p)
        derivde[j][j][p->first] -= fac*(p->second);
      fac=0.5*wgt*sval[j]*dxdsxi;
      for (_CI p=dximax.begin();p!=dximax.end();++p)
        derivde[j][j][p->first] += fac*(p->second);

      // (3) linearization dxdsxi - slave GP jacobian
      fac=wgt*sval[j]*dsxideta;
      for (_CI p=derivjac.begin();p!=derivjac.end();++p)
        derivde[j][j][p->first] += fac*(p->second);

      // (4) linearization dxdsxi - slave GP coordinates
      fac=wgt*sval[j]*dsxideta*dxdsxidsxi;
      for (_CI p=dsxigp.begin();p!=dsxigp.end();++p)
        derivde[j][j][p->first] += fac*(p->second);

      // me and linearization
      for (int k=0; k<nnodes; ++k)
      {
        me(j,k) += wgt * sval[j] * sval[k] *dxdsxi*dsxideta;

        // (1) linearization of slave gp coordinates in ansatz function for derivative of me
        fac=wgt * sval [k] * dxdsxi*dsxideta * sderiv(j,0);
        for (_CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
        {
          derivme[j][k][p->first] += fac*(p->second);
          derivme[k][j][p->first] += fac*(p->second);
        }

        // (2) linearization dsxideta - segment end coordinates
        fac=0.5*wgt*sval[j]*sval[k]*dxdsxi;
        for (_CI p=dximin.begin();p!=dximin.end();++p)
          derivme[j][k][p->first] -= fac*(p->second);
        fac=0.5*wgt*sval[j]*sval[k]*dxdsxi;
        for (_CI p=dximax.begin();p!=dximax.end();++p)
          derivme[j][k][p->first] += fac*(p->second);

        // (3) linearization dxdsxi - slave GP jacobian
        fac=wgt*sval[j]*sval[k]*dsxideta;
        for (_CI p=derivjac.begin();p!=derivjac.end();++p)
          derivme[j][k][p->first] += fac*(p->second);

        // (4) linearization dxdsxi - slave GP coordinates
        fac=wgt*sval[j]*sval[k]*dsxideta*dxdsxidsxi;
        for (_CI p=dsxigp.begin();p!=dsxigp.end();++p)
          derivme[j][k][p->first] += fac*(p->second);
      }
    }
  }
  // build ae matrix
  // invert bi-ortho matrix me
  LINALG::SerialDenseMatrix meinv = LINALG::InvertAndMultiplyByCholesky(me,de,ae);

  // build linearization of ae and store in derivdual
  // (this is done according to a quite complex formula, which
  // we get from the linearization of the biorthogonality condition:
  // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )

  // loop over all entries of ae (index i,j)
  for (int i=0;i<nnodes;++i)
  {
    for (int j=0;j<nnodes;++j)
    {
      // compute Lin(Ae) according to formula above
      for (int l=0;l<nnodes;++l)// loop over sum l
      {
        // part1: Lin(De)*Inv(Me)
        for (_CI p=derivde[i][l].begin();p!=derivde[i][l].end();++p)
          derivae[p->first](i,j) += meinv(l,j)*(p->second);

        // part2: Ae*Lin(Me)*Inv(Me)
        for (int k=0;k<nnodes;++k)// loop over sum k
          for (_CI p=derivme[k][l].begin();p!=derivme[k][l].end();++p)
            derivae[p->first](i,j) -= ae(i,k)*meinv(l,j)*(p->second);
      }
    }
  }

  // store ae matrix in slave element data container
  SlaveElement().MoData().DualShape() = Teuchos::rcp(new LINALG::SerialDenseMatrix(ae));

  return;
}
