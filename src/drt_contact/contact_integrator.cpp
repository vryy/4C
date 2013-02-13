/*!----------------------------------------------------------------------
\file contact_integrator.cpp
\brief A class to perform integrations of Mortar matrices on the overlap
       of two MortarElements in 1D and 2D (derived version for contact)

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_integrator.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 08/08|
 *----------------------------------------------------------------------*/
CONTACT::CoIntegrator::CoIntegrator(DRT::Element::DiscretizationType eletype) :
MORTAR::MortarIntegrator(eletype)
{
  // empty constructor
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 07/09|
 *----------------------------------------------------------------------*/
CONTACT::CoIntegrator::CoIntegrator(const INPAR::MORTAR::ShapeFcn shapefcn,
                                    DRT::Element::DiscretizationType eletype) :
MORTAR::MortarIntegrator(shapefcn,eletype)
{
  // empty constructor
    
  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize D on a slave element (2D / 3D)    popp 02/09|
 |  This method integrates the element D matrix and stores it in dseg.  |
 |  Moreover, derivatives LinD are built and stored directly into the   |
 |  adajcent nodes. (Thus this method combines EVERYTHING before done   |
 |  aperarately in IntegrateD and DerivD!)                              |
 |  ********** modified version: responds to shapefcn_ ***    popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivSlave2D3D(
      MORTAR::MortarElement& sele, double* sxia, double* sxib,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg)
{
  //**********************************************************************
  dserror("ERROR: IntegrateDerivSlave2D3D method is outdated!");
  //**********************************************************************

  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSlave2D3D called without specific shape function defined!");
    
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: IntegrateDerivSlave2D3D called on a non-slave MortarElement!");
  if ((sxia[0]<-1.0) || (sxia[1]<-1.0) || (sxib[0]>1.0) || (sxib[1]>1.0))
    dserror("ERROR: IntegrateDerivSlave2D3D called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();
  int ndof = Dim();
  int ncol = nrow;

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  // get slave element nodes themselves for GP normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivSlave2D3D: Null pointer!");

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except line2 (1D) & tri3 (2D)
  bool duallin = false;
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
  {
    if (sele.Shape()!=MORTAR::MortarElement::line2 && sele.Shape()!=MORTAR::MortarElement::tri3)
    {
      duallin = true;
      sele.DerivShapeDual(dualmap);
    }
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), 0.0};
    if (Dim()==3) eta[1] = Coordinate(gp,1);
    double wgt = Weight(gp);

    // evaluate trace space and Lagrange multiplier shape functions
    sele.EvaluateShape(eta,val,deriv,nrow);
    sele.EvaluateShapeLagMult(shapefcn_,eta,lmval,lmderiv,nrow);

    // evaluate the Jacobian det
    double dxdsxi = sele.Jacobian(eta);

    // compute element D matrix ******************************************
    // loop over all dtemp matrix entries
    // nrow represents the Lagrange multipliers !!!
    // ncol represents the dofs !!!
    for (int j=0;j<nrow*ndof;++j)
    {
      for (int k=0;k<ncol*ndof;++k)
      {
        int jindex = (int)(j/ndof);
        int kindex = (int)(k/ndof);

        // multiply the two shape functions
        double prod = lmval[jindex]*val[kindex];

        // isolate the dseg entries to be filled
        // (both the main diagonal and every other secondary diagonal)
        // and add current Gauss point's contribution to dseg
        if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          (*dseg)(j,k) += prod*dxdsxi*wgt;
      }
    }
    // compute element D matrix ******************************************

    // evaluate the Jacobian derivative
    std::map<int,double> derivjac;
    sele.DerivJacobian(eta,derivjac);
        
    // compute element D linearization ***********************************
    
    // here we have to adapt the alogrithm to support arbitrary shape functions
    // because we can't rely on the diagonality of D anymore
    // this was built analogous to building derivM from IntegrateDerivSegment2D
    
    for (int i=0;i<nrow;++i)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[i]);
      if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSlave: Null pointer!");
      bool bound = mymrtrnode->IsOnBound();

      //******************************************************************
      // standard case (node i requires no edge node modification)
      //******************************************************************
      // this is the standard case
      if (!bound)
      {
        // loop over all slave nodes again
        for (int k=0; k<nrow; ++k)
        {
          // slave node ID
          int sgid = sele.Nodes()[k]->Id();

          // get the correct map as a reference
          std::map<int,double>& ddmap_ik = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];
          
          // multiply the corresponding two shape functions
          double prod = wgt*val[k]*lmval[i];

          // derivative of Jacobian
          for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
            ddmap_ik[p->first] += prod*(p->second);
          
          // derivative of dual shape function
          if (duallin)
          {
            for (int j=0;j<nrow;++j)
            {
              double fac = wgt*val[j]*val[k]*dxdsxi;
              for (CI p=dualmap[i][j].begin();p!=dualmap[i][j].end();++p)
                ddmap_ik[p->first] += fac*(p->second);
            }
          }
        } 
      }

      //******************************************************************
      // edge node case (node i requires edge node modification)
      //******************************************************************
      // all nodes on edges cause additional entries in derivM
      // this is a special case depending on the bound-flag of mrtrnode
      // at the moment this is only possible for dual shape functions in 2D
      else
      {
        // check the shape function type
        if (shapefcn_ == INPAR::MORTAR::shape_standard)
          dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for standard shape functions");
                
        //check for problem dimension
        if (Dim()==3)
          dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for 3D");

        // get gid of current boundary node
        int bgid = mymrtrnode->Id();

        // loop over other nodes (interior nodes)
        for (int k=0;k<nrow;++k)
        {
          MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[k]);
          if (!mymrtrnode2) dserror("ERROR: IntegrateAndDerivSlave: Null pointer!");
          bool bound2 = mymrtrnode2->IsOnBound();
          if (bound2) continue;
          std::map<int,double>& nodemmap = static_cast<CONTACT::CoNode*>(mymrtrnode2)->CoData().GetDerivM()[bgid];

          // derivative of Jacobian
          double fac = wgt*val[i]*lmval[k];
          for (CI p=derivjac.begin();p!=derivjac.end();++p)
            nodemmap[p->first] -= fac*(p->second);
          
          // derivative of dual shape functions
          if (duallin)
          {
            for (int j=0;j<nrow;++j)
            {
              LINALG::SerialDenseVector vallin(nrow-1);
              LINALG::SerialDenseMatrix derivlin(nrow-1,1);
              if (i==0) sele.ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge0,eta,vallin,derivlin);
              else if (i==1) sele.ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge1,eta,vallin,derivlin);
              double fac = wgt*val[i]*vallin[j]*dxdsxi;
              for (CI p=dualmap[k][j].begin();p!=dualmap[k][j].end();++p)
                nodemmap[p->first] -= fac*(p->second);
            }
          }
        }
      }
    }
    // compute element D linearization ***********************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master overlap (2D)  popp 02/09|
 |  This method integrates the overlap M matrix and weighted gap g~     |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM, IntegrateG, DerivM and DerivG!)                         |
 |  Also wear is integrated.                                            |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivSegment2D(
     MORTAR::MortarElement& sele, double& sxia, double& sxib,
     MORTAR::MortarElement& mele, double& mxia, double& mxib,
     INPAR::MORTAR::LagMultQuad lmtype,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg,
     Teuchos::RCP<Epetra_SerialDenseVector> wseg)
{ 
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  // Petrov-Galerkin approach for LM not yet implemented for quadratic FE
  if (sele.Shape()==MORTAR::MortarElement::line3 && shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("ERROR: Petrov-Galerkin approach not yet implemented for quadratic FE interpolation");

  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

  // contact with wear
  bool wear = false;  
  if(DRT::Problem::Instance()->ContactDynamicParams().get<double>("WEARCOEFF")!= 0.0)
    wear = true;
  
  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
  
  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,1);
  
  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  mele.GetNodalCoords(mcoord);
  
  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;
  
  if(wear)
  { 
    scoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalCoordsOld(*scoordold);
    mcoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,mele.NumNode()));
    mele.GetNodalCoordsOld(*mcoordold);
    lagmult = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalLagMult(*lagmult);
  }
  
  // map iterator
  typedef std::map<int,double>::const_iterator CI;

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
  if (sxib==1.0) endslave = true;
  else           endslave = false;

  // get directional derivatives of sxia, sxib, mxia, mxib
  std::vector<std::map<int,double> > ximaps(4);
  DerivXiAB2D(sele,sxia,sxib,mele,mxia,mxib,ximaps,startslave,endslave);

  // prepare directional derivative of dual shape functions
  // this is only necessary for quadratic dual shape functions in 2D
  bool duallin = false;
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape()==MORTAR::MortarElement::line3))
  {
    duallin=true;
    sele.DerivShapeDual(dualmap);
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  if (lmtype == INPAR::MORTAR::lagmult_lin_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false; // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
  }

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
    MORTAR::MortarProjector projector(2);
    projector.ProjectGaussPoint(sele,sxi,mele,mxi);

    // simple version (no Gauss point projection)
    // double test[2] = {0.0, 0.0};
    // test[0] = 0.5*(1-eta[0])*mxib + 0.5*(1+eta[0])*mxia;
    // mxi[0]=test[0];

    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
      dserror("ERROR: IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d",mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (linlm)
      sele.EvaluateShapeLagMultLin(shapefcn_,sxi,lmval,lmderiv,nrow);
    else
      sele.EvaluateShapeLagMult(shapefcn_,sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);
    
    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    // compute segment D/M matrix ****************************************
    // standard shape functions
    if (shapefcn_ == INPAR::MORTAR::shape_standard)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to mseg and dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
          }
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to mseg and dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*dseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
          }
        }
      } // nrow*ndof loop
    }
    
    // dual shape functions
    else if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
    { 
      // loop over all mseg matrix entries
      // nrow represents the slave Lagrange multipliers !!!
      // ncol represents the master dofs !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0;j<nrow*ndof;++j)
      {
        // integrate mseg and dseg (NO boundary modifiation)
        for (int k=0;k<ncol*ndof;++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal for mseg)
          // (only the main diagonal for dseg)
          // and add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
            if (!bound) (*dseg)(j,j) += prod*dxdsxi*dsxideta*wgt;
          }
        }

        // integrate dseg (boundary modification)
        if (bound)
        {
          MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[(int)(j/ndof)]);
          if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
          bool j_boundnode = mymrtrnode->IsOnBound();
          
          for (int k=0;k<nrow*ndof;++k)
          {
            MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[(int)(k/ndof)]);
            if (!mymrtrnode2) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
            bool k_boundnode = mymrtrnode2->IsOnBound();
            
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);
            
            // do not assemble off-diagonal terms if j,k are both non-boundary nodes
            if (!j_boundnode && !k_boundnode && (jindex!=kindex)) continue;
              
            // multiply the two shape functions 
            double prod = lmval[jindex]*sval[kindex];
    
            // isolate the dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to dseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
          }
        }
      }
    } // shapefcn_ switch
    // compute segment D/M matrix ****************************************

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);
    
    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
      gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
      gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

      sgpx[0]+=sval[i]*scoord(0,i);
      sgpx[1]+=sval[i]*scoord(1,i);
      sgpx[2]+=sval[i]*scoord(2,i);
    }

    // normalize interpolated GP normal back to length 1.0 !!!
    double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
    if (length<1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

    for (int i=0;i<3;++i)
      gpn[i]/=length;

    // build interpolation of master GP coordinates
    double mgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<ncol;++i)
    {
      mgpx[0]+=mval[i]*mcoord(0,i);
      mgpx[1]+=mval[i]*mcoord(1,i);
      mgpx[2]+=mval[i]*mcoord(2,i);
    }

    // build gap function at current GP
    double gap = 0.0;
    for (int i=0;i<3;++i)
      gap+=(mgpx[i]-sgpx[i])*gpn[i];
    
    // evaluate wear at current GP
    // not including wearcoefficient
    double wearval = 0;
    if(wear)
    {
      double gpt[3] = {0.0,0.0,0.0};
      double gplm[3] = {0.0, 0.0, 0.0};
      double sgpjump[3] = {0.0, 0.0, 0.0};
      for (int i=0;i<nrow;++i)
      {
         CONTACT::CoNode* myconode = static_cast<CONTACT::CoNode*> (mynodes[i]);
         gpt[0]+=sval[i]*myconode->CoData().txi()[0];
         gpt[1]+=sval[i]*myconode->CoData().txi()[1];
         gpt[2]+=sval[i]*myconode->CoData().txi()[2];
         
         sgpjump[0]+=sval[i]*(scoord(0,i)-(*scoordold)(0,i));
         sgpjump[1]+=sval[i]*(scoord(1,i)-(*scoordold)(1,i));
         sgpjump[2]+=sval[i]*(scoord(2,i)-(*scoordold)(2,i));
       
         gplm[0]+=lmval[i]*(*lagmult)(0,i);
         gplm[1]+=lmval[i]*(*lagmult)(1,i);
         gplm[2]+=lmval[i]*(*lagmult)(2,i);
      }
      
      // normalize interpolated GP tangent back to length 1.0 !!!
      length = sqrt(gpt[0]*gpt[0]+gpt[1]*gpt[1]+gpt[2]*gpt[2]);
      if (length<1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");
      
      for (int i=0;i<3;i++)
        gpt[i]/=length;
      
      // interpolation of master GP jumps (relative displacement increment)
      double mgpjump[3] = {0.0, 0.0, 0.0};
      for (int i=0;i<ncol;++i)
      {
        mgpjump[0]+=mval[i]*(mcoord(0,i)-(*mcoordold)(0,i));
        mgpjump[1]+=mval[i]*(mcoord(1,i)-(*mcoordold)(1,i));
        mgpjump[2]+=mval[i]*(mcoord(2,i)-(*mcoordold)(2,i));
      }
        
      // jump
      double jump[3] = {0.0, 0.0, 0.0};
      jump[0] = sgpjump[0] - mgpjump[0]; 
      jump[1] = sgpjump[1] - mgpjump[1]; 
      jump[2] = sgpjump[2] - mgpjump[2]; 
        
      // evaluate wear   
      // normal contact stress 
      for (int i=0;i<3;++i)
        wearval+=gpn[i]*gplm[i];
 
      // value of relative tangential jump
      double jumpval = 0.0;
      for (int i=0;i<3;++i)
        jumpval+=gpt[i]*jump[i];
      
      // product
      wearval=abs(wearval)*abs(jumpval);
    }
    
    // evaluate linearizations *******************************************
    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    static_cast<CONTACT::CoElement&>(sele).DJacDXi(djacdxi,sxi,ssecderiv);
    double dxdsxidsxi=djacdxi[0]; // only 2D here

    // evalute the GP slave coordinate derivatives
    std::map<int,double> dsxigp;
    for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
      dsxigp[p->first] += 0.5*(1-eta[0])*(p->second);
    for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
      dsxigp[p->first] += 0.5*(1+eta[0])*(p->second);

    // evalute the GP master coordinate derivatives
    std::map<int,double> dmxigp;
    DerivXiGP2D(sele,mele,sxi[0],mxi[0],dsxigp,dmxigp);

    // simple version (no Gauss point projection)
    //for (CI p=ximaps[2].begin();p!=ximaps[2].end();++p)
    //  dmxigp[p->first] += 0.5*(1+eta[0])*(p->second);
    //for (CI p=ximaps[3].begin();p!=ximaps[3].end();++p)
    //  dmxigp[p->first] += 0.5*(1-eta[0])*(p->second);

    // evaluate the Jacobian derivative
    std::map<int,double> derivjac;
    sele.DerivJacobian(sxi,derivjac);

    // evaluate the GP gap function derivatives
    std::map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
      if (!smrtrnodes[i]) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
      if (!mmrtrnodes[i]) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    std::map<int,double> dmap_nxsl_gp;
    std::map<int,double> dmap_nysl_gp;

    for (int i=0;i<nrow;++i)
    {
      std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
      {
        double valx =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
      }
    }

    // build directional derivative of slave GP normal (unit)
    std::map<int,double> dmap_nxsl_gp_unit;
    std::map<int,double> dmap_nysl_gp_unit;

    double ll = length*length;
    double sxsx = gpn[0]*gpn[0]*ll;
    double sxsy = gpn[0]*gpn[1]*ll;
    double sysy = gpn[1]*gpn[1]*ll;

    for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    {
      dmap_nxsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsx*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
    }

    for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
    {
      dmap_nysl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysy*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
    }

    // add everything to dgapgp
    for (CI p=dmap_nxsl_gp_unit.begin();p!=dmap_nxsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

    for (CI p=dmap_nysl_gp_unit.begin();p!=dmap_nysl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

    for (int z=0;z<nrow;++z)
    {
      for (int k=0;k<2;++k)
      {
        dgapgp[smrtrnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * smrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<2;++k)
      {
        dgapgp[mmrtrnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mmrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************
    
    // compute segment gap vector ****************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!
    for (int j=0;j<nrow;++j)
    {
      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
        prod = sval[j]*gap;
      // usual standard or dual LM approach
      else
        prod = lmval[j]*gap;

      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*dxdsxi*dsxideta*wgt;
    }
    // compute segment gap vector ****************************************
    
    // compute segment wear vector ***************************************
    // loop over all wseg vector entries
    // nrow represents the slave side dofs !!!
    if (wear)
    {
      for (int j=0;j<nrow;++j)
      {
        double prod = lmval[j]*wearval;

        // add current Gauss point's contribution to wseg
        (*wseg)(j) += prod*dxdsxi*dsxideta*wgt;
      }
    }

    // compute segment D/M linearization *********************************

    // **************** edge modification ********************************
    if (bound)
    {
      // check the shape function type (not really necessary because only dual shape functions arrive here)
      if (shapefcn_ == INPAR::MORTAR::shape_standard)
        dserror("ERROR: IntegrateDerivSegment2D: Edge node mod. called for standard shape functions");

      // check for Petrov-Galerkin interpolation
      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
        dserror("ERROR: IntegrateDerivSegment2D: Petrov-Galerkin and boundary modification not compatible");

      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
        bool boundnode = mymrtrnode->IsOnBound();
        int sgid = mymrtrnode->Id();
        std::map<int,double>& nodemap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];
        double fac = 0.0;

        //******************************************************************
        // standard case (node j is NO boundary node)
        //******************************************************************
        // only process the entry D_jj, the entried D_jk will be moved to M_jk
        if (!boundnode)
        {
          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0;m<nrow;++m)
            {
              fac = wgt*sval[j]*sval[m]*dsxideta*dxdsxi;
              for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
                nodemap[p->first] += fac*(p->second);
            }

          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*lmderiv(j,0)*sval[j]*dsxideta*dxdsxi;
          for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
            nodemap[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmval[j]*sderiv(j,0)*dsxideta*dxdsxi;
          for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
            nodemap[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - segment end coordinates
          fac = wgt*lmval[j]*sval[j]*dxdsxi;
          for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
            nodemap[p->first] -= 0.5*fac*(p->second);
          for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
            nodemap[p->first] += 0.5*fac*(p->second);

          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*lmval[j]*sval[j]*dsxideta;
          for (CI p=derivjac.begin();p!=derivjac.end();++p)
            nodemap[p->first] += fac*(p->second);

          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*lmval[j]*sval[j]*dsxideta*dxdsxidsxi;
          for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
            nodemap[p->first] += fac*(p->second);
        }

        //******************************************************************
        // edge case (node j is a boundary node)
        //******************************************************************
        else
        {
          // get gid of current boundary node
          int bgid = mymrtrnode->Id();

          // loop over other nodes (interior nodes)
          for (int k=0;k<nrow;++k)
          {
            MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[k]);
            if (!mymrtrnode2) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
            bool boundnode2 = mymrtrnode2->IsOnBound();
            if (boundnode2) continue;
            std::map<int,double>& nodemmap = static_cast<CONTACT::CoNode*>(mymrtrnode2)->CoData().GetDerivM()[bgid];

            // (1) Lin(Phi) - dual shape functions
            if (duallin)
              for (int m=0;m<nrow;++m)
              {
                LINALG::SerialDenseVector vallin(nrow-1);
                LINALG::SerialDenseMatrix derivlin(nrow-1,1);
                if (j==0) sele.ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge0,sxi,vallin,derivlin);
                else if (j==1) sele.ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge1,sxi,vallin,derivlin);
                double fac = wgt*sval[j]*vallin[m]*dsxideta*dxdsxi;
                for (CI p=dualmap[k][m].begin();p!=dualmap[k][m].end();++p)
                  nodemmap[p->first] -= fac*(p->second);
              }

            // (2) Lin(Phi) - slave GP coordinates
            fac = wgt*lmderiv(k,0)*sval[j]*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
              nodemmap[p->first] -= fac*(p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*lmval[k]*sderiv(j,0)*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
              nodemmap[p->first] -= fac*(p->second);

            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*lmval[k]*sval[j]*dxdsxi;
            for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
              nodemmap[p->first] += 0.5*fac*(p->second);
            for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
              nodemmap[p->first] -= 0.5*fac*(p->second);

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*lmval[k]*sval[j]*dsxideta;
            for (CI p=derivjac.begin();p!=derivjac.end();++p)
              nodemmap[p->first] -= fac*(p->second);

            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*lmval[k]*sval[j]*dsxideta*dxdsxidsxi;
            for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
              nodemmap[p->first] -= fac*(p->second);
          }
        }
      }
    }

    // **************** no edge modification *****************************
    // (and LinM also for edge node modification case)

    // check for linear LM interpolation in quadratic FE
    if (linlm)
    {
      if (shapefcn_ != INPAR::MORTAR::shape_standard)
        dserror("ERROR: No linear dual LM interpolation for 2D quadratic FE");

      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
        dserror("ERROR: Petrov-Galerkin and linear LM for 2D quadratic FE not compatible");

      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
        bool jbound = mymrtrnode->IsOnBound();

        // node j is boundary node
        if (jbound)
        {
          // do nothing as respective D and M entries are zero anyway
        }

        // node j is NO boundary node
        else
        {
          // integrate LinM
          for (int k=0; k<ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*lmderiv(j,0)*mval[k]*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt*lmval[j]*mderiv(k, 0)*dsxideta*dxdsxi;
            for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*lmval[j]*mval[k]*dxdsxi;
            for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
              dmmap_jk[p->first] -= 0.5*fac*(p->second);
            for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
              dmmap_jk[p->first] += 0.5*fac*(p->second);

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*lmval[j]*mval[k]*dsxideta;
            for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*lmval[j]*mval[k]*dsxideta*dxdsxidsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);
          } // loop over master nodes

          // integrate LinD
          for (int k=0; k<nrow; ++k)
          {
            MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[k]);
            if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
            bool kbound = mymrtrnode2->IsOnBound();

            // global master node ID
            int sgid = mymrtrnode2->Id();
            double fac = 0.0;

            // node k is boundary node
            if (kbound)
            {
              // move entry to derivM (with minus sign)
              // get the correct map as a reference
              std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

              // (1) Lin(Phi) - dual shape functions
              // this vanishes here since there are no deformation-dependent dual functions

              // (2) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmderiv(j,0)*sval[k]*dsxideta*dxdsxi;
              for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              // (3) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmval[j]*sderiv(k, 0)*dsxideta*dxdsxi;
              for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              // (4) Lin(dsxideta) - segment end coordinates
              fac = wgt*lmval[j]*sval[k]*dxdsxi;
              for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
                dmmap_jk[p->first] += 0.5*fac*(p->second);
              for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
                dmmap_jk[p->first] -= 0.5*fac*(p->second);

              // (5) Lin(dxdsxi) - slave GP Jacobian
              fac = wgt*lmval[j]*sval[k]*dsxideta;
              for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              // (6) Lin(dxdsxi) - slave GP coordinates
              fac = wgt*lmval[j]*sval[k]*dsxideta*dxdsxidsxi;
              for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);
            }

            // node k is NO boundary node
            else
            {
              // get the correct map as a reference
              std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

              // (1) Lin(Phi) - dual shape functions
              // this vanishes here since there are no deformation-dependent dual functions

              // (2) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmderiv(j,0)*sval[k]*dsxideta*dxdsxi;
              for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              // (3) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmval[j]*sderiv(k, 0)*dsxideta*dxdsxi;
              for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              // (4) Lin(dsxideta) - segment end coordinates
              fac = wgt*lmval[j]*sval[k]*dxdsxi;
              for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
                ddmap_jk[p->first] -= 0.5*fac*(p->second);
              for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
                ddmap_jk[p->first] += 0.5*fac*(p->second);

              // (5) Lin(dxdsxi) - slave GP Jacobian
              fac = wgt*lmval[j]*sval[k]*dsxideta;
              for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              // (6) Lin(dxdsxi) - slave GP coordinates
              fac = wgt*lmval[j]*sval[k]*dsxideta*dxdsxidsxi;
              for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);
            }
          } // loop over slave nodes
        }
      }
    }

    // no linear LM interpolation for quadratic FE
    else
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

        int sgid = mymrtrnode->Id();

        // standard shape functions
        if (shapefcn_ == INPAR::MORTAR::shape_standard)
        {
          // integrate LinM
          for (int k=0; k<ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*lmderiv(j,0)*mval[k]*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt*lmval[j]*mderiv(k, 0)*dsxideta*dxdsxi;
            for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*lmval[j]*mval[k]*dxdsxi;
            for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
              dmmap_jk[p->first] -= 0.5*fac*(p->second);
            for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
              dmmap_jk[p->first] += 0.5*fac*(p->second);

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*lmval[j]*mval[k]*dsxideta;
            for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*lmval[j]*mval[k]*dsxideta*dxdsxidsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);
          } // loop over master nodes

          // integrate LinD
          for (int k=0; k<nrow; ++k)
          {
            // global slave node ID
            int sgid = sele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*lmderiv(j,0)*sval[k]*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*lmval[j]*sderiv(k, 0)*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*lmval[j]*sval[k]*dxdsxi;
            for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
              ddmap_jk[p->first] -= 0.5*fac*(p->second);
            for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
              ddmap_jk[p->first] += 0.5*fac*(p->second);

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*lmval[j]*sval[k]*dsxideta;
            for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*lmval[j]*sval[k]*dsxideta*dxdsxidsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
          } // loop over slave nodes
        }

        // dual shape functions
        else if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
        {
          // get the D-map as a reference
          std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // integrate LinM and LinD (NO boundary modification)
          for (int k=0; k<ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

            // (1) Lin(Phi) - dual shape functions
            for (int m=0; m<nrow; ++m)
            {
              fac = wgt*sval[m]*mval[k]*dsxideta*dxdsxi;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                if (!bound) ddmap_jk[p->first] += fac*(p->second);
              }
            }

            // (2) Lin(Phi) - slave GP coordinates
            fac = wgt*lmderiv(j, 0)*mval[k]*dsxideta*dxdsxi;

            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
            {
              dmmap_jk[p->first] += fac*(p->second);
              if (!bound) ddmap_jk[p->first] += fac*(p->second);
            }

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt*lmval(j, 0)*mderiv(k, 0)*dsxideta*dxdsxi;

            for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
            {
              dmmap_jk[p->first] += fac*(p->second);
              if (!bound) ddmap_jk[p->first] += fac*(p->second);
            }

            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*lmval[j]*mval[k]*dxdsxi;

            for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
            {
              dmmap_jk[p->first] -= 0.5*fac*(p->second);
              if (!bound) ddmap_jk[p->first] -= 0.5*fac*(p->second);
            }
            for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
            {
              dmmap_jk[p->first] += 0.5*fac*(p->second);
              if (!bound) ddmap_jk[p->first] += 0.5*fac*(p->second);
            }

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*lmval[j]*mval[k]*dsxideta;

            for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
            {
              dmmap_jk[p->first] += fac*(p->second);
              if (!bound) ddmap_jk[p->first] += fac*(p->second);
            }

            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*lmval[j]*mval[k]*dsxideta*dxdsxidsxi;

            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
            {
              dmmap_jk[p->first] += fac*(p->second);
              if (!bound) ddmap_jk[p->first] += fac*(p->second);
            }
          } // loop over master nodes
        } // shapefcn_ switch
      } // linlm or not
    }
    // compute segment D/M linearization *********************************

    // compute segment gap linearization *********************************
    for (int j=0;j<nrow;++j)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

      // switch if Petrov-Galerkin approach for LM is applied
      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // (1) Lin(Phi) - does not exist in gap for Petrov-Galerkin interpolation
        // as std shape functions are used here

        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*sderiv(j,0)*gap*dsxideta*dxdsxi;
        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*sval[j]*dsxideta*dxdsxi;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - segment end coordinates
        fac = wgt*sval[j]*gap*dxdsxi;
        for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
          dgmap[p->first] -= 0.5*fac*(p->second);
        for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
          dgmap[p->first] += 0.5*fac*(p->second);

        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*sval[j]*gap*dsxideta;
        for (CI p=derivjac.begin();p!=derivjac.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*sval[j]*gap*dsxideta*dxdsxidsxi;
        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dgmap[p->first] += fac*(p->second);
      }

      // the usual standard or dual LM interpolation 
      else
      {
        // (1) Lin(Phi) - dual shape functions
        if (shapefcn_ == INPAR::MORTAR::shape_dual)
        {
          for (int m=0;m<nrow;++m)
          {
            fac = wgt*sval[m]*gap*dsxideta*dxdsxi;
            for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
              dgmap[p->first] += fac*(p->second);
          }
        }

        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmderiv(j,0)*gap*dsxideta*dxdsxi;
        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*lmval[j]*dsxideta*dxdsxi;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - segment end coordinates
        fac = wgt*lmval[j]*gap*dxdsxi;
        for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
          dgmap[p->first] -= 0.5*fac*(p->second);
        for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
          dgmap[p->first] += 0.5*fac*(p->second);

        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*lmval[j]*gap*dsxideta;
        for (CI p=derivjac.begin();p!=derivjac.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*lmval[j]*gap*dsxideta*dxdsxidsxi;
        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dgmap[p->first] += fac*(p->second);
      }
    }
    // compute segment gap linearization *********************************
  }
  //**********************************************************************

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate and linearize                                  farah 01/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3D_Fast(
     MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg,
     bool *boundary_ele,
     bool *proj_)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType dt = meles[0]->Shape();
  DRT::Element::DiscretizationType dt_s = sele.Shape();

  // check input data
  for (int test=0;test<(int)meles.size();++test)
  {
  if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
    dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  }

  int msize   = meles.size();
  int nrow    = sele.NumNode();
  int nmnode  = meles[0]->NumNode();
  //int ncol    = msize*nmnode;
  int ndof    = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(nmnode);
  LINALG::SerialDenseMatrix mderiv(nmnode,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual) && (sele.Shape()!=MORTAR::MortarElement::tri3))
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  //************************************************
  //  Boundary_segmentation test -- HasProj() check
  //************************************************
  const Teuchos::ParameterList& smortar  = DRT::Problem::Instance()->MortarCouplingParams();

  INPAR::MORTAR::IntType integrationtype =
    DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(smortar,"INTTYPE");

  double sxi_test[2] = {0.0, 0.0};
  double alpha_test=0.0;
  bool proj_test=false;
  double glob_test[3] = {0.0, 0.0, 0.0};

  DRT::Node** mynodes_test = sele.Nodes();
  if (!mynodes_test) dserror("ERROR: HasProjStatus: Null pointer!");

  if (dt_s==DRT::Element::quad4 )//|| dt_s==DRT::Element::quad8 || dt_s==DRT::Element::quad9)
  {
    for (int s_test=0;s_test<4;++s_test)
    {
      if (s_test==0) {sxi_test[0]=-1.0;sxi_test[1]=-1.0;}
      else if (s_test==1){sxi_test[0]=-1.0;sxi_test[1]=1.0;}
      else if (s_test==2){sxi_test[0]=1.0;sxi_test[1]=-1.0;}
      else if (s_test==3){sxi_test[0]=1.0;sxi_test[1]=1.0;}

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(3);
        projector_test.ProjectGaussPoint3D(sele,sxi_test,*meles[bs_test],mxi_test,alpha_test);

        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi_test[0]>=-1.0 && mxi_test[1]>=-1.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0)
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else if(dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
        {
          if (mxi_test[0]>=0.0 && mxi_test[1]>=0.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0 && mxi_test[0]+mxi_test[1]<=1.0)
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }
            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  else if (dt_s==DRT::Element::quad9 )//|| dt_s==DRT::Element::quad8 || dt_s==DRT::Element::quad9)
  {
    for (int s_test=0;s_test<9;++s_test)
    {
      if (s_test==0) {sxi_test[0]=-1.0;sxi_test[1]=-1.0;}
      else if (s_test==1){sxi_test[0]=0.0;sxi_test[1]=-1.0;}
      else if (s_test==2){sxi_test[0]=1.0;sxi_test[1]=-1.0;}
      else if (s_test==3){sxi_test[0]=-1.0;sxi_test[1]=0.0;}
      else if (s_test==4){sxi_test[0]=0.0;sxi_test[1]=0.0;}
      else if (s_test==5){sxi_test[0]=1.0;sxi_test[1]=0.0;}
      else if (s_test==6){sxi_test[0]=-1.0;sxi_test[1]=1.0;}
      else if (s_test==7){sxi_test[0]=0.0;sxi_test[1]=1.0;}
      else if (s_test==8){sxi_test[0]=1.0;sxi_test[1]=1.0;}

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(3);
        projector_test.ProjectGaussPoint3D(sele,sxi_test,*meles[bs_test],mxi_test,alpha_test);

        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi_test[0]>=-1.0 && mxi_test[1]>=-1.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0)
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else if(dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
        {
          if (mxi_test[0]>=0.0 && mxi_test[1]>=0.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0 && mxi_test[0]+mxi_test[1]<=1.0)
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }
            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  else if (dt_s==DRT::Element::quad8 )//|| dt_s==DRT::Element::quad8 || dt_s==DRT::Element::quad9)
  {
    for (int s_test=0;s_test<8;++s_test)
    {
      if (s_test==0) {sxi_test[0]=-1.0;sxi_test[1]=-1.0;}
      else if (s_test==1){sxi_test[0]=0.0;sxi_test[1]=-1.0;}
      else if (s_test==2){sxi_test[0]=1.0;sxi_test[1]=-1.0;}
      else if (s_test==3){sxi_test[0]=-1.0;sxi_test[1]=0.0;}
      else if (s_test==4){sxi_test[0]=1.0;sxi_test[1]=0.0;}
      else if (s_test==5){sxi_test[0]=-1.0;sxi_test[1]=1.0;}
      else if (s_test==6){sxi_test[0]=0.0;sxi_test[1]=1.0;}
      else if (s_test==7){sxi_test[0]=1.0;sxi_test[1]=1.0;}

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(3);
        projector_test.ProjectGaussPoint3D(sele,sxi_test,*meles[bs_test],mxi_test,alpha_test);

        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi_test[0]>=-1.0 && mxi_test[1]>=-1.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0)
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else if(dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
        {
          if (mxi_test[0]>=0.0 && mxi_test[1]>=0.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0 && mxi_test[0]+mxi_test[1]<=1.0)
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }
            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  //TRI-ELE
  else if (dt_s==DRT::Element::tri3)
  {
    for (int s_test=0;s_test<3;++s_test)
    {
      if (s_test==0) {sxi_test[0]=0.0;sxi_test[1]=0.0;}
      else if (s_test==1){sxi_test[0]=1.0;sxi_test[1]=0.0;}
      else if (s_test==2){sxi_test[0]=0.0;sxi_test[1]=1.0;}

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(3);
        projector_test.ProjectGaussPoint3D(sele,sxi_test,*meles[bs_test],mxi_test,alpha_test);

        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi_test[0]>=-1.0 && mxi_test[1]>=-1.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0) //Falls Position auf Element
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
        {
          if (mxi_test[0]>=0.0 && mxi_test[1]>=0.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0 && mxi_test[0]+mxi_test[1]<=1.0) //Falls Position auf Element
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  else if (dt_s==DRT::Element::tri6)
  {
    for (int s_test=0;s_test<6;++s_test)
    {
      if (s_test==0) {sxi_test[0]=0.0;sxi_test[1]=0.0;}
      else if (s_test==1){sxi_test[0]=0.5;sxi_test[1]=0.0;}
      else if (s_test==2){sxi_test[0]=1.0;sxi_test[1]=0.0;}
      else if (s_test==3){sxi_test[0]=0.0;sxi_test[1]=0.5;}
      else if (s_test==4){sxi_test[0]=0.5;sxi_test[1]=0.5;}
      else if (s_test==5){sxi_test[0]=0.0;sxi_test[1]=1.0;}

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(3);
        projector_test.ProjectGaussPoint3D(sele,sxi_test,*meles[bs_test],mxi_test,alpha_test);

        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi_test[0]>=-1.0 && mxi_test[1]>=-1.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0) //Falls Position auf Element
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
        {
          if (mxi_test[0]>=0.0 && mxi_test[1]>=0.0 && mxi_test[0]<=1.0 && mxi_test[1]<=1.0 && mxi_test[0]+mxi_test[1]<=1.0) //Falls Position auf Element
          {
            //get hasproj
            sele.LocalToGlobal(sxi_test,glob_test,0);
            for (int ii=0;ii<sele.NumNode();++ii)
            {
              MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
              if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

              if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
                mycnode_test->HasProj()=true;
            }

            glob_test[0]=0.0;
            glob_test[1]=0.0;
            glob_test[2]=0.0;

            proj_test=true;
          }
        }
        else
        {
          dserror("Non valid element type for master discretization!");
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  else
  {
    dserror("Non valid element type for slave discretization!");
  }


  if (*boundary_ele==false || integrationtype==INPAR::MORTAR::inttype_fast)
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

      //**********************************************************************
      // loop over all mele
      //**********************************************************************
      for(int nummaster=0;nummaster<msize;++nummaster)
      {
        LINALG::SerialDenseMatrix mcoord(3,meles[nummaster]->NumNode());
        meles[nummaster]->GetNodalCoords(mcoord);

        // project Gauss point onto master element
        MORTAR::MortarProjector projector(3);
        projector.ProjectGaussPoint3D(sele,sxi,*meles[nummaster],mxi,projalpha);

        bool is_on_mele=true;

        // check GP projection
        double tol = 0.00;
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

        // evaluate Lagrange multiplier shape functions (on slave element)
        sele.EvaluateShapeLagMult(shapefcn_,sxi,lmval,lmderiv,nrow);

        // evaluate trace space shape functions (on both elements)
        sele.EvaluateShape(sxi,sval,sderiv,nrow);
        meles[nummaster]->EvaluateShape(mxi,mval,mderiv,nmnode);

        // evaluate the two Jacobians (int. cell and slave element)
        double jacslave = sele.Jacobian(sxi);

        if (is_on_mele==true)
        {
          *proj_=true;
          iter_proj+=1;
          // compute cell D/M matrix *******************************************
          // standard shape functions
          if (shapefcn_ == INPAR::MORTAR::shape_standard)
          {
            // loop over all mseg matrix entries
            // !!! nrow represents the slave Lagrange multipliers !!!
            // !!! ncol represents the master dofs                !!!
            // (this DOES matter here for mseg, as it might
            // sometimes be rectangular, not quadratic!)
            for (int j=0; j<nrow*ndof; ++j)
            {
              // integrate LinM
              for (int k=0; k<nmnode*ndof; ++k)
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/ndof);

                // multiply the two shape functions
                double prod = lmval[jindex]*mval[kindex];

                // isolate the mseg entries to be filled and
                // add current Gauss point's contribution to mseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                  (*mseg)(j, k+nummaster*nmnode*ndof) += prod*jacslave*wgt;
              }

              // integrate LinD
              for (int k=0; k<nrow*ndof; ++k)
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/ndof);

                // multiply the two shape functions
                double prod = lmval[jindex]*sval[kindex];

                // isolate the mseg entries to be filled and
                // add current Gauss point's contribution to mseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                    (*dseg)(j, k) += prod*jacslave*wgt;
              }
            }
          }

          // dual shape functions
          else if (shapefcn_ == INPAR::MORTAR::shape_dual)
          {
            // loop over all mseg matrix entries
            // !!! nrow represents the slave Lagrange multipliers !!!
            // !!! ncol represents the master dofs                !!!
            // (this DOES matter here for mseg, as it might
            // sometimes be rectangular, not quadratic!)
            for (int j=0; j<nrow*ndof; ++j)
            {
              // for dual shape functions we can make use
              // of the row summing lemma: D_jj = Sum(k) M_jk
              // hence, they can be combined into one single loop

              // integrate LinM and LinD
              for (int k=0; k<nmnode*ndof; ++k)
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/ndof);

                // multiply the two shape functions
                double prod = lmval[jindex]*mval[kindex];

                // isolate the mseg entries to be filled and
                // add current Gauss point's contribution to mseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                {
                  (*mseg)(j, k+nummaster*nmnode*ndof)+= prod*jacslave*wgt;
                  (*dseg)(j, j) += prod*jacslave*wgt;
                }
              }
            }
          }
          // compute cell D/M matrix *******************************************

          // evaluate 2nd deriv of trace space shape functions (on slave element)
          sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);

          // build interpolation of slave GP normal and coordinates
          double gpn[3] = {0.0,0.0,0.0};
          double sgpx[3] = {0.0, 0.0, 0.0};
          for (int i=0;i<nrow;++i)
          {
            MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*> (mynodes[i]);
            gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
            gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
            gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

            sgpx[0]+=sval[i]*scoord(0,i);
            sgpx[1]+=sval[i]*scoord(1,i);
            sgpx[2]+=sval[i]*scoord(2,i);
          }

          // normalize interpolated GP normal back to length 1.0 !!!
          double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
          if (length<1.0e-12) dserror("ERROR: IntegrateDerivCell3D: Divide by zero!");

          for (int i=0;i<3;++i)
            gpn[i]/=length;

          // build interpolation of master GP coordinates
          double mgpx[3] = {0.0, 0.0, 0.0};
          for (int i=0;i<nmnode;++i)
          {
            mgpx[0]+=mval[i]*mcoord(0,i);
            mgpx[1]+=mval[i]*mcoord(1,i);
            mgpx[2]+=mval[i]*mcoord(2,i);
          }

          // build normal gap at current GP
          double gap = 0.0;
          for (int i=0;i<3;++i)
            gap+=(mgpx[i]-sgpx[i])*gpn[i];

          // evaluate linearizations *******************************************
          // evaluate the derivative djacdxi = (Jac,xi , Jac,eta)
          double djacdxi[2] = {0.0, 0.0};
          static_cast<CONTACT::CoElement&>(sele).DJacDXi(djacdxi,sxi,ssecderiv);

          /* finite difference check of DJacDXi
          double inc = 1.0e-8;
          double jacnp1 = 0.0;
          double fdres[2] = {0.0, 0.0};
          sxi[0] += inc;
          jacnp1 = sele.Jacobian(sxi);
          fdres[0] = (jacnp1-jacslave)/inc;
          sxi[0] -= inc;
          sxi[1] += inc;
          jacnp1 = sele.Jacobian(sxi);
          fdres[1] = (jacnp1-jacslave)/inc;
          sxi[1] -= inc;
          std::cout << "DJacDXi: " << scientific << djacdxi[0] << " " << djacdxi[1] << endl;
          std::cout << "FD-DJacDXi: " << scientific << fdres[0] << " " << fdres[1] << endl << endl;*/

          // evaluate the slave Jacobian derivative
          std::map<int,double> jacslavemap;
          sele.DerivJacobian(sxi,jacslavemap);

          // evaluate the GP slave coordinate derivatives
          std::vector<std::map<int,double> > dsxigp(2);

          // evalute the GP master coordinate derivatives
          std::vector<std::map<int,double> > dmxigp(2);
          DerivXiGP3D(sele,*meles[nummaster],sxi,mxi,dsxigp,dmxigp,projalpha);

          // evaluate the GP gap function derivatives
          std::map<int,double> dgapgp;

          // we need the participating slave and master nodes
          DRT::Node** snodes = sele.Nodes();
          DRT::Node** mnodes = meles[nummaster]->Nodes();
          std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
          std::vector<MORTAR::MortarNode*> mmrtrnodes(meles[nummaster]->NumNode());

          for (int i=0;i<nrow;++i)
          {
            smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
            if (!smrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
          }

          for (int i=0;i<nmnode;++i)
          {
            mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
            if (!mmrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
          }

          // build directional derivative of slave GP normal (non-unit)
          std::map<int,double> dmap_nxsl_gp;
          std::map<int,double> dmap_nysl_gp;
          std::map<int,double> dmap_nzsl_gp;

          for (int i=0;i<nrow;++i)
          {
            std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
            std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
            std::map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

            for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
              dmap_nxsl_gp[p->first] += sval[i]*(p->second);
            for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
              dmap_nysl_gp[p->first] += sval[i]*(p->second);
            for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
              dmap_nzsl_gp[p->first] += sval[i]*(p->second);
          }

          // build directional derivative of slave GP normal (unit)
          std::map<int,double> dmap_nxsl_gp_unit;
          std::map<int,double> dmap_nysl_gp_unit;
          std::map<int,double> dmap_nzsl_gp_unit;

          double ll = length*length;
          double sxsx = gpn[0]*gpn[0]*ll;
          double sxsy = gpn[0]*gpn[1]*ll;
          double sxsz = gpn[0]*gpn[2]*ll;
          double sysy = gpn[1]*gpn[1]*ll;
          double sysz = gpn[1]*gpn[2]*ll;
          double szsz = gpn[2]*gpn[2]*ll;

          for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
          {
            dmap_nxsl_gp_unit[p->first] += 1/length*(p->second);
            dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsx*(p->second);
            dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
            dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
          }

          for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
          {
            dmap_nysl_gp_unit[p->first] += 1/length*(p->second);
            dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysy*(p->second);
            dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
            dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
          }

          for (CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
          {
            dmap_nzsl_gp_unit[p->first] += 1/length*(p->second);
            dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*szsz*(p->second);
            dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
            dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
          }

          // add everything to dgapgp
          for (CI p=dmap_nxsl_gp_unit.begin();p!=dmap_nxsl_gp_unit.end();++p)
            dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

          for (CI p=dmap_nysl_gp_unit.begin();p!=dmap_nysl_gp_unit.end();++p)
            dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

          for (CI p=dmap_nzsl_gp_unit.begin();p!=dmap_nzsl_gp_unit.end();++p)
            dgapgp[p->first] += (mgpx[2]-sgpx[2]) *(p->second);

          for (int z=0;z<nrow;++z)
          {
            for (int k=0;k<3;++k)
            {
              dgapgp[smrtrnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];
            }
          }

          for (int z=0;z<nmnode;++z)
          {
            for (int k=0;k<3;++k)
            {
              dgapgp[mmrtrnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

              for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
                dgapgp[p->first] += gpn[k] * mderiv(z,0) * mmrtrnodes[z]->xspatial()[k] * (p->second);

              for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
                dgapgp[p->first] += gpn[k] * mderiv(z,1) * mmrtrnodes[z]->xspatial()[k] * (p->second);
            }
          }
          // evaluate linearizations *******************************************

          // compute cell gap vector *******************************************
          // loop over all gseg vector entries
          // nrow represents the slave side dofs !!!  */
          for (int j=0;j<nrow;++j)
          {
            double prod = lmval[j]*gap;

            // add current Gauss point's contribution to gseg
            (*gseg)(j) += prod*jacslave*wgt;
          }
          // compute cell gap vector *******************************************

          // compute cell D/M linearization ************************************
          // loop over all slave nodes
          for (int j=0; j<nrow; ++j)
          {
            MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
            if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

            int sgid = mymrtrnode->Id();

            // standard shape functions
            if (shapefcn_ == INPAR::MORTAR::shape_standard)
            {
              // integrate LinM
              for (int k=0; k<nmnode; ++k)
              {
                // global master node ID
                int mgid = meles[nummaster]->Nodes()[k]->Id();
                double fac = 0.0;

                // get the correct map as a reference
                std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

                // (1) Lin(Phi) - dual shape functions --> 0
                // this vanishes here since there are no deformation-dependent dual functions

                // (2) Lin(NSlave) - slave GP coordinates --> 0

                // (3) Lin(NMaster) - master GP coordinates
                fac = wgt*lmval[j]*mderiv(k, 0)*jacslave;
                for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
                  dmmap_jk[p->first] += fac*(p->second);

                fac = wgt*lmval[j]*mderiv(k, 1)*jacslave;
                for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
                  dmmap_jk[p->first] += fac*(p->second);

                // (4) Lin(dsxideta) - intcell GP Jacobian --> 0

                // (5) Lin(dxdsxi) - slave GP Jacobian
                fac = wgt*lmval[j]*mval[k];
                for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
                  dmmap_jk[p->first] += fac*(p->second);

                // (6) Lin(dxdsxi) - slave GP coordinates --> 0
              } // loop over master nodes

              // integrate LinD
              for (int k=0; k<nrow; ++k)
              {
                // global master node ID -- evtl anderes sgid!!!
                int ssgid = sele.Nodes()[k]->Id();//int sgid = meles[nummaster]->Nodes()[k]->Id();
                double fac = 0.0;

                // get the correct map as a reference
                std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[ssgid];

                // (1) Lin(Phi) - dual shape functions  --> 0
                // this vanishes here since there are no deformation-dependent dual functions

                // (2) Lin(NSlave) - slave GP coordinates --> 0

                // (3) Lin(NSlave) - slave GP coordinates --> 0

                // (4) Lin(dsxideta) - intcell GP Jacobian --> 0

                // (5) Lin(dxdsxi) - slave GP Jacobian
                fac = wgt*lmval[j]*sval[k];
                for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
                  ddmap_jk[p->first] += fac*(p->second);

                // (6) Lin(dxdsxi) - slave GP coordinates --> 0
              } // loop over slave nodes
            }

            //************************************
            // dual shape functions
            //************************************
            else if (shapefcn_ == INPAR::MORTAR::shape_dual)
            {
              // get the D-map as a reference
              std::map<int,double>& ddmap_jj = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

              // integrate LinM and LinD
              for (int k=0; k<nmnode; ++k)
              {
                // global master node ID
                int mgid = meles[nummaster]->Nodes()[k]->Id();
                double fac = 0.0;

                // get the correct map as a reference
                std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

                // (1) Lin(Phi) - dual shape functions
                if (duallin)
                  for (int m=0; m<nrow; ++m)
                  {
                    fac = wgt*sval[m]*mval[k]*jacslave;
                    for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
                    {
                      dmmap_jk[p->first] += fac*(p->second);
                      ddmap_jj[p->first] += fac*(p->second);
                    }
                  }

                // (2) Lin(Phi) - slave GP coordinates --> 0

                // (3) Lin(NMaster) - master GP coordinates
                fac = wgt*lmval[j]*mderiv(k, 0)*jacslave;
                for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
                {
                  dmmap_jk[p->first] += fac*(p->second);
                  ddmap_jj[p->first] += fac*(p->second);
                }
                fac = wgt*lmval[j]*mderiv(k, 1)*jacslave;
                for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
                {
                  dmmap_jk[p->first] += fac*(p->second);
                  ddmap_jj[p->first] += fac*(p->second);
                }

                // (4) Lin(dsxideta) - intcell GP Jacobian --> 0

                // (5) Lin(dxdsxi) - slave GP Jacobian
                fac = wgt*lmval[j]*mval[k];
                for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
                {
                  dmmap_jk[p->first] += fac*(p->second);
                  ddmap_jj[p->first] += fac*(p->second);
                }

                // (6) Lin(dxdsxi) - slave GP coordinates --> 0
              } // loop over master nodes
            } // shapefcn_ switch
          } // loop over slave nodes
          // compute cell D/M linearization ************************************

          // compute cell gap linearization ************************************
          for (int j=0;j<nrow;++j)
          {
            MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
            if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

            double fac = 0.0;

            // get the corresponding map as a reference
            std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

            // standard shape functions
            if( shapefcn_ == INPAR::MORTAR::shape_standard )
            {
              // (1) Lin(Phi) - dual shape functions --> 0
              // this vanishes here since there are no deformation-dependent dual functions

               // (2) Lin(NSlave) - slave GP coordinates --> 0

               // (3) Lin(g) - gap function
               fac = wgt*lmval[j]*jacslave;
               for (CI p=dgapgp.begin(); p!=dgapgp.end(); ++p)
                 dgmap[p->first] += fac*(p->second);

               // (4) Lin(dsxideta) - intcell GP Jacobian --> 0

               // (5) Lin(dxdsxi) - slave GP Jacobian
               fac = wgt*lmval[j]*gap;
               for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
                 dgmap[p->first] += fac*(p->second);

               // (6) Lin(dxdsxi) - slave GP coordinates --> 0
            }

            // dual shape functions
            else if( shapefcn_ == INPAR::MORTAR::shape_dual )
            {
              // (1) Lin(Phi) - dual shape functions
              if (duallin)
                for (int m=0; m<nrow; ++m)
                {
                  fac = wgt*sval[m]*gap*jacslave;
                  for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
                    dgmap[p->first] += fac*(p->second);
                }

              // (2) Lin(Phi) - slave GP coordinates --> 0

              // (3) Lin(g) - gap function
              fac = wgt*lmval[j]*jacslave;
              for (CI p=dgapgp.begin(); p!=dgapgp.end(); ++p)
                dgmap[p->first] += fac*(p->second);

              // (4) Lin(dsxideta) - intcell GP Jacobian --> 0

              // (5) Lin(dxdsxi) - slave GP Jacobian
              fac = wgt*lmval[j]*gap;
              for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
                dgmap[p->first] += fac*(p->second);

              // (6) Lin(dxdsxi) - slave GP coordinates --> 0
            }
          }
          // compute cell gap linearization ************************************
        }//is_on_mele
        if (is_on_mele==true) break;
      }//mele loop
      //if one gp has counterparts on 2 elements --> non-uniqueness
      if (iter_proj>1)
      {
        dserror("Multiple feasible projections of one integration point!");
      }
    }//GP-loop
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 02/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3D(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType dt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3D called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

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
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  mele.GetNodalCoords(mcoord);

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape()!=MORTAR::MortarElement::tri3))
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // note that the third component of sxi is necessary!
    // (although it will always be 0.0 of course)
    double tempsxi[3] = {0.0, 0.0, 0.0};
    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};
    double projalpha = 0.0;

    // get Gauss point in slave element coordinates
    cell->LocalToGlobal(eta,tempsxi,0);
    sxi[0] = tempsxi[0];
    sxi[1] = tempsxi[1];

    // project Gauss point onto master element
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPoint3D(sele,sxi,mele,mxi,projalpha);

    // check GP projection
    double tol = 0.01;
    if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
        std::cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
        std::cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(shapefcn_,sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the two Jacobians (int. cell and slave element)
    double jaccell = cell->Jacobian(eta);
    double jacslave = sele.Jacobian(sxi);

    // compute cell D/M matrix *******************************************
    // standard shape functions
    if (shapefcn_ == INPAR::MORTAR::shape_standard)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate LinM
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg  
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jaccell*jacslave*wgt;
        }
        
        // integrate LinD
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j, k) += prod*jaccell*jacslave*wgt;
        }
      }
    }
    
    // dual shape functions
    else if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(k) M_jk
        // hence, they can be combined into one single loop
        
        // integrate LinM and LinD
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg  
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jaccell*jacslave*wgt;
            (*dseg)(j, j) += prod*jaccell*jacslave*wgt;
          }
        }
      }
    }
    // compute cell D/M matrix *******************************************

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);

    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
      gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
      gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

      sgpx[0]+=sval[i]*scoord(0,i);
      sgpx[1]+=sval[i]*scoord(1,i);
      sgpx[2]+=sval[i]*scoord(2,i);
    }

    // normalize interpolated GP normal back to length 1.0 !!!
    double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
    if (length<1.0e-12) dserror("ERROR: IntegrateDerivCell3D: Divide by zero!");

    for (int i=0;i<3;++i)
      gpn[i]/=length;

    // build interpolation of master GP coordinates
    double mgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<ncol;++i)
    {
      mgpx[0]+=mval[i]*mcoord(0,i);
      mgpx[1]+=mval[i]*mcoord(1,i);
      mgpx[2]+=mval[i]*mcoord(2,i);
    }

    // build normal gap at current GP
    double gap = 0.0;
    for (int i=0;i<3;++i)
      gap+=(mgpx[i]-sgpx[i])*gpn[i];

    // evaluate linearizations *******************************************
    // evaluate the derivative djacdxi = (Jac,xi , Jac,eta)
    double djacdxi[2] = {0.0, 0.0};
    static_cast<CONTACT::CoElement&>(sele).DJacDXi(djacdxi,sxi,ssecderiv);

    /* finite difference check of DJacDXi
    double inc = 1.0e-8;
    double jacnp1 = 0.0;
    double fdres[2] = {0.0, 0.0};
    sxi[0] += inc;
    jacnp1 = sele.Jacobian(sxi);
    fdres[0] = (jacnp1-jacslave)/inc;
    sxi[0] -= inc;
    sxi[1] += inc;
    jacnp1 = sele.Jacobian(sxi);
    fdres[1] = (jacnp1-jacslave)/inc;
    sxi[1] -= inc;
    std::cout << "DJacDXi: " << std::scientific << djacdxi[0] << " " << djacdxi[1] << endl;
    std::cout << "FD-DJacDXi: " << std::scientific << fdres[0] << " " << fdres[1] << endl << endl;*/

    // evaluate the slave Jacobian derivative
    std::map<int,double> jacslavemap;
    sele.DerivJacobian(sxi,jacslavemap);

    // evaluate the intcell Jacobian derivative
    // these are pre-factors for intcell vertex coordinate linearizations
    std::vector<double> jacintcellvec(2*(cell->NumVertices()));
    cell->DerivJacobian(eta,jacintcellvec);

    // evalute the GP slave coordinate derivatives
    std::vector<std::map<int,double> > dsxigp(2);
    int nvcell = cell->NumVertices();
    LINALG::SerialDenseVector svalcell(nvcell);
    LINALG::SerialDenseMatrix sderivcell(nvcell,2,true);
    cell->EvaluateShape(eta,svalcell,sderivcell);

    for (int v=0;v<nvcell;++v)
    {
      for (CI p=(cell->GetDerivVertex(v))[0].begin();p!=(cell->GetDerivVertex(v))[0].end();++p)
        dsxigp[0][p->first] += svalcell[v] * (p->second);
      for (CI p=(cell->GetDerivVertex(v))[1].begin();p!=(cell->GetDerivVertex(v))[1].end();++p)
        dsxigp[1][p->first] += svalcell[v] * (p->second);
    }

    // evalute the GP master coordinate derivatives
    std::vector<std::map<int,double> > dmxigp(2);
    DerivXiGP3D(sele,mele,sxi,mxi,dsxigp,dmxigp,projalpha);

    // evaluate the GP gap function derivatives
    std::map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
      if (!smrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
      if (!mmrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    std::map<int,double> dmap_nxsl_gp;
    std::map<int,double> dmap_nysl_gp;
    std::map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
      std::map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }
    }

    // build directional derivative of slave GP normal (unit)
    std::map<int,double> dmap_nxsl_gp_unit;
    std::map<int,double> dmap_nysl_gp_unit;
    std::map<int,double> dmap_nzsl_gp_unit;

    double ll = length*length;
    double sxsx = gpn[0]*gpn[0]*ll;
    double sxsy = gpn[0]*gpn[1]*ll;
    double sxsz = gpn[0]*gpn[2]*ll;
    double sysy = gpn[1]*gpn[1]*ll;
    double sysz = gpn[1]*gpn[2]*ll;
    double szsz = gpn[2]*gpn[2]*ll;

    for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    {
      dmap_nxsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsx*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
    }

    for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
    {
      dmap_nysl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysy*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
    }

    for (CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
    {
      dmap_nzsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*szsz*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
    }

    // add everything to dgapgp
    for (CI p=dmap_nxsl_gp_unit.begin();p!=dmap_nxsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

    for (CI p=dmap_nysl_gp_unit.begin();p!=dmap_nysl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

    for (CI p=dmap_nzsl_gp_unit.begin();p!=dmap_nzsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[2]-sgpx[2]) *(p->second);

    for (int z=0;z<nrow;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[smrtrnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * smrtrnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,1) * smrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[mmrtrnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mmrtrnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,1) * mmrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************
    
    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
        prod = sval[j]*gap;
      // usual standard or dual LM approach
      else
        prod = lmval[j]*gap;

      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jaccell*jacslave*wgt;
    }
    // compute cell gap vector *******************************************

    // compute cell D/M linearization ************************************
    // loop over all slave nodes
    for (int j=0; j<nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

      int sgid = mymrtrnode->Id();

      // standard shape functions
      if (shapefcn_ == INPAR::MORTAR::shape_standard)
      {   
        // integrate LinM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;
  
          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];
  
          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions
  
          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmderiv(j, 1)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
  
          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval[j]*mderiv(k, 0)*jaccell*jacslave;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*mderiv(k, 1)*jaccell*jacslave;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
  
          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*mval[k]*jacslave;
          for (int m=0; m<(int)jacintcellvec.size(); ++m)
          {
            int v = m/2; // which vertex?
            int dof = m%2; // which dof?
            for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            {
              dmmap_jk[p->first] += fac * jacintcellvec[m] * (p->second);
            }
          }
  
          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*lmval[j]*mval[k]*jaccell;
          for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
  
          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*lmval[j]*mval[k]*jaccell*djacdxi[0];
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*mval[k]*jaccell*djacdxi[1];
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes
        
        // integrate LinD
        for (int k=0; k<nrow; ++k)
        {
          // global master node ID
          int sgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*sval[k]*jaccell*jacslave;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmderiv(j, 1)*sval[k]*jaccell*jacslave;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmval[j]*sderiv(k, 0)*jaccell*jacslave;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*sderiv(k, 1)*jaccell*jacslave;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*sval[k]*jacslave;
          for (int m=0; m<(int)jacintcellvec.size(); ++m)
          {
            int v = m/2; // which vertex?
            int dof = m%2; // which dof?
            for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            {
              ddmap_jk[p->first] += fac * jacintcellvec[m] * (p->second);
            }
          }

          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*lmval[j]*sval[k]*jaccell;
          for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*lmval[j]*sval[k]*jaccell*djacdxi[0];
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*sval[k]*jaccell*djacdxi[1];
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);
        } // loop over slave nodes
      }
      
      // dual shape functions
      else if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // get the D-map as a reference
        std::map<int,double>& ddmap_jj = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];
  
        // integrate LinM and LinD
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;
  
          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];
  
          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0; m<nrow; ++m)
            {
              fac = wgt*sval[m]*mval[k]*jaccell*jacslave;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                ddmap_jj[p->first] += fac*(p->second);
              }
            }
  
          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*lmderiv(j, 1)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
  
          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval[j]*mderiv(k, 0)*jaccell*jacslave;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*lmval[j]*mderiv(k, 1)*jaccell*jacslave;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
  
          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*mval[k]*jacslave;
          for (int m=0; m<(int)jacintcellvec.size(); ++m)
          {
            int v = m/2; // which vertex?
            int dof = m%2; // which dof?
            for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            {
              dmmap_jk[p->first] += fac * jacintcellvec[m] * (p->second);
              ddmap_jj[p->first] += fac * jacintcellvec[m] *(p->second);
            }
          }
  
          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*lmval[j]*mval[k]*jaccell;
          for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
  
          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*lmval[j]*mval[k]*jaccell*djacdxi[0];
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*lmval[j]*mval[k]*jaccell*djacdxi[1];
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    } // loop over slave nodes
    // compute cell D/M linearization ************************************

    // compute cell gap linearization ************************************
    for (int j=0;j<nrow;++j)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

      // standard shape functions
      if( shapefcn_ == INPAR::MORTAR::shape_standard )
      {
        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions
   
         // (2) Lin(NSlave) - slave GP coordinates
         fac = wgt*lmderiv(j, 0)*gap*jaccell*jacslave;
         for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
           dgmap[p->first] += fac*(p->second);
         fac = wgt*lmderiv(j, 1)*gap*jaccell*jacslave;
         for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
           dgmap[p->first] += fac*(p->second);
   
         // (3) Lin(g) - gap function
         fac = wgt*lmval[j]*jaccell*jacslave;
         for (CI p=dgapgp.begin(); p!=dgapgp.end(); ++p)
           dgmap[p->first] += fac*(p->second);
   
         // (4) Lin(dsxideta) - intcell GP Jacobian
         fac = wgt*lmval[j]*gap*jacslave;
         for (int m=0; m<(int)jacintcellvec.size(); ++m)
         {
           int v = m/2; // which vertex?
           int dof = m%2; // which dof?
           for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
             dgmap[p->first] += fac * jacintcellvec[m] * (p->second);
         }
   
         // (5) Lin(dxdsxi) - slave GP Jacobian
         fac = wgt*lmval[j]*gap*jaccell;
         for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
           dgmap[p->first] += fac*(p->second);
   
         // (6) Lin(dxdsxi) - slave GP coordinates
         fac = wgt*lmval[j]*gap*jaccell*djacdxi[0];
         for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
           dgmap[p->first] += fac*(p->second);
         fac = wgt*lmval[j]*gap*jaccell*djacdxi[1];
         for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
           dgmap[p->first] += fac*(p->second);        
      }
      
      // dual shape functions
      else if( shapefcn_ == INPAR::MORTAR::shape_dual)
      {
        // (1) Lin(Phi) - dual shape functions
        if (duallin)
          for (int m=0; m<nrow; ++m)
          {
            fac = wgt*sval[m]*gap*jaccell*jacslave;
            for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              dgmap[p->first] += fac*(p->second);
          }
  
        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmderiv(j, 0)*gap*jaccell*jacslave;
        for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          dgmap[p->first] += fac*(p->second);
        fac = wgt*lmderiv(j, 1)*gap*jaccell*jacslave;
        for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          dgmap[p->first] += fac*(p->second);
  
        // (3) Lin(g) - gap function
        fac = wgt*lmval[j]*jaccell*jacslave;
        for (CI p=dgapgp.begin(); p!=dgapgp.end(); ++p)
          dgmap[p->first] += fac*(p->second);
  
        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*lmval[j]*gap*jacslave;
        for (int m=0; m<(int)jacintcellvec.size(); ++m)
        {
          int v = m/2; // which vertex?
          int dof = m%2; // which dof?
          for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            dgmap[p->first] += fac * jacintcellvec[m] * (p->second);
        }
  
        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*lmval[j]*gap*jaccell;
        for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
          dgmap[p->first] += fac*(p->second);
  
        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*lmval[j]*gap*jaccell*djacdxi[0];
        for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          dgmap[p->first] += fac*(p->second);
        fac = wgt*lmval[j]*gap*jaccell*djacdxi[1];
        for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          dgmap[p->first] += fac*(p->second);
      }

      // Petrov-Galerkin approach for LM
      else if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // (1) Lin(Phi) - does not exist in gap for Petrov-Galerkin interpolation
        // as std shape functions are used here

        // (2) Lin(N) - slave GP coordinates
        fac = wgt*sderiv(j, 0)*gap*jaccell*jacslave;
        for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          dgmap[p->first] += fac*(p->second);
        fac = wgt*sderiv(j, 1)*gap*jaccell*jacslave;
        for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*sval[j]*jaccell*jacslave;
        for (CI p=dgapgp.begin(); p!=dgapgp.end(); ++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*sval[j]*gap*jacslave;
        for (int m=0; m<(int)jacintcellvec.size(); ++m)
        {
          int v = m/2; // which vertex?
          int dof = m%2; // which dof?
          for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            dgmap[p->first] += fac * jacintcellvec[m] * (p->second);
        }

        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*sval[j]*gap*jaccell;
        for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
          dgmap[p->first] += fac*(p->second);

        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*sval[j]*gap*jaccell*djacdxi[0];
        for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          dgmap[p->first] += fac*(p->second);
        fac = wgt*sval[j]*gap*jaccell*djacdxi[1];
        for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          dgmap[p->first] += fac*(p->second);
      }
    }
    // compute cell gap linearization ************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 |  This is the auxiliary plane coupling version!!!                     |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg,
     Teuchos::RCP<Epetra_SerialDenseVector> mdisssegs,
     Teuchos::RCP<Epetra_SerialDenseVector> mdisssegm,
     Teuchos::RCP<Epetra_SerialDenseMatrix> aseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> bseg,
     Teuchos::RCP<Epetra_SerialDenseVector> wseg)
{
  
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
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
  
  // flags for thermo-structure-interaction with contact
  bool tsi = false;
    if (DRT::Problem::Instance()->ProblemType()==prb_tsi) tsi=true;
  bool friction = false;     // friction 
  bool thermolagmult = true; // thermal contact with or without LM

  if(tsi)
  {
    const Teuchos::ParameterList& input = DRT::Problem::Instance()->ContactDynamicParams();

    if(DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_none)
      friction = true;
    if (DRT::INPUT::IntegralValue<int>(input,"THERMOLAGMULT")==false)
      thermolagmult = false;
  }
  
  // contact with wear
  bool wear = false;  
  if(DRT::Problem::Instance()->ContactDynamicParams().get<double>("WEARCOEFF")!= 0.0)
    wear = true;

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

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
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  mele.GetNodalCoords(mcoord);
  
  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;
  
  // get them in the case of tsi
  if ((tsi and friction) or wear)
  {
    scoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalCoordsOld(*scoordold);
    mcoordold = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,mele.NumNode()));
    mele.GetNodalCoordsOld(*mcoordold);
    lagmult = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalLagMult(*lagmult);
  }

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape()!=MORTAR::MortarElement::tri3))
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
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
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,sele,sxi,sprojalpha);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,mele,mxi,mprojalpha);
    
    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt==DRT::Element::quad4 || sdt==DRT::Element::quad8 || sdt==DRT::Element::quad9)
    {
      if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(shapefcn_,sxi,lmval,lmderiv,nrow);
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // compute cell D/M matrix *******************************************
    // standard shape functions
    if (shapefcn_ == INPAR::MORTAR::shape_standard)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate LinM
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jac*wgt;
        }

        // integrate LinD
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j, k) += prod*jac*wgt;
        }
      }
    }
    
    // dual shape functions
    else if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(k) M_jk
        // hence, they can be combined into one single loop

        // integrate LinM and LinD
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jac*wgt;
            (*dseg)(j, j) += prod*jac*wgt;
          }
        }
      } // nrow*ndof loop
    } // shapefcn_ switch
    // compute cell D/M matrix *******************************************

    
    // evaluate matrix A for thermo-structure-interaction
    // in the case of using thermal lagrange multipliers
    if ((tsi and friction) and thermolagmult == true)
    {
      // loop over all aseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the dofs                       !!!
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate LinM and LinD
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*lmval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to aseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*aseg)(j, k) += prod*jac*wgt;
          }
        }
      } // nrow*ndof loop
    }

    // evaluate matrix B for thermo-structure-interaction
    // in the case of NOT using thermal lagrange multipliers
    if ((tsi and friction) and thermolagmult == false)
    {
      // loop over all bseg matrix entries
      // !!! nrow represents the master shape functions     !!!
      // !!! ncol represents the dofs                       !!!
      for (int j=0; j<ncol*ndof; ++j)
      {
        // integrate bseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = mval[jindex]*mval[kindex];

          // isolate the bseg entries to be filled and
          // add current Gauss point's contribution to bseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*bseg)(j, k) += prod*jac*wgt;
          }
        }
      } // nrow*ndof loop
    } // tsi
    
    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);
    
    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
      gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
      gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

      sgpx[0]+=sval[i]*scoord(0,i);
      sgpx[1]+=sval[i]*scoord(1,i);
      sgpx[2]+=sval[i]*scoord(2,i);
    }

    // normalize interpolated GP normal back to length 1.0 !!!
    double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
    if (length<1.0e-12) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Divide by zero!");

    for (int i=0;i<3;++i)
      gpn[i]/=length;

    // build interpolation of master GP coordinates
    double mgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<ncol;++i)
    {
      mgpx[0]+=mval[i]*mcoord(0,i);
      mgpx[1]+=mval[i]*mcoord(1,i);
      mgpx[2]+=mval[i]*mcoord(2,i);
    }

    // build normal gap at current GP
    double gap = 0.0;
    for (int i=0;i<3;++i)
      gap+=(mgpx[i]-sgpx[i])*gpn[i];
    
    
    // evaluate mechanical dissipation in the case of thermo-structure-
    // -interaction with contact

    // mechanical dissipation and wear
    double mechdiss = 0;
    double wearval = 0.0;
    
    if((tsi and friction) or wear)
    {
      // tangent plane  
      Epetra_SerialDenseMatrix tanplane(3,3);
      tanplane(0,0)= 1-(gpn[0]*gpn[0]);
      tanplane(0,1)=  -(gpn[0]*gpn[1]);
      tanplane(0,2)=  -(gpn[0]*gpn[2]);
      tanplane(1,0)=  -(gpn[1]*gpn[0]);
      tanplane(1,1)= 1-(gpn[1]*gpn[1]);
      tanplane(1,2)=  -(gpn[1]*gpn[2]);
      tanplane(2,0)=  -(gpn[2]*gpn[0]);
      tanplane(2,1)=  -(gpn[2]*gpn[1]);
      tanplane(2,2)= 1-(gpn[2]*gpn[2]);
      
      // interpolation of slave GP jumps (relative displacement increment)
      double sgpjump[3] = {0.0, 0.0, 0.0};
      for (int i=0;i<nrow;++i)
      {
        sgpjump[0]+=sval[i]*(scoord(0,i)-(*scoordold)(0,i));
        sgpjump[1]+=sval[i]*(scoord(1,i)-(*scoordold)(1,i));
        sgpjump[2]+=sval[i]*(scoord(2,i)-(*scoordold)(2,i));
      }
      
      // interpolation of master GP jumps (relative displacement increment)
      double mgpjump[3] = {0.0, 0.0, 0.0};
      for (int i=0;i<ncol;++i)
      {
        mgpjump[0]+=mval[i]*(mcoord(0,i)-(*mcoordold)(0,i));
        mgpjump[1]+=mval[i]*(mcoord(1,i)-(*mcoordold)(1,i));
        mgpjump[2]+=mval[i]*(mcoord(2,i)-(*mcoordold)(2,i));
      }
   
      // build jump (relative displacement increment) at current GP
      Epetra_SerialDenseMatrix jump (3,1);
      jump(0,0)=(sgpjump[0]-mgpjump[0]); 
      jump(1,0)=(sgpjump[1]-mgpjump[1]); 
      jump(2,0)=(sgpjump[2]-mgpjump[2]); 
      
      // build lagrange multiplier at current GP
      Epetra_SerialDenseMatrix lm (3,1);
      for (int i=0;i<nrow;++i)
      {
        lm(0,0)+=lmval[i]*(*lagmult)(0,i);
        lm(1,0)+=lmval[i]*(*lagmult)(1,i);
        lm(2,0)+=lmval[i]*(*lagmult)(2,i);
      }
       
      // build tangential jump
      Epetra_SerialDenseMatrix jumptan(3,1);
      jumptan.Multiply('N','N',1,tanplane,jump,0.0);
  
      // build tangential lm
      Epetra_SerialDenseMatrix lmtan(3,1);
      lmtan.Multiply('N','N',1,tanplane,lm,0.0);

      // evaluate mechanical dissipation
      if (tsi and friction)
      {  
        for (int i=0;i<3;i++)
          mechdiss += lmtan(i,0)*jumptan(i,0);
      }

      // evaluate wear
      // not including wearcoefficient
      if (wear)
      { 
        // normal contact stress
        for (int i=0;i<3;++i)
          wearval+=lm(i,0)*gpn[i];
        
        // absolute value of relative tangential jump
        double jumpval = 0.0;
        jumpval = sqrt(jumptan(0,0)*jumptan(0,0)+jumptan(1,0)*jumptan(1,0)+jumptan(2,0)*jumptan(2,0));
        
        // product
        wearval = abs(wearval)*jumpval;
      }
    }  
    
    // evaluate linearizations *******************************************
    // evaluate the intcell Jacobian derivative
    std::map<int,double> jacintcellmap;
    cell->DerivJacobian(eta,jacintcellmap);

    // evaluate global GP coordinate derivative
    std::vector<std::map<int,double> > lingp(3);
    int nvcell = cell->NumVertices();
    LINALG::SerialDenseVector svalcell(nvcell);
    LINALG::SerialDenseMatrix sderivcell(nvcell,2,true);
    cell->EvaluateShape(eta,svalcell,sderivcell);

    for (int v=0;v<nvcell;++v)
    {
      for (CI p=(cell->GetDerivVertex(v))[0].begin();p!=(cell->GetDerivVertex(v))[0].end();++p)
        lingp[0][p->first] += svalcell[v] * (p->second);
      for (CI p=(cell->GetDerivVertex(v))[1].begin();p!=(cell->GetDerivVertex(v))[1].end();++p)
        lingp[1][p->first] += svalcell[v] * (p->second);
      for (CI p=(cell->GetDerivVertex(v))[2].begin();p!=(cell->GetDerivVertex(v))[2].end();++p)
        lingp[2][p->first] += svalcell[v] * (p->second);
    }

    // evalute the GP slave coordinate derivatives
    std::vector<std::map<int,double> > dsxigp(2);
    DerivXiGP3DAuxPlane(sele,sxi,cell->Auxn(),dsxigp,sprojalpha,cell->GetDerivAuxn(),lingp);

    // evalute the GP master coordinate derivatives
    std::vector<std::map<int,double> > dmxigp(2);
    DerivXiGP3DAuxPlane(mele,mxi,cell->Auxn(),dmxigp,mprojalpha,cell->GetDerivAuxn(),lingp);

    // evaluate the GP gap function derivatives
    std::map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
      if (!smrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
      if (!mmrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    std::map<int,double> dmap_nxsl_gp;
    std::map<int,double> dmap_nysl_gp;
    std::map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
      std::map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }
    }

    // build directional derivative of slave GP normal (unit)
    std::map<int,double> dmap_nxsl_gp_unit;
    std::map<int,double> dmap_nysl_gp_unit;
    std::map<int,double> dmap_nzsl_gp_unit;

    double ll = length*length;
    double sxsx = gpn[0]*gpn[0]*ll;
    double sxsy = gpn[0]*gpn[1]*ll;
    double sxsz = gpn[0]*gpn[2]*ll;
    double sysy = gpn[1]*gpn[1]*ll;
    double sysz = gpn[1]*gpn[2]*ll;
    double szsz = gpn[2]*gpn[2]*ll;

    for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    {
      dmap_nxsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsx*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
    }

    for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
    {
      dmap_nysl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysy*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
    }

    for (CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
    {
      dmap_nzsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*szsz*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
    }

    // add everything to dgapgp
    for (CI p=dmap_nxsl_gp_unit.begin();p!=dmap_nxsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

    for (CI p=dmap_nysl_gp_unit.begin();p!=dmap_nysl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

    for (CI p=dmap_nzsl_gp_unit.begin();p!=dmap_nzsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[2]-sgpx[2]) *(p->second);

    for (int z=0;z<nrow;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[smrtrnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * smrtrnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,1) * smrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[mmrtrnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mmrtrnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,1) * mmrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************
    
    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = 0.0;
      // Petrov-Galerkin approach (dual LM for D/M but standard LM for gap)
      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
        prod = sval[j]*gap;
      // usual standard or dual LM approach
      else
        prod = lmval[j]*gap;
      
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jac*wgt;
    }
    // compute cell gap vector *******************************************
    
    // tsi with contact
    if(tsi and friction)
    {
      // compute cell mechanical dissipation / slave *********************
      // loop over all mdisssegs vector entries
      // nrow represents the slave side dofs !!!
      for (int j=0;j<nrow;++j)
      {
        double prod;
        if(thermolagmult==true) prod = lmval[j]*mechdiss;
        else                    prod =  sval[j]*mechdiss;
      
        // add current Gauss point's contribution to mdisssegs
        (*mdisssegs)(j) += prod*jac*wgt;
      }
      // compute cell mechanical dissipation vector / slave **************

      // compute cell jump vector / master *******************************
      // loop over all mdisssegm vector entries
      // ncol represents the master side dofs !!!
      for (int j=0;j<ncol;++j)
      {
        double prod = mval[j]*mechdiss;
      
        // add current Gauss point's contribution to mdisssegm
        (*mdisssegm)(j) += prod*jac*wgt;
      }
      // compute cell mechanical dissipation vector / master *************
    } // tsi with contact
    
    // compute cell wear vector ******************************************
    // loop over all wseg vector entries
    // nrow represents the slave side dofs !!!
    if (wear)
    {
      for (int j=0;j<nrow;++j)
      {
        double prod = lmval[j]*wearval;

        // add current Gauss point's contribution to gseg
        (*wseg)(j) += prod*jac*wgt;
      }
    }
    // compute cell wear vector ******************************************
      
    // compute cell D/M linearization ************************************
    // loop over slave nodes
    for (int j=0; j<nrow; ++j)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      int sgid = mymrtrnode->Id();

      // standard shape functions
      if (shapefcn_ == INPAR::MORTAR::shape_standard)
      {
        // integrate LinM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmderiv(j, 1)*mval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval[j]*mderiv(k, 0)*jac;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes

        // integrate LinD
        for (int k=0; k<nrow; ++k)
        {
          // global slave node ID
          int sgid = sele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*sval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmderiv(j, 1)*sval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmval[j]*sderiv(k, 0)*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*sderiv(k, 1)*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*sval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);
        } // loop over slave nodes
      }
      
      // dual shape functions
      else if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // get the D-map as a reference
        std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

        // integrate LinM and LinD
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0; m<nrow; ++m)
            {
              fac = wgt*sval[m]*mval[k]*jac;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                ddmap_jk[p->first] += fac*(p->second);
              }
            }

          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jk[p->first] += fac*(p->second);
          }
          fac = wgt*lmderiv(j, 1)*mval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jk[p->first] += fac*(p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval[j]*mderiv(k, 0)*jac;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jk[p->first] += fac*(p->second);
          }
          fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jk[p->first] += fac*(p->second);
          }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jk[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    } // loop over slave nodes
    // compute cell D/M linearization ************************************

    // compute cell gap linearization ************************************
    for (int j=0;j<nrow;++j)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

      // switch if Petrov-Galerkin approach for LM is applied
      if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
      {
        // (1) Lin(Phi) - does not exist here for Petrov-Galerkin approach

        // (2) Lin(N) - slave GP coordinates
        fac = wgt*sderiv(j,0)*gap*jac;
        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgmap[p->first] += fac*(p->second);

        fac = wgt*sderiv(j,1)*gap*jac;
        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*sval[j]*jac;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*sval[j]*gap;
        for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
          dgmap[p->first] += fac*(p->second);
      }
      
      // the usual standard or dual LM approach
      else
      {
        // (1) Lin(Phi) - dual shape functions
        if (duallin)
          for (int m=0;m<nrow;++m)
          {
            fac = wgt*sval[m]*gap*jac;
            for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
              dgmap[p->first] += fac*(p->second);
          }

        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmderiv(j,0)*gap*jac;
        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgmap[p->first] += fac*(p->second);

        fac = wgt*lmderiv(j,1)*gap*jac;
        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*lmval[j]*jac;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*lmval[j]*gap;
        for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
          dgmap[p->first] += fac*(p->second);
      }
    }
    // compute cell gap linearization ************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 |  This is the QUADRATIC auxiliary plane coupling version!!!           |
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivCell3DAuxPlaneQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     INPAR::MORTAR::LagMultQuad lmtype,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
  
  // Petrov-Galerkin approach for LM not yet implemented
  if (shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("ERROR: Petrov-Galerkin approach not yet implemented for quadratic FE interpolation");

  /*std::cout << endl;
  std::cout << "Slave type: " << sele.Shape() << endl;
  std::cout << "SlaveElement Nodes:";
  for (int k=0;k<sele.NumNode();++k) std::cout << " " << sele.NodeIds()[k];
  std::cout << "\nMaster type: " << mele.Shape() << endl;
  std::cout << "MasterElement Nodes:";
  for (int k=0;k<mele.NumNode();++k) std::cout << " " << mele.NodeIds()[k];
  std::cout << endl;
  std::cout << "SlaveSub type: " << sintele.Shape() << endl;
  std::cout << "SlaveSubElement Nodes:";
  for (int k=0;k<sintele.NumNode();++k) std::cout << " " << sintele.NodeIds()[k];
  std::cout << "\nMasterSub type: " << mintele.Shape() << endl;
  std::cout << "MasterSubElement Nodes:";
  for (int k=0;k<mintele.NumNode();++k) std::cout << " " << mintele.NodeIds()[k];  */

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of slave and master IntElement
  DRT::Element::DiscretizationType sdt = sintele.Shape();
  DRT::Element::DiscretizationType mdt = mintele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();
  int nintrow = sintele.NumNode();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector svalmod(nrow);
  LINALG::SerialDenseMatrix sderivmod(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);
  LINALG::SerialDenseVector lmintval(nintrow);
  LINALG::SerialDenseMatrix lmintderiv(nintrow,2,true);

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);
  
  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  mele.GetNodalCoords(mcoord);

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
  DRT::Node** myintnodes = sintele.Nodes();
  if(!myintnodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
      (sele.Shape()!=MORTAR::MortarElement::tri3))
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }
  
  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic Lagrange multipliers on quad8 and tri6 elements
  bool dualquad3d = false;
  if ( (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
       (lmtype == INPAR::MORTAR::lagmult_quad_quad) &&
       (sele.Shape() == DRT::Element::quad8 || sele.Shape() == DRT::Element::tri6) )
  {
    dualquad3d = true;
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

    // project Gauss point onto slave integration element
    // project Gauss point onto master integration element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,sintele,sxi,sprojalpha);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,mintele,mxi,mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt==DRT::Element::quad4 || sdt==DRT::Element::quad8 || sdt==DRT::Element::quad9)
    {
      if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }

    // map Gauss point back to slave element (affine map)
    // map Gauss point back to master element (affine map)
    double psxi[2] = {0.0, 0.0};
    double pmxi[2] = {0.0, 0.0};
    sintele.MapToParent(sxi,psxi);
    mintele.MapToParent(mxi,pmxi);

    //std::cout << "SInt-GP:    " << sxi[0] << " " << sxi[1] << endl;
    //std::cout << "MInt-GP:    " << mxi[0] << " " << mxi[1] << endl;
    //std::cout << "SParent-GP: " << psxi[0] << " " << psxi[1] << endl;
    //std::cout << "MParent-GP: " << pmxi[0] << " " << pmxi[1] << endl;

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (bound)
      sele.EvaluateShapeLagMultLin(shapefcn_,psxi,lmval,lmderiv,nrow);
    else
    {
      sele.EvaluateShapeLagMult(shapefcn_,psxi,lmval,lmderiv,nrow);
      sintele.EvaluateShapeLagMult(shapefcn_,sxi,lmintval,lmintderiv,nintrow);
    }
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(psxi,sval,sderiv,nrow);
    if (dualquad3d) sele.EvaluateShape(psxi,svalmod,sderivmod,nrow,true);
    mele.EvaluateShape(pmxi,mval,mderiv,ncol);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // compute cell D/M matrix *******************************************
    // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
    if (shapefcn_ == INPAR::MORTAR::shape_standard &&
        (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin))
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate LinM
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);
  
          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }
  
        // integrate LinD
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }
    
    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    else if (shapefcn_ == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nintrow*ndof; ++j)
      {
        // integrate LinM
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);
  
          // multiply the two shape functions
          double prod = lmintval[jindex]*mval[lindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }
  
        // integrate LinD
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmintval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }
    
    // CASE 4: Dual LM shape functions and quadratic interpolation
    else if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(l) M_jl
        // hence, they can be combined into one single loop

        // integrate LinM and LinD
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];

          // isolate the mseg/dseg entries to be filled and
          // add current Gauss point's contribution to mseg/dseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
          {
            (*mseg)(j,l) += prod*jac*wgt;
            (*dseg)(j,j) += prod*jac*wgt;
          }
        }
      }
    }
    
    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic contact!");
    }
    // compute cell D/M matrix *******************************************

    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.Evaluate2ndDerivShape(psxi,ssecderiv,nrow);
    
    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
      gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
      gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

      sgpx[0]+=sval[i]*scoord(0,i);
      sgpx[1]+=sval[i]*scoord(1,i);
      sgpx[2]+=sval[i]*scoord(2,i);
    }

    // normalize interpolated GP normal back to length 1.0 !!!
    double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
    if (length<1.0e-12) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Divide by zero!");

    for (int i=0;i<3;++i)
      gpn[i]/=length;

    // build interpolation of master GP coordinates
    double mgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<ncol;++i)
    {
      mgpx[0]+=mval[i]*mcoord(0,i);
      mgpx[1]+=mval[i]*mcoord(1,i);
      mgpx[2]+=mval[i]*mcoord(2,i);
    }

    // build normal gap at current GP
    double gap = 0.0;
    for (int i=0;i<3;++i)
      gap+=(mgpx[i]-sgpx[i])*gpn[i];

    // evaluate linearizations *******************************************
    // evaluate the intcell Jacobian derivative
    std::map<int,double> jacintcellmap;
    cell->DerivJacobian(eta,jacintcellmap);

    // evaluate global GP coordinate derivative
    std::vector<std::map<int,double> > lingp(3);
    int nvcell = cell->NumVertices();
    LINALG::SerialDenseVector svalcell(nvcell);
    LINALG::SerialDenseMatrix sderivcell(nvcell,2,true);
    cell->EvaluateShape(eta,svalcell,sderivcell);

    for (int v=0;v<nvcell;++v)
    {
      for (CI p=(cell->GetDerivVertex(v))[0].begin();p!=(cell->GetDerivVertex(v))[0].end();++p)
        lingp[0][p->first] += svalcell[v] * (p->second);
      for (CI p=(cell->GetDerivVertex(v))[1].begin();p!=(cell->GetDerivVertex(v))[1].end();++p)
        lingp[1][p->first] += svalcell[v] * (p->second);
      for (CI p=(cell->GetDerivVertex(v))[2].begin();p!=(cell->GetDerivVertex(v))[2].end();++p)
        lingp[2][p->first] += svalcell[v] * (p->second);
    }

    // evalute the GP slave coordinate derivatives
    std::vector<std::map<int,double> > dsxigp(2);
    DerivXiGP3DAuxPlane(sintele,sxi,cell->Auxn(),dsxigp,sprojalpha,cell->GetDerivAuxn(),lingp);

    // evalute the GP master coordinate derivatives
    std::vector<std::map<int,double> > dmxigp(2);
    DerivXiGP3DAuxPlane(mintele,mxi,cell->Auxn(),dmxigp,mprojalpha,cell->GetDerivAuxn(),lingp);

    // map GP coordinate derivatives back to slave element (affine map)
    // map GP coordinate derivatives back to master element (affine map)
    std::vector<std::map<int,double> > dpsxigp(2);
    std::vector<std::map<int,double> > dpmxigp(2);
    sintele.MapToParent(dsxigp,dpsxigp);
    mintele.MapToParent(dmxigp,dpmxigp);

    // evaluate the GP gap function derivatives
    std::map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
      if (!smrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
      if (!mmrtrnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    std::map<int,double> dmap_nxsl_gp;
    std::map<int,double> dmap_nysl_gp;
    std::map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
      std::map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*smrtrnodes[i]->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*smrtrnodes[i]->MoData().n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }
    }

    // build directional derivative of slave GP normal (unit)
    std::map<int,double> dmap_nxsl_gp_unit;
    std::map<int,double> dmap_nysl_gp_unit;
    std::map<int,double> dmap_nzsl_gp_unit;

    double ll = length*length;
    double sxsx = gpn[0]*gpn[0]*ll;
    double sxsy = gpn[0]*gpn[1]*ll;
    double sxsz = gpn[0]*gpn[2]*ll;
    double sysy = gpn[1]*gpn[1]*ll;
    double sysz = gpn[1]*gpn[2]*ll;
    double szsz = gpn[2]*gpn[2]*ll;

    for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    {
      dmap_nxsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsx*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
    }

    for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
    {
      dmap_nysl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysy*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
    }

    for (CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
    {
      dmap_nzsl_gp_unit[p->first] += 1/length*(p->second);
      dmap_nzsl_gp_unit[p->first] -= 1/(length*length*length)*szsz*(p->second);
      dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsz*(p->second);
      dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysz*(p->second);
    }

    // add everything to dgapgp
    for (CI p=dmap_nxsl_gp_unit.begin();p!=dmap_nxsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

    for (CI p=dmap_nysl_gp_unit.begin();p!=dmap_nysl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

    for (CI p=dmap_nzsl_gp_unit.begin();p!=dmap_nzsl_gp_unit.end();++p)
      dgapgp[p->first] += (mgpx[2]-sgpx[2]) *(p->second);

    for (int z=0;z<nrow;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[smrtrnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * smrtrnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,1) * smrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[mmrtrnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dpmxigp[0].begin();p!=dpmxigp[0].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mmrtrnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dpmxigp[1].begin();p!=dpmxigp[1].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,1) * mmrtrnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************
    
    // compute cell gap vector *******************************************
    // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
    if (shapefcn_ == INPAR::MORTAR::shape_standard &&
        (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin))
    {
      for (int j=0;j<nrow;++j)
      {
        // add current Gauss point's contribution to gseg
        (*gseg)(j) += lmval[j]*gap*jac*wgt;
      }
    }
    
    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    else if (shapefcn_ == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      for (int j=0;j<nintrow;++j)
      {
        // add current Gauss point's contribution to gseg
        (*gseg)(j) += lmintval[j]*gap*jac*wgt;
      }
    }
    
    // CASE 4: Dual LM shape functions and quadratic interpolation
    else if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      for (int j=0;j<nrow;++j)
      {
        // add current Gauss point's contribution to gseg
        (*gseg)(j) += lmval[j]*gap*jac*wgt;
      }
    }
    
    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic contact!");
    }
    // compute cell gap vector *******************************************

    // compute cell D/M linearization ************************************
    // CASE 1: Standard LM shape functions and quadratic interpolation
    if (shapefcn_ == INPAR::MORTAR::shape_standard &&
        lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        // integrate LinM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmderiv(j, 1)*mval[k]*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval[j]*mderiv(k, 0)*jac;
          for (CI p=dpmxigp[0].begin(); p!=dpmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          for (CI p=dpmxigp[1].begin(); p!=dpmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes

        // integrate LinD
        for (int k=0; k<nrow; ++k)
        {
          // global master node ID
          int sgid = sele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*sval[k]*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmderiv(j, 1)*sval[k]*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmval[j]*sderiv(k, 0)*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmval[j]*sderiv(k, 1)*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*sval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);
        } // loop over slave nodes
      }
    }
    
    // CASE 2: Standard LM shape functions and linear interpolation
    // (this has to be treated seperately here for LinDM because of bound)
    else if (shapefcn_ == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_lin_lin)
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
        bool jbound = mymrtrnode->IsOnBound();
        
        // node j is boundary node
        if (jbound)
        {
          // do nothing as respective D and M entries are zero anyway
        }
        
        // node j is NO boundary node
        else
        {
          // integrate LinM
          for (int k=0; k<ncol; ++k)
          {
            // global master node ID
            int mgid = mele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*lmderiv(j, 0)*mval[k]*jac;
            for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            fac = wgt*lmderiv(j, 1)*mval[k]*jac;
            for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NMaster) - master GP coordinates
            fac = wgt*lmval[j]*mderiv(k, 0)*jac;
            for (CI p=dpmxigp[0].begin(); p!=dpmxigp[0].end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            fac = wgt*lmval[j]*mderiv(k, 1)*jac;
            for (CI p=dpmxigp[1].begin(); p!=dpmxigp[1].end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt*lmval[j]*mval[k];
            for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
              dmmap_jk[p->first] += fac*(p->second);
          } // loop over master nodes

          // integrate LinD
          for (int k=0; k<nrow; ++k)
          {
            MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[k]);
            if (!mymrtrnode2) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
            bool kbound = mymrtrnode2->IsOnBound();

            // global master node ID
            int sgid = mymrtrnode2->Id();
            double fac = 0.0;

            // node k is boundary node
            if (kbound)
            {
              // move entry to derivM (with minus sign)
              // get the correct map as a reference
              std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

              // (1) Lin(Phi) - dual shape functions
              // this vanishes here since there are no deformation-dependent dual functions

              // (2) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmderiv(j, 0)*sval[k]*jac;
              for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              fac = wgt*lmderiv(j, 1)*sval[k]*jac;
              for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              // (3) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmval[j]*sderiv(k, 0)*jac;
              for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              fac = wgt*lmval[j]*sderiv(k, 1)*jac;
              for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);

              // (4) Lin(dsxideta) - intcell GP Jacobian
              fac = wgt*lmval[j]*sval[k];
              for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
                dmmap_jk[p->first] -= fac*(p->second);
            }

            // node k is NO boundary node
            else
            {
              // get the correct map as a reference
              std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

              // (1) Lin(Phi) - dual shape functions
              // this vanishes here since there are no deformation-dependent dual functions

              // (2) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmderiv(j, 0)*sval[k]*jac;
              for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              fac = wgt*lmderiv(j, 1)*sval[k]*jac;
              for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              // (3) Lin(NSlave) - slave GP coordinates
              fac = wgt*lmval[j]*sderiv(k, 0)*jac;
              for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              fac = wgt*lmval[j]*sderiv(k, 1)*jac;
              for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);

              // (4) Lin(dsxideta) - intcell GP Jacobian
              fac = wgt*lmval[j]*sval[k];
              for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
                ddmap_jk[p->first] += fac*(p->second);
            }

          } // loop over slave nodes
        }
      }
    }
    
    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    else if (shapefcn_ == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      for (int j=0;j<nintrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(myintnodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

        // integrate LinM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmintderiv(j, 0)*mval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmintderiv(j, 1)*mval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmintval[j]*mderiv(k, 0)*jac;
          for (CI p=dpmxigp[0].begin(); p!=dpmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmintval[j]*mderiv(k, 1)*jac;
          for (CI p=dpmxigp[1].begin(); p!=dpmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmintval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes

        // integrate LinD
        for (int k=0; k<nrow; ++k)
        {
          // global master node ID
          int sgid = sele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmintderiv(j, 0)*sval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmintderiv(j, 1)*sval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmintval[j]*sderiv(k, 0)*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          fac = wgt*lmintval[j]*sderiv(k, 1)*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmintval[j]*sval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            ddmap_jk[p->first] += fac*(p->second);
        } // loop over slave nodes
      }
    }
    
    // CASE 4: Dual LM shape functions and quadratic interpolation
    else if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
        int sgid = mymrtrnode->Id();

        // for dual shape functions ddmap_jj and dmmap_jk can be calculated together
        std::map<int,double>& ddmap_jj = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

        // integrate LinM and LinD
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0; m<nrow; ++m)
            {
              if (dualquad3d) fac = wgt*svalmod[m]*mval[k]*jac;
              else            fac = wgt*sval[m]*mval[k]*jac;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                ddmap_jj[p->first] += fac*(p->second);
              }
            }

          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*lmderiv(j, 1)*mval[k]*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval[j]*mderiv(k, 0)*jac;
          for (CI p=dpmxigp[0].begin(); p!=dpmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*lmval[j]*mderiv(k, 1)*jac;
          for (CI p=dpmxigp[1].begin(); p!=dpmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*lmval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            ddmap_jj[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // loop over slave nodes
    }
    
    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic contact!");
    }
    // compute cell D/M linearization ************************************

    // compute cell gap linearization ************************************
    // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
    if (shapefcn_ == INPAR::MORTAR::shape_standard &&
        (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin))
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
        double fac = 0.0;

        // get the corresponding map as a reference
        std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions

        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmderiv(j,0)*gap*jac;
        for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
          dgmap[p->first] += fac*(p->second);
        
        fac = wgt*lmderiv(j,1)*gap*jac;
        for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*lmval[j]*jac;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*lmval[j]*gap;
        for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
          dgmap[p->first] += fac*(p->second);
      }
    }
    
    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    else if (shapefcn_ == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      for (int j=0;j<nintrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(myintnodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
        double fac = 0.0;
  
        // get the corresponding map as a reference
        std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();
  
        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions
        
        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmintderiv(j,0)*gap*jac;
        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgmap[p->first] += fac*(p->second);
        
        fac = wgt*lmintderiv(j,1)*gap*jac;
        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgmap[p->first] += fac*(p->second);
  
        // (3) Lin(g) - gap function
        fac = wgt*lmintval[j]*jac;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);
  
        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*lmintval[j]*gap;
        for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
          dgmap[p->first] += fac*(p->second);
      }
    }
    
    // CASE 4: Dual LM shape functions and quadratic interpolation
    else if ((shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin) &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
        double fac = 0.0;

        // get the corresponding map as a reference
        std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

        // (1) Lin(Phi) - dual shape functions
        if (duallin)
          for (int m=0;m<nrow;++m)
          {
            if (dualquad3d) fac = wgt*svalmod[m]*gap*jac;
            else            fac = wgt*sval[m]*gap*jac;

            for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
              dgmap[p->first] += fac*(p->second);
          }

        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*lmderiv(j,0)*gap*jac;
        for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
          dgmap[p->first] += fac*(p->second);
        
        fac = wgt*lmderiv(j,1)*gap*jac;
        for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
          dgmap[p->first] += fac*(p->second);

        // (3) Lin(g) - gap function
        fac = wgt*lmval[j]*jac;
        for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
          dgmap[p->first] += fac*(p->second);

        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*lmval[j]*gap;
        for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
          dgmap[p->first] += fac*(p->second);
      }
    }
    
    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic contact!");
    }
    // compute cell gap linearization ************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize D                                farah 01/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::FastIntegration(
     MORTAR::MortarElement& sele,
     std::vector<MORTAR::MortarElement*> meles,
     INPAR::MORTAR::LagMultQuad lmtype,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg,
     Teuchos::RCP<Epetra_SerialDenseVector> wseg,
     bool *boundary_ele)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");

  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  for (int i=0;i<(int)meles.size();++i)
  {
    if ((!sele.IsSlave()) || (meles[i]->IsSlave()))
      dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  }

  //consider entire slave element --> parameter space [-1,1]
  double sxia=-1.0;
  double sxib=1.0;

  // number of nodes (slave, master)
  int nrow      =   sele.NumNode();
  int ncol      =   meles[0]->NumNode();
  int ndof      =   Dim();
  int nodemaster  = meles[0]->NumNode();
  int mndof     =   static_cast<MORTAR::MortarNode*>(meles[0]->Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseMatrix ssecderiv(nrow,1);

  // get slave nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);

  // nodal coords from previous time step and lagrange mulitplier
  Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold;
  Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult;

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is only necessary for quadratic dual shape functions in 2D
  //bool duallin = false; // --> coming soon
  std::vector<std::vector<std::map<int,double> > > dualmap(nrow,std::vector<std::map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual) && (sele.Shape()==MORTAR::MortarElement::line3))
  {
    //duallin=true;
    sele.DerivShapeDual(dualmap);
  }

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  if (lmtype == INPAR::MORTAR::lagmult_lin_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false; // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
  }

  // get numerical integration type
  const Teuchos::ParameterList& smortar  = DRT::Problem::Instance()->MortarCouplingParams();
  INPAR::MORTAR::IntType inttype =
    DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(smortar,"INTTYPE");

  //************************************************************************
  //Boundary Segmentation check -- HasProj()-check
  //************************************************************************
  double sxi_test[2] = {0.0, 0.0};
  bool proj_test=false;

  double glob_test[3] = {0.0, 0.0, 0.0};

  DRT::Node** mynodes_test = sele.Nodes();
  if (!mynodes_test) dserror("ERROR: HasProjStatus: Null pointer!");

  if(sele.Shape()==DRT::Element::line2)
  {
    for (int s_test=0;s_test<2;++s_test)
    {
      if (s_test==0) sxi_test[0]=-1.0;
      else sxi_test[0]=1.0;

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(2);
        projector_test.ProjectGaussPoint(sele,sxi_test,*meles[bs_test],mxi_test);

        if ((mxi_test[0]>=-1.0) && (mxi_test[0]<=1.0))
        {
          //get hasproj
          sele.LocalToGlobal(sxi_test,glob_test,0);
          for (int ii=0;ii<sele.NumNode();++ii)
          {
            MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
            if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

            if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
              mycnode_test->HasProj()=true;
          }

          glob_test[0]=0.0;
          glob_test[1]=0.0;
          glob_test[2]=0.0;

          proj_test=true;
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  else if (sele.Shape()==DRT::Element::line3)
  {
    for (int s_test=0;s_test<3;++s_test)
    {
      if (s_test==0) sxi_test[0]=-1.0;
      else if (s_test==1) sxi_test[0]=0.0;
      else if (s_test==2) sxi_test[0]=1.0;

      proj_test=false;
      for (int bs_test=0;bs_test<(int)meles.size();++bs_test)
      {
        double mxi_test[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector_test(2);
        projector_test.ProjectGaussPoint(sele,sxi_test,*meles[bs_test],mxi_test);

        if ((mxi_test[0]>=-1.0) && (mxi_test[0]<=1.0))
        {
          //get hasproj
          sele.LocalToGlobal(sxi_test,glob_test,0);
          for (int ii=0;ii<sele.NumNode();++ii)
          {
            MORTAR::MortarNode* mycnode_test = static_cast<MORTAR::MortarNode*> (mynodes_test[ii]);
            if (!mycnode_test) dserror("ERROR: HasProjStatus: Null pointer!");

            if (glob_test[0]==mycnode_test->xspatial()[0] && glob_test[1]==mycnode_test->xspatial()[1] && glob_test[2]==mycnode_test->xspatial()[2])
              mycnode_test->HasProj()=true;
          }

          glob_test[0]=0.0;
          glob_test[1]=0.0;
          glob_test[2]=0.0;

          proj_test=true;
        }
      }
      if(proj_test==false) *boundary_ele=true;
    }
  }
  else
  {
    dserror("No valid element type for slave discretization!");
  }



  if (*boundary_ele==false || inttype==INPAR::MORTAR::inttype_fast)
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

      // evaluate Lagrange multiplier shape functions (on slave element)
      if (linlm)
        sele.EvaluateShapeLagMultLin(shapefcn_,sxi,lmval,lmderiv,nrow);
      else
        sele.EvaluateShapeLagMult(shapefcn_,sxi,lmval,lmderiv,nrow);

      // evaluate trace space shape functions
      sele.EvaluateShape(sxi,sval,sderiv,nrow);

      //****************************************************************************************************************
      //                loop over all Master Elements
      //****************************************************************************************************************
      for (int nummaster=0;nummaster<(int)meles.size();++nummaster)
      {
        // get master nodal coords for Jacobian / GP evaluation
        LINALG::SerialDenseMatrix mcoord(3,meles[nummaster]->NumNode());
        meles[nummaster]->GetNodalCoords(mcoord);

        // project Gauss point onto master element
        double mxi[2] = {0.0, 0.0};
        MORTAR::MortarProjector projector(2);
        projector.ProjectGaussPoint(sele,sxi,*meles[nummaster],mxi);

        // evaluate trace space shape functions
        meles[nummaster]->EvaluateShape(mxi,mval,mderiv,ncol);

        // evaluate the two slave side Jacobians
        double dxdsxi = sele.Jacobian(sxi);
        double dsxideta = -0.5*sxia + 0.5*sxib;

        if ((mxi[0]>=-1.0) && (mxi[0]<=1.0) && (kink_projection==false))
        {
          kink_projection=true;

          // compute segment D/M matrix ****************************************
          // standard shape functions
          if (shapefcn_ == INPAR::MORTAR::shape_standard)
          {
            // loop over all mseg matrix entries
            // !!! nrow represents the slave Lagrange multipliers !!!
            // !!! ncol represents the master dofs                !!!
            // (this DOES matter here for mseg, as it might
            // sometimes be rectangular, not quadratic!)
            for (int j=0; j<nrow*ndof; ++j)
            {
              // integrate mseg
              for (int k=0; k<nodemaster*mndof; ++k)    //loop over dofs of the considered master element
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/mndof);

                DRT::Node** mnodes = meles[nummaster]->Nodes();
                MORTAR::MortarNode* mymnode = static_cast<MORTAR::MortarNode*> (mnodes[kindex]);
                if (!mymnode) dserror("ERROR: Null pointer!");

                // multiply the two shape functions
                double prod = lmval[jindex]*mval[kindex];   //Lagrange-Shapefunction * Masterelement-Shapefunction

                // isolate the mseg and dseg entries to be filled
                // (both the main diagonal and every other secondary diagonal)
                // and add current Gauss point's contribution to mseg and dseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*mndof)))
                {
                  (*mseg)(j, k+nummaster*nodemaster*mndof) += prod*dxdsxi*dsxideta*wgt;
                }
              }
            }
          }
          //***********************
          // dual shape functions
          //***********************
          else if (shapefcn_ == INPAR::MORTAR::shape_dual)
          {
            //dserror("ERROR: noch nicht vervollständigt");
            // loop over all mseg matrix entries
            // nrow represents the slave Lagrange multipliers !!!
            // ncol represents the master dofs !!!
            // (this DOES matter here for mseg, as it might
            // sometimes be rectangular, not quadratic!)
            for (int j=0;j<nrow*ndof;++j)
            {
              // for dual shape functions we can make use
              // of the row summing lemma: D_jj = Sum(k) M_jk
              // hence, they can be combined into one single loop

              // integrate mseg and dseg (no boundary modification)
              for (int k=0;k<nodemaster*mndof;++k)
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/ndof);

                // multiply the two shape functions
                double prod = lmval[jindex]*mval[kindex];

                // isolate the mseg and dseg entries to be filled
                // (both the main diagonal and every other secondary diagonal for mseg)
                // (only the main diagonal for dseg)
                // and add current Gauss point's contribution to mseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                {
                  (*mseg)(j, k+nummaster*nodemaster*mndof) += prod*dxdsxi*dsxideta*wgt;
                  if(!bound)
                  {
                    (*dseg)(j,j) += prod*dxdsxi*dsxideta*wgt;
                  }
                }
              }
            }
          }//end - dualshapefunction // shapefcn_ switch

          if (shapefcn_ == INPAR::MORTAR::shape_standard)
          {
            //integrate dseg
            for (int j=0; j<nrow*ndof; ++j)
            {
              for (int k=0; k<nrow*ndof; ++k)
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/ndof);

                // multiply the two shape functions
                double prod = lmval[jindex]*sval[kindex];

                // isolate the mseg and dseg entries to be filled
                // (both the main diagonal and every other secondary diagonal)
                // and add current Gauss point's contribution to mseg and dseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                {
                  (*dseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
                }
              }
            }
          }
          else if (shapefcn_ == INPAR::MORTAR::shape_dual)
          {
            for (int j=0;j<nrow*ndof;++j)
            {
              // integrate dseg (boundary modification)
              if (bound)
              {
                MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[(int)(j/ndof)]);
                if (!mymrtrnode) dserror("ERROR: Fastintegration: Null pointer!");
                bool j_boundnode = mymrtrnode->IsOnBound();

                for (int k=0;k<nrow*ndof;++k)
                {
                  MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[(int)(k/ndof)]);
                  if (!mymrtrnode2) dserror("ERROR: FastIntegration: Null pointer!");
                  bool k_boundnode = mymrtrnode2->IsOnBound();

                  int jindex = (int)(j/ndof);
                  int kindex = (int)(k/ndof);

                  // do not assemble off-diagonal terms if j,k are both non-boundary nodes
                  if (!j_boundnode && !k_boundnode && (jindex!=kindex)) continue;

                  // multiply the two shape functions
                  double prod = lmval[jindex]*sval[kindex];

                  // isolate the dseg entries to be filled
                  // (both the main diagonal and every other secondary diagonal)
                  // and add current Gauss point's contribution to dseg
                  if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                    (*dseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
                }
              }
            }
          }

          // evaluate 2nd deriv of trace space shape functions (on slave element)
          sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);  //--> lin: 0

          // build interpolation of slave GP normal and coordinates
          double gpn[3] = {0.0,0.0,0.0};
          double sgpx[3] = {0.0, 0.0, 0.0};
          for (int i=0;i<nrow;++i)  //loop over all slave nodes
          {
            MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*> (mynodes[i]);
            gpn[0]+=sval[i]*mymrtrnode->MoData().n()[0];
            gpn[1]+=sval[i]*mymrtrnode->MoData().n()[1];
            gpn[2]+=sval[i]*mymrtrnode->MoData().n()[2];

            sgpx[0]+=sval[i]*scoord(0,i);
            sgpx[1]+=sval[i]*scoord(1,i);
            sgpx[2]+=sval[i]*scoord(2,i);
          }

          // normalize interpolated GP normal back to length 1.0 !!!
          double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
          if (length<1.0e-12) dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

          for (int i=0;i<3;++i)
            gpn[i]/=length;

          // build interpolation of master GP coordinates
          double mgpx[3] = {0.0, 0.0, 0.0};
          for (int i=0;i<ncol;++i)
          {
            mgpx[0]+=mval[i]*mcoord(0,i);
            mgpx[1]+=mval[i]*mcoord(1,i);
            mgpx[2]+=mval[i]*mcoord(2,i);
          }

          // build gap function at current GP
          double gap = 0.0;
          for (int i=0;i<3;++i)
            gap+=(mgpx[i]-sgpx[i])*gpn[i];


          // get directional derivatives of sxia, sxib, mxia, mxib --> derivatives of mxia/mxib not required
          std::vector<std::map<int,double> > ximaps(4);
          bool startslave = true;
          bool endslave = true;
          double mxia=-0.1; //--> arbitrary value
          double mxib=0.1;  //--> arbitrary value
          DerivXiAB2D(sele,sxia,sxib,*meles[nummaster],mxia,mxib,ximaps,startslave,endslave);

          // evalute the GP slave coordinate derivatives --> no entries
          std::map<int,double> dsxigp;
          for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
            dsxigp[p->first]    = 0.0;

          // evalute the GP master coordinate derivatives
          std::map<int,double> dmxigp;
          DerivXiGP2D(sele,*meles[nummaster],sxi[0],mxi[0],dsxigp,dmxigp);

          // evaluate the Jacobian derivative
          std::map<int,double> derivjac;
          sele.DerivJacobian(sxi,derivjac); //direct derivative if xi^1_g does not change

          // evaluate the GP gap function derivatives
          std::map<int,double> dgapgp;

          // we need the participating slave and master nodes
          DRT::Node** snodes = sele.Nodes();
          DRT::Node** mnodes = meles[nummaster]->Nodes();
          std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
          std::vector<MORTAR::MortarNode*> mmrtrnodes(meles[nummaster]->NumNode());

          //check whether the pointer to the nodes exist
          for (int i=0;i<nrow;++i)
          {
            smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
            if (!smrtrnodes[i]) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
          }

          for (int i=0;i<ncol;++i)
          {
            mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
            if (!mmrtrnodes[i]) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
          }

          // build directional derivative of slave GP normal (non-unit)
          std::map<int,double> dmap_nxsl_gp;
          std::map<int,double> dmap_nysl_gp;

          for (int i=0;i<nrow;++i)
          {
            std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
            std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

            //build nc --> gp nomal
            for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
              dmap_nxsl_gp[p->first] += sval[i]*(p->second);
            for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
              dmap_nysl_gp[p->first] += sval[i]*(p->second);
          }

          // build directional derivative of slave GP normal (unit)
          std::map<int,double> dmap_nxsl_gp_unit;
          std::map<int,double> dmap_nysl_gp_unit;

          double ll = length*length;
          double sxsx = gpn[0]*gpn[0]*ll;
          double sxsy = gpn[0]*gpn[1]*ll;
          double sysy = gpn[1]*gpn[1]*ll;

          for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
          {
            dmap_nxsl_gp_unit[p->first] += 1/length*(p->second);
            dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsx*(p->second);
            dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
          }

          for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
          {
            dmap_nysl_gp_unit[p->first] += 1/length*(p->second);
            dmap_nysl_gp_unit[p->first] -= 1/(length*length*length)*sysy*(p->second);
            dmap_nxsl_gp_unit[p->first] -= 1/(length*length*length)*sxsy*(p->second);
          }

          // add everything to dgapgp
          for (CI p=dmap_nxsl_gp_unit.begin();p!=dmap_nxsl_gp_unit.end();++p)
            dgapgp[p->first] += (mgpx[0]-sgpx[0]) * (p->second);

          for (CI p=dmap_nysl_gp_unit.begin();p!=dmap_nysl_gp_unit.end();++p)
            dgapgp[p->first] += (mgpx[1]-sgpx[1]) * (p->second);

          for (int z=0;z<nrow;++z)
          {
            for (int k=0;k<2;++k)
            {
              dgapgp[smrtrnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];
            }
          }

          for (int z=0;z<ncol;++z)
          {
            for (int k=0;k<2;++k)
            {
              dgapgp[mmrtrnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

              for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
                dgapgp[p->first] += gpn[k] * mderiv(z,0) * mmrtrnodes[z]->xspatial()[k] * (p->second);
            }
          }

          // evaluate linearizations *******************************************

          // compute gap vector ****************************************
          // loop over all gseg vector entries
          // nrow represents the slave side dofs !!!
          for (int j=0;j<nrow;++j)
          {
            double prod = lmval[j]*gap;

            // add current Gauss point's contribution to gseg
            (*gseg)(j) += prod*dxdsxi*wgt;
          }

          // compute segment D/M linearization *********************************

          // no linear LM interpolation for quadratic FE
          for (int j=0;j<nrow;++j)
          {
            MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
            if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

            int sgid = mymrtrnode->Id();

            // standard shape functions
            if (shapefcn_ == INPAR::MORTAR::shape_standard)
            {
              // integrate LinM
              for (int k=0; k<ncol; ++k)
              {
                // global master node ID
                int mgid = meles[nummaster]->Nodes()[k]->Id();
                double fac = 0.0;

                // get the correct map as a reference
                std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

                // (1) Lin(Phi) - dual shape functions    --> 0

                // (2) Lin(NSlave) - slave GP coordinates --> 0


                // (3) Lin(NMaster) - master GP coordinates
                fac = wgt*lmval[j]*mderiv(k, 0)*dsxideta*dxdsxi;
                for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
                  dmmap_jk[p->first] += fac*(p->second);

                // (4) Lin(dsxideta) - segment end coordinates --> 0


                // (5) Lin(dxdsxi) - slave GP Jacobian
                fac = wgt*lmval[j]*mval[k]*dsxideta;
                for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
                  dmmap_jk[p->first] += fac*(p->second);

                // (6) Lin(dxdsxi) - slave GP coordinates --> 0
              } // loop over master nodes

              // integrate LinD
              for (int k=0; k<nrow; ++k)
              {
                // global slave node ID
                int sgid = sele.Nodes()[k]->Id();
                double fac = 0.0;

                // get the correct map as a reference
                std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

                // (1) Lin(Phi) - dual shape functions --> 0

                // (2) Lin(NSlave) - slave GP coordinates --> 0

                // (3) Lin(NSlave) - slave GP coordinates --> 0

                // (4) Lin(dsxideta) - segment end coordinates --> 0

                // (5) Lin(dxdsxi) - slave GP Jacobian
                fac = wgt*lmval[j]*sval[k]*dsxideta;
                for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
                  ddmap_jk[p->first] += fac*(p->second);

                // (6) Lin(dxdsxi) - slave GP coordinates --> 0
              } // loop over slave nodes
            }

            // dual shape functions
            else if (shapefcn_ == INPAR::MORTAR::shape_dual)
            {
              // get the D-map as a reference
              std::map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

              // integrate LinM and LinD (NO boundary modification)
              for (int k=0; k<ncol; ++k)
              {
                // global master node ID
                int mgid = meles[nummaster]->Nodes()[k]->Id();
                double fac = 0.0;

                // get the correct map as a reference
                std::map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

                // (1) Lin(Phi) - dual shape functions
                for (int m=0; m<nrow; ++m)
                {
                  fac = wgt*sval[m]*mval[k]*dsxideta*dxdsxi;
                  for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
                  {
                    dmmap_jk[p->first] += fac*(p->second);
                    if (!bound) ddmap_jk[p->first] += fac*(p->second);
                  }
                }

                // (2) Lin(Phi) - slave GP coordinates --> 0

                // (3) Lin(NMaster) - master GP coordinates
                fac = wgt*lmval(j, 0)*mderiv(k, 0)*dsxideta*dxdsxi;
                for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
                {
                  dmmap_jk[p->first] += fac*(p->second);
                  if (!bound) ddmap_jk[p->first] += fac*(p->second);
                }

                // (4) Lin(dsxideta) - segment end coordinates --> 0


                // (5) Lin(dxdsxi) - slave GP Jacobian
                fac = wgt*lmval[j]*mval[k]*dsxideta;
                for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
                {
                  dmmap_jk[p->first] += fac*(p->second);
                  if (!bound) ddmap_jk[p->first] += fac*(p->second);
                }

                // (6) Lin(dxdsxi) - slave GP coordinates --> 0
              } // loop over master nodes
            } // shapefcn_ switch
          } // linlm or not

          // compute segment gap linearization *********************************
          for (int j=0;j<nrow;++j)
          {
            MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
            if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

            double fac = 0.0;

            // get the corresponding map as a reference
            std::map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

            // (1) Lin(Phi) - dual shape functions
            if (shapefcn_ == INPAR::MORTAR::shape_dual)
            {
              for (int m=0;m<nrow;++m)
              {
                fac = wgt*sval[m]*gap*dsxideta*dxdsxi;
                for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
                {
                  dgmap[p->first] += fac*(p->second);
                }
              }
            }

            // (2) Lin(Phi) - slave GP coordinates --> 0

            // (3) Lin(g) - gap function
            fac = wgt*lmval[j]*dsxideta*dxdsxi;
            for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
              dgmap[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - segment end coordinates --> 0

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*lmval[j]*gap*dsxideta;
            for (CI p=derivjac.begin();p!=derivjac.end();++p)
              dgmap[p->first] += fac*(p->second);

            // (6) Lin(dxdsxi) - slave GP coordinates --> 0
          }
        } //Abfrage ob GP auf Mele
      }//End Loop over all Master Elements
    } // End Loop over all GP
  }//boundary_ele abfrage

  return;
}

/*----------------------------------------------------------------------*
 |  Compute penalty scaling factor kappa                      popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateKappaPenalty(MORTAR::MortarElement& sele,
                                                  double* sxia, double* sxib,
                                                  Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateKappaPenalty called without specific shape function defined!");
    
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: IntegrateKappaPenalty called on a non-slave MortarElement!");
  if ((sxia[0]<-1.0) || (sxia[1]<-1.0) || (sxib[0]>1.0) || (sxib[1]>1.0))
    dserror("ERROR: IntegrateKappaPenalty called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,2,true);
  
  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateKappaPenalty: Null pointer!");
    
  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateKappaPenalty: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }
    
  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), 0.0};
    if (Dim()==3) eta[1] = Coordinate(gp,1);
    double wgt = Weight(gp);

    // evaluate shape functions
    if (bound) sele.EvaluateShapeLagMultLin(shapefcn_,eta,val,deriv,nrow);
    else       sele.EvaluateShapeLagMult(shapefcn_,eta,val,deriv,nrow);

    // evaluate the Jacobian det
    double jac = sele.Jacobian(eta);
        
    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += val[j]*jac*wgt;
    }
    // compute cell gap vector *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute penalty scaling factor kappa (3D piecewise lin)   popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateKappaPenalty(MORTAR::MortarElement& sele,
                                                  MORTAR::IntElement& sintele,
                                                  double* sxia, double* sxib,
                                                  Teuchos::RCP<Epetra_SerialDenseVector> gseg,
                                                  INPAR::MORTAR::LagMultQuad lmtype)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ != INPAR::MORTAR::shape_standard)
    dserror("ERROR: IntegrateKappaPenalty -> you should not be here!");
    
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: IntegrateKappaPenalty called on a non-slave MortarElement!");
  if ((sxia[0]<-1.0) || (sxia[1]<-1.0) || (sxib[0]>1.0) || (sxib[1]>1.0))
    dserror("ERROR: IntegrateKappaPenalty called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();
  int nintrow = sintele.NumNode();

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,2,true);
  LINALG::SerialDenseVector intval(nintrow);
  LINALG::SerialDenseMatrix intderiv(nintrow,2,true);

  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), 0.0};
    if (Dim()==3) eta[1] = Coordinate(gp,1);
    double wgt = Weight(gp);

    // map Gauss point back to slave element (affine map)
    double psxi[2] = {0.0, 0.0};
    sintele.MapToParent(eta,psxi);
    
    // evaluate shape functions
    sele.EvaluateShape(psxi,val,deriv,nrow);
    sintele.EvaluateShape(eta,intval,intderiv,nintrow);

    // evaluate the Jacobian det
    double jac = sintele.Jacobian(eta);
    
    // compute cell gap vector *******************************************
    if (lmtype==INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      for (int j=0;j<nintrow;++j)
      {
        // add current Gauss point's contribution to gseg
        (*gseg)(j) += intval[j]*jac*wgt;
      }
    }
    
    else
    {
      dserror("ERROR: Invalid LM interpolation case!");
    }
    // compute cell gap vector *******************************************
  }
  //**********************************************************************
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiAB (2D)               popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiAB2D(MORTAR::MortarElement& sele,
                                    double& sxia, double& sxib,
                                    MORTAR::MortarElement& mele,
                                    double& mxia, double& mxib,
                                    std::vector<std::map<int,double> >& derivxi,
                                    bool& startslave, bool& endslave)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  for (int i=0;i<numsnode;++i)
  {
    smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  for (int i=0;i<nummnode;++i)
  {
    mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
    if (!mmrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  // we also need shape function derivs in A and B
  double psxia[2] = {sxia, 0.0};
  double psxib[2] = {sxib, 0.0};
  double pmxia[2] = {mxia, 0.0};
  double pmxib[2] = {mxib, 0.0};
  LINALG::SerialDenseVector valsxia(numsnode);
  LINALG::SerialDenseVector valsxib(numsnode);
  LINALG::SerialDenseVector valmxia(nummnode);
  LINALG::SerialDenseVector valmxib(nummnode);
  LINALG::SerialDenseMatrix derivsxia(numsnode,1);
  LINALG::SerialDenseMatrix derivsxib(numsnode,1);
  LINALG::SerialDenseMatrix derivmxia(nummnode,1);
  LINALG::SerialDenseMatrix derivmxib(nummnode,1);

  sele.EvaluateShape(psxia,valsxia,derivsxia,numsnode);
  sele.EvaluateShape(psxib,valsxib,derivsxib,numsnode);
  mele.EvaluateShape(pmxia,valmxia,derivmxia,nummnode);
  mele.EvaluateShape(pmxib,valmxib,derivmxib,nummnode);

  // compute factors and leading constants for master
  double cmxia = 0.0;
  double cmxib = 0.0;
  double fac_dxm_a = 0.0;
  double fac_dym_a = 0.0;
  double fac_xmsl_a = 0.0;
  double fac_ymsl_a = 0.0;
  double fac_dxm_b = 0.0;
  double fac_dym_b = 0.0;
  double fac_xmsl_b = 0.0;
  double fac_ymsl_b = 0.0;

  // compute leading constant for DerivXiBMaster if start node = slave node
  if (startslave==true)
  {
    for (int i=0;i<nummnode;++i)
    {
      fac_dxm_b += derivmxib(i,0)*(mmrtrnodes[i]->xspatial()[0]);
      fac_dym_b += derivmxib(i,0)*(mmrtrnodes[i]->xspatial()[1]);
      fac_xmsl_b += valmxib[i]*(mmrtrnodes[i]->xspatial()[0]);
      fac_ymsl_b += valmxib[i]*(mmrtrnodes[i]->xspatial()[1]);
    }

    cmxib = -1/(fac_dxm_b*(smrtrnodes[0]->MoData().n()[1])-fac_dym_b*(smrtrnodes[0]->MoData().n()[0]));
    //std::cout << "cmxib: " << cmxib << endl;

    fac_xmsl_b -= smrtrnodes[0]->xspatial()[0];
    fac_ymsl_b -= smrtrnodes[0]->xspatial()[1];
  }

  // compute leading constant for DerivXiAMaster if end node = slave node
  if (endslave==true)
  {
    for (int i=0;i<nummnode;++i)
    {
      fac_dxm_a += derivmxia(i,0)*(mmrtrnodes[i]->xspatial()[0]);
      fac_dym_a += derivmxia(i,0)*(mmrtrnodes[i]->xspatial()[1]);
      fac_xmsl_a += valmxia[i]*(mmrtrnodes[i]->xspatial()[0]);
      fac_ymsl_a += valmxia[i]*(mmrtrnodes[i]->xspatial()[1]);
    }

    cmxia = -1/(fac_dxm_a*(smrtrnodes[1]->MoData().n()[1])-fac_dym_a*(smrtrnodes[1]->MoData().n()[0]));
    //std::cout << "cmxia: " << cmxia << endl;

    fac_xmsl_a -= smrtrnodes[1]->xspatial()[0];
    fac_ymsl_a -= smrtrnodes[1]->xspatial()[1];
  }

  // compute factors and leading constants for slave
  double csxia = 0.0;
  double csxib = 0.0;
  double fac_dxsl_a = 0.0;
  double fac_dysl_a = 0.0;
  double fac_xslm_a = 0.0;
  double fac_yslm_a = 0.0;
  double fac_dnx_a = 0.0;
  double fac_dny_a = 0.0;
  double fac_nx_a = 0.0;
  double fac_ny_a = 0.0;
  double fac_dxsl_b = 0.0;
  double fac_dysl_b = 0.0;
  double fac_xslm_b = 0.0;
  double fac_yslm_b = 0.0;
  double fac_dnx_b = 0.0;
  double fac_dny_b = 0.0;
  double fac_nx_b = 0.0;
  double fac_ny_b = 0.0;

  // compute leading constant for DerivXiASlave if start node = master node
  if (startslave==false)
  {
    for (int i=0;i<numsnode;++i)
    {
      fac_dxsl_a += derivsxia(i,0)*(smrtrnodes[i]->xspatial()[0]);
      fac_dysl_a += derivsxia(i,0)*(smrtrnodes[i]->xspatial()[1]);
      fac_xslm_a += valsxia[i]*(smrtrnodes[i]->xspatial()[0]);
      fac_yslm_a += valsxia[i]*(smrtrnodes[i]->xspatial()[1]);
      fac_dnx_a  += derivsxia(i,0)*(smrtrnodes[i]->MoData().n()[0]);
      fac_dny_a  += derivsxia(i,0)*(smrtrnodes[i]->MoData().n()[1]);
      fac_nx_a   += valsxia[i]*(smrtrnodes[i]->MoData().n()[0]);
      fac_ny_a   += valsxia[i]*(smrtrnodes[i]->MoData().n()[1]);
    }

    fac_xslm_a -= mmrtrnodes[1]->xspatial()[0];
    fac_yslm_a -= mmrtrnodes[1]->xspatial()[1];

    csxia = -1/(fac_dxsl_a*fac_ny_a - fac_dysl_a*fac_nx_a + fac_xslm_a*fac_dny_a - fac_yslm_a*fac_dnx_a);
    //std::cout << "csxia: " << csxia << endl;
  }

  // compute leading constant for DerivXiBSlave if end node = master node
  if (endslave==false)
  {
    for (int i=0;i<numsnode;++i)
    {
      fac_dxsl_b  += derivsxib(i,0)*(smrtrnodes[i]->xspatial()[0]);
      fac_dysl_b  += derivsxib(i,0)*(smrtrnodes[i]->xspatial()[1]);
      fac_xslm_b += valsxib[i]*(smrtrnodes[i]->xspatial()[0]);
      fac_yslm_b += valsxib[i]*(smrtrnodes[i]->xspatial()[1]);
      fac_dnx_b  += derivsxib(i,0)*(smrtrnodes[i]->MoData().n()[0]);
      fac_dny_b  += derivsxib(i,0)*(smrtrnodes[i]->MoData().n()[1]);
      fac_nx_b   += valsxib[i]*(smrtrnodes[i]->MoData().n()[0]);
      fac_ny_b   += valsxib[i]*(smrtrnodes[i]->MoData().n()[1]);
    }

    fac_xslm_b -= mmrtrnodes[0]->xspatial()[0];
    fac_yslm_b -= mmrtrnodes[0]->xspatial()[1];

    csxib = -1/(fac_dxsl_b*fac_ny_b - fac_dysl_b*fac_nx_b + fac_xslm_b*fac_dny_b - fac_yslm_b*fac_dnx_b);
    //std::cout << "csxib: " << csxib << endl;
  }

  // prepare linearizations
  typedef std::map<int,double>::const_iterator CI;

  // *********************************************************************
  // finally compute Lin(XiAB_master)
  // *********************************************************************
  // build DerivXiBMaster if start node = slave node
  if (startslave==true)
  {
    std::map<int,double> dmap_mxib;
    std::map<int,double>& nxmap_b = static_cast<CONTACT::CoNode*>(smrtrnodes[0])->CoData().GetDerivN()[0];
    std::map<int,double>& nymap_b = static_cast<CONTACT::CoNode*>(smrtrnodes[0])->CoData().GetDerivN()[1];

    // add derivative of slave node coordinates
    dmap_mxib[smrtrnodes[0]->Dofs()[0]] -= smrtrnodes[0]->MoData().n()[1];
    dmap_mxib[smrtrnodes[0]->Dofs()[1]] += smrtrnodes[0]->MoData().n()[0];

    // add derivatives of master node coordinates
    for (int i=0;i<nummnode;++i)
    {
      dmap_mxib[mmrtrnodes[i]->Dofs()[0]] += valmxib[i]*(smrtrnodes[0]->MoData().n()[1]);
      dmap_mxib[mmrtrnodes[i]->Dofs()[1]] -= valmxib[i]*(smrtrnodes[0]->MoData().n()[0]);
    }

    // add derivative of slave node normal
    for (CI p=nxmap_b.begin();p!=nxmap_b.end();++p)
      dmap_mxib[p->first] -= fac_ymsl_b*(p->second);
    for (CI p=nymap_b.begin();p!=nymap_b.end();++p)
      dmap_mxib[p->first] += fac_xmsl_b*(p->second);

    // multiply all entries with cmxib
    for (CI p=dmap_mxib.begin();p!=dmap_mxib.end();++p)
      dmap_mxib[p->first] = cmxib*(p->second);

    // return map to DerivM() method
    derivxi[3] = dmap_mxib;
  }

  // build DerivXiAMaster if end node = slave node
  if (endslave==true)
  {
    std::map<int,double> dmap_mxia;
    std::map<int,double>& nxmap_a = static_cast<CONTACT::CoNode*>(smrtrnodes[1])->CoData().GetDerivN()[0];
    std::map<int,double>& nymap_a = static_cast<CONTACT::CoNode*>(smrtrnodes[1])->CoData().GetDerivN()[1];

    // add derivative of slave node coordinates
    dmap_mxia[smrtrnodes[1]->Dofs()[0]] -= smrtrnodes[1]->MoData().n()[1];
    dmap_mxia[smrtrnodes[1]->Dofs()[1]] += smrtrnodes[1]->MoData().n()[0];

    // add derivatives of master node coordinates
    for (int i=0;i<nummnode;++i)
    {
      dmap_mxia[mmrtrnodes[i]->Dofs()[0]] += valmxia[i]*(smrtrnodes[1]->MoData().n()[1]);
      dmap_mxia[mmrtrnodes[i]->Dofs()[1]] -= valmxia[i]*(smrtrnodes[1]->MoData().n()[0]);
    }

    // add derivative of slave node normal
    for (CI p=nxmap_a.begin();p!=nxmap_a.end();++p)
      dmap_mxia[p->first] -= fac_ymsl_a*(p->second);
    for (CI p=nymap_a.begin();p!=nymap_a.end();++p)
      dmap_mxia[p->first] += fac_xmsl_a*(p->second);

    // multiply all entries with cmxia
    for (CI p=dmap_mxia.begin();p!=dmap_mxia.end();++p)
      dmap_mxia[p->first] = cmxia*(p->second);

    // return map to DerivM() method
    derivxi[2] = dmap_mxia;
  }

  // *********************************************************************
  // finally compute Lin(XiAB_slave)
  // *********************************************************************
  // build DerivXiASlave if start node = master node
  if (startslave==false)
  {
    std::map<int,double> dmap_sxia;

    // add derivative of master node coordinates
    dmap_sxia[mmrtrnodes[1]->Dofs()[0]] -= fac_ny_a;
    dmap_sxia[mmrtrnodes[1]->Dofs()[1]] += fac_nx_a;

    // add derivatives of slave node coordinates
    for (int i=0;i<numsnode;++i)
    {
      dmap_sxia[smrtrnodes[i]->Dofs()[0]] += valsxia[i]*fac_ny_a;
      dmap_sxia[smrtrnodes[i]->Dofs()[1]] -= valsxia[i]*fac_nx_a;
    }

    // add derivatives of slave node normals
    for (int i=0;i<numsnode;++i)
    {
      std::map<int,double>& nxmap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      std::map<int,double>& nymap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

      for (CI p=nxmap_curr.begin();p!=nxmap_curr.end();++p)
        dmap_sxia[p->first] -= valsxia[i]*fac_yslm_a*(p->second);
      for (CI p=nymap_curr.begin();p!=nymap_curr.end();++p)
        dmap_sxia[p->first] += valsxia[i]*fac_xslm_a*(p->second);
    }

    // multiply all entries with csxia
    for (CI p=dmap_sxia.begin();p!=dmap_sxia.end();++p)
      dmap_sxia[p->first] = csxia*(p->second);

    // return map to DerivM() method
    derivxi[0] = dmap_sxia;
  }

  // build DerivXiBSlave if end node = master node
  if (endslave==false)
  {
    std::map<int,double> dmap_sxib;

    // add derivative of master node coordinates
    dmap_sxib[mmrtrnodes[0]->Dofs()[0]] -= fac_ny_b;
    dmap_sxib[mmrtrnodes[0]->Dofs()[1]] += fac_nx_b;

    // add derivatives of slave node coordinates
    for (int i=0;i<numsnode;++i)
    {
      dmap_sxib[smrtrnodes[i]->Dofs()[0]] += valsxib[i]*fac_ny_b;
      dmap_sxib[smrtrnodes[i]->Dofs()[1]] -= valsxib[i]*fac_nx_b;
    }

    // add derivatives of slave node normals
    for (int i=0;i<numsnode;++i)
    {
      std::map<int,double>& nxmap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      std::map<int,double>& nymap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

      for (CI p=nxmap_curr.begin();p!=nxmap_curr.end();++p)
        dmap_sxib[p->first] -= valsxib[i]*fac_yslm_b*(p->second);
      for (CI p=nymap_curr.begin();p!=nymap_curr.end();++p)
        dmap_sxib[p->first] += valsxib[i]*fac_xslm_b*(p->second);
    }

    // multiply all entries with csxib
    for (CI p=dmap_sxib.begin();p!=dmap_sxib.end();++p)
      dmap_sxib[p->first] = csxib*(p->second);

    // return map to DerivM() method
    derivxi[1] = dmap_sxib;
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiGP master (2D)        popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP2D(MORTAR::MortarElement& sele,
                                    MORTAR::MortarElement& mele,
                                    double& sxigp, double& mxigp,
                                    const std::map<int,double>& derivsxi,
                                    std::map<int,double>& derivmxi)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  for (int i=0;i<numsnode;++i)
  {
    smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  for (int i=0;i<nummnode;++i)
  {
    mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
    if (!mmrtrnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  // we also need shape function derivs in A and B
  double psxigp[2] = {sxigp, 0.0};
  double pmxigp[2] = {mxigp, 0.0};
  LINALG::SerialDenseVector valsxigp(numsnode);
  LINALG::SerialDenseVector valmxigp(nummnode);
  LINALG::SerialDenseMatrix derivsxigp(numsnode,1);
  LINALG::SerialDenseMatrix derivmxigp(nummnode,1);

  sele.EvaluateShape(psxigp,valsxigp,derivsxigp,numsnode);
  mele.EvaluateShape(pmxigp,valmxigp,derivmxigp,nummnode);

  // we also need the GP slave coordinates + normal
  double sgpn[3] = {0.0,0.0,0.0};
  double sgpx[3] = {0.0,0.0,0.0};
  for (int i=0;i<numsnode;++i)
  {
    sgpn[0]+=valsxigp[i]*smrtrnodes[i]->MoData().n()[0];
    sgpn[1]+=valsxigp[i]*smrtrnodes[i]->MoData().n()[1];
    sgpn[2]+=valsxigp[i]*smrtrnodes[i]->MoData().n()[2];

    sgpx[0]+=valsxigp[i]*smrtrnodes[i]->xspatial()[0];
    sgpx[1]+=valsxigp[i]*smrtrnodes[i]->xspatial()[1];
    sgpx[2]+=valsxigp[i]*smrtrnodes[i]->xspatial()[2];
  }

  // FIXME: This does not have to be the UNIT normal (see 3D)!
  // The reason for this is that we linearize the Gauss point
  // projection from slave to master side here and this condition
  // only includes the Gauss point normal in a cross product.
  // When looking at MortarProjector::ProjectGaussPoint, one can see
  // that we do NOT use a unit normal there, either. Thus, why here?
  // First results suggest that it really makes no difference!

  // normalize interpolated GP normal back to length 1.0 !!!
  double length = sqrt(sgpn[0]*sgpn[0]+sgpn[1]*sgpn[1]+sgpn[2]*sgpn[2]);
  if (length<1.0e-12) dserror("ERROR: DerivXiGP2D: Divide by zero!");
  for (int i=0;i<3;++i) sgpn[i]/=length;

  // compute factors and leading constants for master
  double cmxigp = 0.0;
  double fac_dxm_gp = 0.0;
  double fac_dym_gp = 0.0;
  double fac_xmsl_gp = 0.0;
  double fac_ymsl_gp = 0.0;

  for (int i=0;i<nummnode;++i)
  {
    fac_dxm_gp += derivmxigp(i,0)*(mmrtrnodes[i]->xspatial()[0]);
    fac_dym_gp += derivmxigp(i,0)*(mmrtrnodes[i]->xspatial()[1]);

    fac_xmsl_gp += valmxigp[i]*(mmrtrnodes[i]->xspatial()[0]);
    fac_ymsl_gp += valmxigp[i]*(mmrtrnodes[i]->xspatial()[1]);
  }

  cmxigp = -1/(fac_dxm_gp*sgpn[1]-fac_dym_gp*sgpn[0]);
  //std::cout << "cmxigp: " << cmxigp << endl;

  fac_xmsl_gp -= sgpx[0];
  fac_ymsl_gp -= sgpx[1];

  // prepare linearization
  typedef std::map<int,double>::const_iterator CI;

  // build directional derivative of slave GP coordinates
  std::map<int,double> dmap_xsl_gp;
  std::map<int,double> dmap_ysl_gp;

  for (int i=0;i<numsnode;++i)
  {
    dmap_xsl_gp[smrtrnodes[i]->Dofs()[0]] += valsxigp[i];
    dmap_ysl_gp[smrtrnodes[i]->Dofs()[1]] += valsxigp[i];

    for (CI p=derivsxi.begin();p!=derivsxi.end();++p)
    {
      double facx = derivsxigp(i,0)*(smrtrnodes[i]->xspatial()[0]);
      double facy = derivsxigp(i,0)*(smrtrnodes[i]->xspatial()[1]);
      dmap_xsl_gp[p->first] += facx*(p->second);
      dmap_ysl_gp[p->first] += facy*(p->second);
    }
  }

  // build directional derivative of slave GP normal
  std::map<int,double> dmap_nxsl_gp;
  std::map<int,double> dmap_nysl_gp;

  double sgpnmod[3] = {0.0,0.0,0.0};
  for (int i=0;i<3;++i) sgpnmod[i]=sgpn[i]*length;

  std::map<int,double> dmap_nxsl_gp_mod;
  std::map<int,double> dmap_nysl_gp_mod;

  for (int i=0;i<numsnode;++i)
  {
    std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
    std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

    for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
      dmap_nxsl_gp_mod[p->first] += valsxigp[i]*(p->second);
    for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
      dmap_nysl_gp_mod[p->first] += valsxigp[i]*(p->second);

    for (CI p=derivsxi.begin();p!=derivsxi.end();++p)
    {
      double valx =  derivsxigp(i,0)*smrtrnodes[i]->MoData().n()[0];
      dmap_nxsl_gp_mod[p->first] += valx*(p->second);
      double valy =  derivsxigp(i,0)*smrtrnodes[i]->MoData().n()[1];
      dmap_nysl_gp_mod[p->first] += valy*(p->second);
    }
  }

  double sxsx = sgpnmod[0]*sgpnmod[0];
  double sxsy = sgpnmod[0]*sgpnmod[1];
  double sysy = sgpnmod[1]*sgpnmod[1];

  for (CI p=dmap_nxsl_gp_mod.begin();p!=dmap_nxsl_gp_mod.end();++p)
  {
    dmap_nxsl_gp[p->first] += 1/length*(p->second);
    dmap_nxsl_gp[p->first] -= 1/(length*length*length)*sxsx*(p->second);
    dmap_nysl_gp[p->first] -= 1/(length*length*length)*sxsy*(p->second);
  }

  for (CI p=dmap_nysl_gp_mod.begin();p!=dmap_nysl_gp_mod.end();++p)
  {
    dmap_nysl_gp[p->first] += 1/length*(p->second);
    dmap_nysl_gp[p->first] -= 1/(length*length*length)*sysy*(p->second);
    dmap_nxsl_gp[p->first] -= 1/(length*length*length)*sxsy*(p->second);
  }

  // *********************************************************************
  // finally compute Lin(XiGP_master)
  // *********************************************************************

  // add derivative of slave GP coordinates
  for (CI p=dmap_xsl_gp.begin();p!=dmap_xsl_gp.end();++p)
    derivmxi[p->first] -= sgpn[1]*(p->second);
  for (CI p=dmap_ysl_gp.begin();p!=dmap_ysl_gp.end();++p)
    derivmxi[p->first] += sgpn[0]*(p->second);

  // add derivatives of master node coordinates
  for (int i=0;i<nummnode;++i)
  {
    derivmxi[mmrtrnodes[i]->Dofs()[0]] += valmxigp[i]*sgpn[1];
    derivmxi[mmrtrnodes[i]->Dofs()[1]] -= valmxigp[i]*sgpn[0];
  }

  // add derivative of slave GP normal
  for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
    derivmxi[p->first] -= fac_ymsl_gp*(p->second);
  for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
      derivmxi[p->first] += fac_xmsl_gp*(p->second);

  // multiply all entries with cmxigp
  for (CI p=derivmxi.begin();p!=derivmxi.end();++p)
    derivmxi[p->first] = cmxigp*(p->second);
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiGP master (3D)        popp 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP3D(MORTAR::MortarElement& sele,
                                      MORTAR::MortarElement& mele,
                                      double* sxigp, double* mxigp,
                                      const std::vector<std::map<int,double> >& derivsxi,
                                      std::vector<std::map<int,double> >& derivmxi,
                                      double& alpha)
{
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  std::vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  std::vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  for (int i=0;i<numsnode;++i)
  {
    smrtrnodes[i] = static_cast<MORTAR::MortarNode*>(snodes[i]);
    if (!smrtrnodes[i]) dserror("ERROR: DerivXiGP3D: Null pointer!");
  }

  for (int i=0;i<nummnode;++i)
  {
    mmrtrnodes[i] = static_cast<MORTAR::MortarNode*>(mnodes[i]);
    if (!mmrtrnodes[i]) dserror("ERROR: DerivXiGP3D: Null pointer!");
  }

  // we also need shape function derivs at the GP
  LINALG::SerialDenseVector valsxigp(numsnode);
  LINALG::SerialDenseVector valmxigp(nummnode);
  LINALG::SerialDenseMatrix derivsxigp(numsnode,2,true);
  LINALG::SerialDenseMatrix derivmxigp(nummnode,2,true);

  sele.EvaluateShape(sxigp,valsxigp,derivsxigp,numsnode);
  mele.EvaluateShape(mxigp,valmxigp,derivmxigp,nummnode);

  // we also need the GP slave coordinates + normal
  double sgpn[3] = {0.0,0.0,0.0};
  double sgpx[3] = {0.0,0.0,0.0};
  for (int i=0;i<numsnode;++i)
    for (int k=0;k<3;++k)
    {
      sgpn[k]+=valsxigp[i]*smrtrnodes[i]->MoData().n()[k];
      sgpx[k]+=valsxigp[i]*smrtrnodes[i]->xspatial()[k];
    }

  // build 3x3 factor matrix L
  LINALG::Matrix<3,3> lmatrix(true);
  for (int k=0;k<3;++k) lmatrix(k,2) = -sgpn[k];
  for (int z=0;z<nummnode;++z)
    for (int k=0;k<3;++k)
    {
      lmatrix(k,0) += derivmxigp(z,0) * mmrtrnodes[z]->xspatial()[k];
      lmatrix(k,1) += derivmxigp(z,1) * mmrtrnodes[z]->xspatial()[k];
    }

  // get inverse of the 3x3 matrix L (in place)
  lmatrix.Invert();

  // build directional derivative of slave GP normal
  typedef std::map<int,double>::const_iterator CI;
  std::map<int,double> dmap_nxsl_gp;
  std::map<int,double> dmap_nysl_gp;
  std::map<int,double> dmap_nzsl_gp;

  for (int i=0;i<numsnode;++i)
  {
    std::map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
    std::map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
    std::map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

    for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
      dmap_nxsl_gp[p->first] += valsxigp[i]*(p->second);
    for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
      dmap_nysl_gp[p->first] += valsxigp[i]*(p->second);
    for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
      dmap_nzsl_gp[p->first] += valsxigp[i]*(p->second);

    for (CI p=derivsxi[0].begin();p!=derivsxi[0].end();++p)
    {
      double valx =  derivsxigp(i,0)*smrtrnodes[i]->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx*(p->second);
      double valy =  derivsxigp(i,0)*smrtrnodes[i]->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy*(p->second);
      double valz =  derivsxigp(i,0)*smrtrnodes[i]->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz*(p->second);
    }

    for (CI p=derivsxi[1].begin();p!=derivsxi[1].end();++p)
    {
      double valx =  derivsxigp(i,1)*smrtrnodes[i]->MoData().n()[0];
      dmap_nxsl_gp[p->first] += valx*(p->second);
      double valy =  derivsxigp(i,1)*smrtrnodes[i]->MoData().n()[1];
      dmap_nysl_gp[p->first] += valy*(p->second);
      double valz =  derivsxigp(i,1)*smrtrnodes[i]->MoData().n()[2];
      dmap_nzsl_gp[p->first] += valz*(p->second);
    }
  }

  // start to fill linearization maps for master GP
  // (1) all master nodes coordinates part
  for (int z=0;z<nummnode;++z)
  {
    for (int k=0;k<3;++k)
    {
      derivmxi[0][mmrtrnodes[z]->Dofs()[k]] -= valmxigp[z] * lmatrix(0,k);
      derivmxi[1][mmrtrnodes[z]->Dofs()[k]] -= valmxigp[z] * lmatrix(1,k);
    }
  }

  // (2) slave Gauss point coordinates part
  for (int z=0;z<numsnode;++z)
  {
    for (int k=0;k<3;++k)
    {
      derivmxi[0][smrtrnodes[z]->Dofs()[k]] += valsxigp[z] * lmatrix(0,k);
      derivmxi[1][smrtrnodes[z]->Dofs()[k]] += valsxigp[z] * lmatrix(1,k);

      for (CI p=derivsxi[0].begin();p!=derivsxi[0].end();++p)
      {
        derivmxi[0][p->first] += derivsxigp(z,0) * smrtrnodes[z]->xspatial()[k] * lmatrix(0,k) * (p->second);
        derivmxi[1][p->first] += derivsxigp(z,0) * smrtrnodes[z]->xspatial()[k] * lmatrix(1,k) * (p->second);
      }

      for (CI p=derivsxi[1].begin();p!=derivsxi[1].end();++p)
      {
        derivmxi[0][p->first] += derivsxigp(z,1) * smrtrnodes[z]->xspatial()[k] *lmatrix(0,k) * (p->second);
        derivmxi[1][p->first] += derivsxigp(z,1) * smrtrnodes[z]->xspatial()[k] *lmatrix(1,k) * (p->second);
      }
    }
  }

  // (3) slave Gauss point normal part
  for (CI p=dmap_nxsl_gp.begin();p!=dmap_nxsl_gp.end();++p)
  {
    derivmxi[0][p->first] += alpha * lmatrix(0,0) *(p->second);
    derivmxi[1][p->first] += alpha * lmatrix(1,0) *(p->second);
  }
  for (CI p=dmap_nysl_gp.begin();p!=dmap_nysl_gp.end();++p)
  {
    derivmxi[0][p->first] += alpha * lmatrix(0,1) *(p->second);
    derivmxi[1][p->first] += alpha * lmatrix(1,1) *(p->second);
  }
  for (CI p=dmap_nzsl_gp.begin();p!=dmap_nzsl_gp.end();++p)
  {
    derivmxi[0][p->first] += alpha * lmatrix(0,2) *(p->second);
    derivmxi[1][p->first] += alpha * lmatrix(1,2) *(p->second);
  }

  /*
  // check linearization
  typedef std::map<int,double>::const_iterator CI;
  std::cout << "\nLinearization of current master GP:" << endl;
  std::cout << "-> Coordinate 1:" << endl;
  for (CI p=derivmxi[0].begin();p!=derivmxi[0].end();++p)
    std::cout << p->first << " " << p->second << endl;
  std::cout << "-> Coordinate 2:" << endl;
  for (CI p=derivmxi[1].begin();p!=derivmxi[1].end();++p)
      std::cout << p->first << " " << p->second << endl;
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute deriv. of XiGP slave / master AuxPlane (3D)       popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP3DAuxPlane(MORTAR::MortarElement& ele,
                                        double* xigp, double* auxn,
                                        std::vector<std::map<int,double> >& derivxi, double& alpha,
                                        const std::vector<std::map<int,double> >& derivauxn,
                                        const std::vector<std::map<int,double> >& derivgp)
{
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating element nodes
  DRT::Node** nodes = ele.Nodes();
  std::vector<MORTAR::MortarNode*> mrtrnodes(ele.NumNode());
  int numnode = ele.NumNode();

  for (int i=0;i<numnode;++i)
  {
    mrtrnodes[i] = static_cast<MORTAR::MortarNode*>(nodes[i]);
    if (!mrtrnodes[i]) dserror("ERROR: DerivXiGP3DAuxPlane: Null pointer!");
  }

  // we also need shape function derivs at the GP
  LINALG::SerialDenseVector valxigp(numnode);
  LINALG::SerialDenseMatrix derivxigp(numnode,2,true);
  ele.EvaluateShape(xigp,valxigp,derivxigp,numnode);

  // build 3x3 factor matrix L
  LINALG::Matrix<3,3> lmatrix(true);
  for (int k=0;k<3;++k) lmatrix(k,2) = -auxn[k];
  for (int z=0;z<numnode;++z)
    for (int k=0;k<3;++k)
    {
      lmatrix(k,0) += derivxigp(z,0) * mrtrnodes[z]->xspatial()[k];
      lmatrix(k,1) += derivxigp(z,1) * mrtrnodes[z]->xspatial()[k];
    }

  // get inverse of the 3x3 matrix L (in place)
  lmatrix.Invert();

  // start to fill linearization maps for element GP
  typedef std::map<int,double>::const_iterator CI;

  // (1) all nodes coordinates part
  for (int z=0;z<numnode;++z)
  {
    for (int k=0;k<3;++k)
    {
      derivxi[0][mrtrnodes[z]->Dofs()[k]] -= valxigp[z] * lmatrix(0,k);
      derivxi[1][mrtrnodes[z]->Dofs()[k]] -= valxigp[z] * lmatrix(1,k);
    }
  }

  // (2) Gauss point coordinates part
  for (CI p=derivgp[0].begin();p!=derivgp[0].end();++p)
  {
    derivxi[0][p->first] += lmatrix(0,0) * (p->second);
    derivxi[1][p->first] += lmatrix(1,0) * (p->second);
  }
  for (CI p=derivgp[1].begin();p!=derivgp[1].end();++p)
  {
    derivxi[0][p->first] += lmatrix(0,1) * (p->second);
    derivxi[1][p->first] += lmatrix(1,1) * (p->second);
  }
  for (CI p=derivgp[2].begin();p!=derivgp[2].end();++p)
  {
    derivxi[0][p->first] += lmatrix(0,2) * (p->second);
    derivxi[1][p->first] += lmatrix(1,2) * (p->second);
  }

  // (3) AuxPlane normal part
  for (CI p=derivauxn[0].begin();p!=derivauxn[0].end();++p)
  {
    derivxi[0][p->first] += alpha * lmatrix(0,0) * (p->second);
    derivxi[1][p->first] += alpha * lmatrix(1,0) * (p->second);
  }
  for (CI p=derivauxn[1].begin();p!=derivauxn[1].end();++p)
  {
    derivxi[0][p->first] += alpha * lmatrix(0,1) * (p->second);
    derivxi[1][p->first] += alpha * lmatrix(1,1) * (p->second);
  }
  for (CI p=derivauxn[2].begin();p!=derivauxn[2].end();++p)
  {
    derivxi[0][p->first] += alpha * lmatrix(0,2) * (p->second);
    derivxi[1][p->first] += alpha * lmatrix(1,2) * (p->second);
  }

  /*
  // check linearization
  typedef std::map<int,double>::const_iterator CI;
  std::cout << "\nLinearization of current slave / master GP:" << endl;
  std::cout << "-> Coordinate 1:" << endl;
  for (CI p=derivxi[0].begin();p!=derivxi[0].end();++p)
    std::cout << p->first << " " << p->second << endl;
  std::cout << "-> Coordinate 2:" << endl;
  for (CI p=derivxi[1].begin();p!=derivxi[1].end();++p)
    std::cout << p->first << " " << p->second << endl;
  std::cout << "-> Coordinate 3:" << endl;
  for (CI p=derivxi[2].begin();p!=derivxi[2].end();++p)
    std::cout << p->first << " " << p->second << endl;
  */

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave             |
 |  element to the D map of the adjacent slave nodes.                   |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleD(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    Epetra_SerialDenseMatrix& dseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleD: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
    int sndof = snode->NumDof();
    
    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all slave nodes again ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node again ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = dseg(slave*sndof+sdof,master*mndof+mdof);

          // BOUNDARY NODE MODIFICATION **********************************
          // We have modified their neighbors' dual shape functions, so we
          // now have a problem with off-diagonal entries occuring in D.
          // Of course we want to keep the diagonality property of the D
          // matrix, but still we may not modify the whole Mortar coupling
          // setting! We achieve both by appling a quite simple but very
          // effective trick: The boundary nodes have already been defined
          // as being master nodes, so all we have to do here, is to shift
          // the off-diagonal terms from D to the resepective place in M,
          // which is not diagonal anyway! (Mind the MINUS sign!!!)
          // *************************************************************
          if (mnode->IsOnBound())
          {
            double minusval = -val;
            if(abs(val)>1e-12) snode->AddMValue(sdof,col,minusval);   
            if(abs(val)>1e-12) snode->AddMNode(mnode->Id()); // only for friction!
          }
          else
          { 
            if(abs(val)>1e-12) snode->AddDValue(sdof,col,val);
            if(abs(val)>1e-12) snode->AddSNode(mnode->Id()); // only for friction!
          }  
        }
      }
    }
    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    std::map<int, double> nodemap0 = (snode->GetD())[0];
    std::map<int, double> nodemap1 = (snode->GetD())[1];
    typedef std::map<int,double>::const_iterator CI;

    std::cout << "Row dof id: " << sdofs[0] << endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      std::cout << p->first << '\t' << p->second << endl;

    std::cout << "Row dof id: " << sdofs[1] << endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      std::cout << p->first << '\t' << p->second << endl;
#endif // #ifdef DEBUG
    */
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleD(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::IntElement& sintele,
                                    Epetra_SerialDenseMatrix& dseg)
{
  // get adjacent int nodes to assemble to
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: AssembleD: Null pointer for sintnodes!");
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleD: Null pointer for snodes!");

  // loop over all slave int nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* sintnode = static_cast<CONTACT::CoNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all slave nodes ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = dseg(slave*sintndof+sdof,master*mndof+mdof);

          // assembly
          if(abs(val)>1e-12) sintnode->AddDValue(sdof,col,val);
          if(abs(val)>1e-12) sintnode->AddSNode(mnode->Id()); // only for friction!
        }
      }
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave and master  |
 |  overlap pair to the M map of the adjacent slave nodes.              |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleM(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::MortarElement& mele,
                                    Epetra_SerialDenseMatrix& mseg)
{
  // get adjacent slave nodes and master nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleM: Null pointer for snodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleM: Null pointer for mnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
    int sndof = snode->NumDof();

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sndof+sdof,master*mndof+mdof);
          if(abs(val)>1e-12) snode->AddMValue(sdof,col,val);
          if(abs(val)>1e-12) snode->AddMNode(mnode->Id());  // only for friction!
        }         
      }
    }
    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    std::map<int, double> nodemap0 = (snode->GetM())[0];
    std::map<int, double> nodemap1 = (snode->GetM())[1];
    typedef std::map<int,double>::const_iterator CI;

    std::cout << "Row dof id: " << sdofs[0] << endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      std::cout << p->first << '\t' << p->second << endl;

    std::cout << "Row dof id: " << sdofs[1] << endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      std::cout << p->first << '\t' << p->second << endl;
#endif // #ifdef DEBUG
     */
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Assemble M contribution for FastIntegration(2D / 3D)     farah 01/13|
 |  Assemblle also the contribution for friction                        |
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleM_Fast(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         std::vector<MORTAR::MortarElement*> meles,
                                         Epetra_SerialDenseMatrix& mseg)
{
  //problem dimension via numdof of an arbitrary slave node
  DRT::Node** testnodes = sele.Nodes();
  CONTACT::CoNode* tnode = static_cast<CONTACT::CoNode*>(testnodes[0]);
  int dim = tnode->NumDof();

  for ( int nummasterele=0;nummasterele<(int)meles.size();++nummasterele )
  {
    //if all entries zero-valued --> no assembly
    bool nozero=false;
    for ( int i=0;i<dim*meles[nummasterele]->NumNode();++i  )
    {
      for (int j=0;j<dim*sele.NumNode();++j)
      {
        if(mseg(j,i+nummasterele*meles[nummasterele]->NumNode()*dim) != 0.0 )
        {
          nozero=true;
          break;
        }
      }
    }

    if (nozero==true)
    {
      // get adjacent slave nodes and master nodes
      DRT::Node** snodes = sele.Nodes();
      if (!snodes)
      dserror("ERROR: AssembleM: Null pointer for snodes!");
      DRT::Node** mnodes = meles[nummasterele]->Nodes();
      if (!mnodes)
      dserror("ERROR: AssembleM: Null pointer for mnodes!");

      // loop over all slave nodes
      for (int slave=0;slave<sele.NumNode();++slave)
      {
        CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
        int sndof = snode->NumDof();

        // only process slave node rows that belong to this proc
        if (snode->Owner() != comm.MyPID())
          continue;

        // do not process slave side boundary nodes
        // (their row entries would be zero anyway!)
        if (snode->IsOnBound())
          continue;

        // loop over all dofs of the slave node
        for (int sdof=0;sdof<sndof;++sdof)
        {
          // loop over all master nodes
          for (int master=0;master<meles[nummasterele]->NumNode();++master)
          {
            CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
            const int* mdofs = mnode->Dofs(); //global dofs
            int mndof = mnode->NumDof();

            // loop over all dofs of the master node
            for (int mdof=0;mdof<mndof;++mdof)
            {
              int col = mdofs[mdof];

              double val = mseg(slave*sndof+sdof,master*mndof+mdof+nummasterele*mndof*meles[nummasterele]->NumNode());
              //snode->AddMValue(sdof,col,val);

              if(abs(val)>1e-12) snode->AddMValue(sdof,col,val);
              if(abs(val)>1e-12) snode->AddMNode(mnode->Id());  // only for friction!
            }
          }
        }
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution (2D / 3D)                         popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleM(const Epetra_Comm& comm,
                                    MORTAR::IntElement& sintele,
                                    MORTAR::MortarElement& mele,
                                    Epetra_SerialDenseMatrix& mseg)
{
  // get adjacent slave int nodes and master nodes
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: AssembleM: Null pointer for sintnodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: AssembleM: Null pointer for mnodes!");

  // loop over all slave int nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* sintnode = static_cast<CONTACT::CoNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sintndof+sdof,master*mndof+mdof);
          if(abs(val)>1e-12) sintnode->AddMValue(sdof,col,val);
          if(abs(val)>1e-12) sintnode->AddMNode(mnode->Id()); // only for friction!
        }
      }
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 01/08|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleG(const Epetra_Comm& comm,
                                      MORTAR::MortarElement& sele,
                                      Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = gseg(slave);
    snode->AddgValue(val);
    
    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    std::cout << "Weighted gap: " << snode->Getg() << endl;
#endif // #ifdef DEBUG
    */
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleG(const Epetra_Comm& comm,
                                      MORTAR::IntElement& sintele,
                                      Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave int nodes to assemble to
  DRT::Node** snodes = sintele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for sintnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = gseg(slave);
    snode->AddgValue(val);
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble  mechanical dissipation to slave nodes       gitterle 08/10|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the mechanical dissipation of adjacent nodes        |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleMechDissSlave(const Epetra_Comm& comm,
                                                  MORTAR::MortarElement& sele,
                                                  Epetra_SerialDenseVector& mdissseg)
{
  // get adjacent slave node to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleMechDissSlave: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = mdissseg(slave);
    snode->AddMechDissValue(val);
   }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble  mechanical dissipation to master nodes      gitterle 10/10|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the mechanical dissipation of adjacent nodes        |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleMechDissMaster(const Epetra_Comm& comm,
                                                   MORTAR::MortarElement& mele,
                                                   Epetra_SerialDenseVector& mdissseg)
{
  // get adjacent master node to assemble to
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: AssembleMechDissMaster: Null pointer for snodes!");

  // loop over all master nodes
  for (int master=0;master<mele.NumNode();++master)
  {
    CONTACT::FriNode* mnode = static_cast<CONTACT::FriNode*>(mnodes[master]);

    double val = mdissseg(master);
    mnode->AddMechDissValue(val);
   }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble A contribution (2D / 3D)                     gitterle 10/10|
 |  This method assembles the contrubution of a 1D/2D slave             |
 |  element to the A map of the adjacent slave nodes.                   |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleA(const Epetra_Comm& comm,
                                     MORTAR::MortarElement& sele,
                                     Epetra_SerialDenseMatrix& aseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleA: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
    int sndof = snode->NumDof();
    
    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all slave nodes again ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node again ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = aseg(slave*sndof+sdof,master*mndof+mdof);

          if(abs(val)>1e-12) snode->AddAValue(sdof,col,val);
          if(abs(val)>1e-12) snode->AddANode(mnode->Id()); 
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble B contribution (2D / 3D)                     gitterle 10/10|
 |  This method assembles the contribution of a 1D/2D master            |
 |  element to the B map of the adjacent master node.                   |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleB(const Epetra_Comm& comm,
                                     MORTAR::MortarElement& mele,
                                     Epetra_SerialDenseMatrix& bseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleB: Null pointer for mnodes!");

  // loop over all master nodes
  for (int slave=0;slave<mele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(mnodes[slave]);
    int sndof = snode->NumDof();
    
    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all master nodes again ("master nodes")
      for (int master=0;master<mele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node again ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = bseg(slave*sndof+sdof,master*mndof+mdof);

          if(abs(val)>1e-12) snode->AddBValue(sdof,col,val);
          if(abs(val)>1e-12) snode->AddBNode(mnode->Id()); 
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
  |  Assemble wear contribution (2D / 3D)                  gitterle 12/10|
  |  This method assembles the contribution of a 1D/2D slave and master  |
  |  overlap pair to the wear of the adjacent slave nodes.               |
  *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleWear(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         Epetra_SerialDenseVector& wseg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleWear: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = wseg(slave);
    snode->AddDeltaWearValue(val);
  }

  return true;
}

