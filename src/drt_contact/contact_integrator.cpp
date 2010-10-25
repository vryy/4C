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
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

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
      RCP<Epetra_SerialDenseMatrix> dseg)
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
  typedef map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except line2 (1D) & tri3 (2D)
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if (shapefcn_ == INPAR::MORTAR::shape_dual)
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
    map<int,double> derivjac;
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
          map<int,double>& ddmap_ik = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];
          
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
          map<int,double>& nodemmap = static_cast<CONTACT::CoNode*>(mymrtrnode2)->CoData().GetDerivM()[bgid];

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
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateDerivSegment2D(
     MORTAR::MortarElement& sele, double& sxia, double& sxib,
     MORTAR::MortarElement& mele, double& mxia, double& mxib,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{ 
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

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
  
  // map iterator
  typedef map<int,double>::const_iterator CI;

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
  vector<map<int,double> > ximaps(4);
  DerivXiAB2D(sele,sxia,sxib,mele,mxia,mxib,ximaps,startslave,endslave);

  // prepare directional derivative of dual shape functions
  // this is only necessary for quadratic dual shape functions in 2D
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual) && (sele.Shape()==MORTAR::MortarElement::line3))
  {
    duallin=true;
    sele.DerivShapeDual(dualmap);
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
      cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
      cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
      cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
      dserror("ERROR: IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d",mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(shapefcn_,sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);
    
    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    // decide whether boundary modification has to be considered or not
    // this is element-specific (is there a boundary node in this element?)
    bool bound = false;
    for (int k=0;k<nrow;++k)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
      bound += mymrtrnode->IsOnBound();
    }

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
    else if (shapefcn_ == INPAR::MORTAR::shape_dual)
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

    // build interpolation of slave GP coordinates
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
    
    // evaluate linearizations *******************************************
    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    static_cast<CONTACT::CoElement&>(sele).DJacDXi(djacdxi,sxi,ssecderiv);
    double dxdsxidsxi=djacdxi[0]; // only 2D here

    // evalute the GP slave coordinate derivatives
    map<int,double> dsxigp;
    for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
      dsxigp[p->first] += 0.5*(1-eta[0])*(p->second);
    for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
      dsxigp[p->first] += 0.5*(1+eta[0])*(p->second);

    // evalute the GP master coordinate derivatives
    map<int,double> dmxigp;
    DerivXiGP2D(sele,mele,sxi[0],mxi[0],dsxigp,dmxigp);

    // simple version (no Gauss point projection)
    //for (CI p=ximaps[2].begin();p!=ximaps[2].end();++p)
    //  dmxigp[p->first] += 0.5*(1+eta[0])*(p->second);
    //for (CI p=ximaps[3].begin();p!=ximaps[3].end();++p)
    //  dmxigp[p->first] += 0.5*(1-eta[0])*(p->second);

    // evaluate the Jacobian derivative
    map<int,double> derivjac;
    sele.DerivJacobian(sxi,derivjac);

    // evaluate the GP gap function derivatives
    map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

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
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

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
    map<int,double> dmap_nxsl_gp_unit;
    map<int,double> dmap_nysl_gp_unit;

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
      double prod = lmval[j]*gap;

      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*dxdsxi*dsxideta*wgt;
    }
    // compute segment gap vector ****************************************

    // compute segment D/M linearization *********************************

    // **************** edge modification ********************************
    if (bound)
    {
      // check the shape function type (not really necessary because only dual shape functions arrive here)
      if (shapefcn_ == INPAR::MORTAR::shape_standard)
        dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for standard shape functions");

      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
        bool boundnode = mymrtrnode->IsOnBound();
        int sgid = mymrtrnode->Id();
        map<int,double>& nodemap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];
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
            map<int,double>& nodemmap = static_cast<CONTACT::CoNode*>(mymrtrnode2)->CoData().GetDerivM()[bgid];

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
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*lmderiv(j, 0)*mval[k]*dsxideta*dxdsxi;
          for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*lmval(j, 0)*mderiv(k, 0)*dsxideta*dxdsxi;
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
					map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

					// (1) Lin(Phi) - dual shape functions
					// this vanishes here since there are no deformation-dependent dual functions

					// (2) Lin(NSlave) - slave GP coordinates
					fac = wgt*lmderiv(j, 0)*sval[k]*dsxideta*dxdsxi;
					for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
						ddmap_jk[p->first] += fac*(p->second);

					// (3) Lin(NSlave) - slave GP coordinates
					fac = wgt*lmval(j, 0)*sderiv(k, 0)*dsxideta*dxdsxi;
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
      else if (shapefcn_ == INPAR::MORTAR::shape_dual)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

        // integrate LinM and LinD (NO boundary modification)
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
            if (bound) ddmap_jk[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    }
    // compute segment D/M linearization *********************************

    // compute segment gap linearization *********************************
    for (int j=0;j<nrow;++j)
    {
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
      if (!mymrtrnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

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
    // compute segment gap linearization *********************************
  }
  //**********************************************************************

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
     RCP<MORTAR::IntCell> cell,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
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
  if (cell==null)
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
  typedef map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if ((shapefcn_ == INPAR::MORTAR::shape_dual) && (sele.Shape()!=MORTAR::MortarElement::tri3))
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
        cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
        cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
        cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
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
    cout << "DJacDXi: " << scientific << djacdxi[0] << " " << djacdxi[1] << endl;
    cout << "FD-DJacDXi: " << scientific << fdres[0] << " " << fdres[1] << endl << endl;*/

    // evaluate the slave Jacobian derivative
    map<int,double> jacslavemap;
    sele.DerivJacobian(sxi,jacslavemap);

    // evaluate the intcell Jacobian derivative
    // these are pre-factors for intcell vertex coordinate linearizations
    vector<double> jacintcellvec(2*(cell->NumVertices()));
    cell->DerivJacobian(eta,jacintcellvec);

    // evalute the GP slave coordinate derivatives
    vector<map<int,double> > dsxigp(2);
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
    vector<map<int,double> > dmxigp(2);
    DerivXiGP3D(sele,mele,sxi,mxi,dsxigp,dmxigp,projalpha);

    // evaluate the GP gap function derivatives
    map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

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
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;
    map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
      map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

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
    map<int,double> dmap_nxsl_gp_unit;
    map<int,double> dmap_nysl_gp_unit;
    map<int,double> dmap_nzsl_gp_unit;

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
      double prod = lmval[j]*gap;

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
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];
  
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
					map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

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
      else if (shapefcn_ == INPAR::MORTAR::shape_dual)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jj = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];
  
        // integrate LinM and LinD
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;
  
          // get the correct map as a reference
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];
  
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
      map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

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
      else if( shapefcn_ == INPAR::MORTAR::shape_dual )
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
     RCP<MORTAR::IntCell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg,
     RCP<Epetra_SerialDenseVector> mdisssegs,
     RCP<Epetra_SerialDenseVector> mdisssegm,
     RCP<Epetra_SerialDenseMatrix> aseg,
     RCP<Epetra_SerialDenseMatrix> bseg)
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
  if (cell==null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");
  
  // flag for thermo-structure-interaction with contact
  bool tsi = false;
    if (DRT::Problem::Instance()->ProblemType()=="tsi") tsi=true;

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
  RCP<LINALG::SerialDenseMatrix> scoordold;
  RCP<LINALG::SerialDenseMatrix> mcoordold;
  RCP<LINALG::SerialDenseMatrix> lagmult;
  
  // get them in the case of tsi
  if (tsi)
  {
    scoordold = rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalCoordsOld(*scoordold);
    mcoordold = rcp(new LINALG::SerialDenseMatrix(3,mele.NumNode()));
    mele.GetNodalCoordsOld(*mcoordold);
    lagmult = rcp(new LINALG::SerialDenseMatrix(3,sele.NumNode()));
    sele.GetNodalLagMult(*lagmult);
  }

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // map iterator
  typedef map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if ( (shapefcn_ == INPAR::MORTAR::shape_dual) && (sele.Shape()!=MORTAR::MortarElement::tri3))
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
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
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
    if (tsi)
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

    // mechanical dissipation
    double mechdiss = 0;
    
    if(tsi)
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
        lm(1,0)+=lmval[i]*(*lagmult)(0,i);
        lm(2,0)+=lmval[i]*(*lagmult)(0,i);
      }
       
      // build tangential jump
      Epetra_SerialDenseMatrix jumptan(3,1);
      jumptan.Multiply('N','N',1,tanplane,jump,0.0);
  
      // build tangential lm
      Epetra_SerialDenseMatrix lmtan(3,1);
      lmtan.Multiply('N','N',1,tanplane,lm,0.0);
       
      // build mechanical dissipation
      for (int i=0;i<3;i++)
        mechdiss += lmtan(i,0)*jumptan(i,0);
    }  
    
    // evaluate linearizations *******************************************
    // evaluate the intcell Jacobian derivative
    map<int,double> jacintcellmap;
    cell->DerivJacobian(eta,jacintcellmap);

    // evaluate global GP coordinate derivative
    vector<map<int,double> > lingp(3);
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
    vector<map<int,double> > dsxigp(2);
    DerivXiGP3DAuxPlane(sele,sxi,cell->Auxn(),dsxigp,sprojalpha,cell->GetDerivAuxn(),lingp);

    // evalute the GP master coordinate derivatives
    vector<map<int,double> > dmxigp(2);
    DerivXiGP3DAuxPlane(mele,mxi,cell->Auxn(),dmxigp,mprojalpha,cell->GetDerivAuxn(),lingp);

    // evaluate the GP gap function derivatives
    map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

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
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;
    map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
      map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

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
    map<int,double> dmap_nxsl_gp_unit;
    map<int,double> dmap_nysl_gp_unit;
    map<int,double> dmap_nzsl_gp_unit;

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
      double prod = lmval[j]*gap;
      
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jac*wgt;
    }
    // compute cell gap vector *******************************************
    
    // tsi with contact
    if(tsi)
    {
      // compute cell mechanical dissipation / slave *********************
      // loop over all mdisssegs vector entries
      // nrow represents the slave side dofs !!!
      for (int j=0;j<nrow;++j)
      {
        double prod = lmval[j]*mechdiss;
      
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
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
					map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

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
      else if (shapefcn_ == INPAR::MORTAR::shape_dual)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

        // integrate LinM and LinD
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
      map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

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
     RCP<MORTAR::IntCell> cell, double* auxn,
     INPAR::MORTAR::LagMultQuad3D lmtype,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
  
  /*cout << endl;
  cout << "Slave type: " << sele.Shape() << endl;
  cout << "SlaveElement Nodes:";
  for (int k=0;k<sele.NumNode();++k) cout << " " << sele.NodeIds()[k];
  cout << "\nMaster type: " << mele.Shape() << endl;
  cout << "MasterElement Nodes:";
  for (int k=0;k<mele.NumNode();++k) cout << " " << mele.NodeIds()[k];
  cout << endl;
  cout << "SlaveSub type: " << sintele.Shape() << endl;
  cout << "SlaveSubElement Nodes:";
  for (int k=0;k<sintele.NumNode();++k) cout << " " << sintele.NodeIds()[k];
  cout << "\nMasterSub type: " << mintele.Shape() << endl;
  cout << "MasterSubElement Nodes:";
  for (int k=0;k<mintele.NumNode();++k) cout << " " << mintele.NodeIds()[k];  */

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of slave and master IntElement
  DRT::Element::DiscretizationType sdt = sintele.Shape();
  DRT::Element::DiscretizationType mdt = mintele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell==null)
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
  typedef map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if ( (shapefcn_ == INPAR::MORTAR::shape_dual) && (sele.Shape()!=MORTAR::MortarElement::tri3))
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool dualintlin = false;
  vector<vector<map<int,double> > > dualintmap(nintrow,vector<map<int,double> >(nintrow));
  if ( (shapefcn_ == INPAR::MORTAR::shape_dual) && (sintele.Shape()!=MORTAR::MortarElement::tri3))
  {
    dualintlin = true;
    sintele.DerivShapeDual(dualintmap);
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
  if ( (shapefcn_ == INPAR::MORTAR::shape_dual) && (lmtype == INPAR::MORTAR::lagmult_quad_quad)
    && (sele.Shape() == DRT::Element::quad8 || sele.Shape() == DRT::Element::tri6) )
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
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
        cout << "GP local: " << eta[0] << " " << eta[1] << endl;
        cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << endl;
      }
    }

    // map Gauss point back to slave element (affine map)
    // map Gauss point back to master element (affine map)
    double psxi[2] = {0.0, 0.0};
    double pmxi[2] = {0.0, 0.0};
    sintele.MapToParent(sxi,psxi);
    mintele.MapToParent(mxi,pmxi);

    //cout << "SInt-GP:    " << sxi[0] << " " << sxi[1] << endl;
    //cout << "MInt-GP:    " << mxi[0] << " " << mxi[1] << endl;
    //cout << "SParent-GP: " << psxi[0] << " " << psxi[1] << endl;
    //cout << "MParent-GP: " << pmxi[0] << " " << pmxi[1] << endl;

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
    else if (shapefcn_ == INPAR::MORTAR::shape_dual &&
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
    map<int,double> jacintcellmap;
    cell->DerivJacobian(eta,jacintcellmap);

    // evaluate global GP coordinate derivative
    vector<map<int,double> > lingp(3);
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
    vector<map<int,double> > dsxigp(2);
    DerivXiGP3DAuxPlane(sintele,sxi,cell->Auxn(),dsxigp,sprojalpha,cell->GetDerivAuxn(),lingp);

    // evalute the GP master coordinate derivatives
    vector<map<int,double> > dmxigp(2);
    DerivXiGP3DAuxPlane(mintele,mxi,cell->Auxn(),dmxigp,mprojalpha,cell->GetDerivAuxn(),lingp);

    // map GP coordinate derivatives back to slave element (affine map)
    // map GP coordinate derivatives back to master element (affine map)
    vector<map<int,double> > dpsxigp(2);
    vector<map<int,double> > dpmxigp(2);
    sintele.MapToParent(dsxigp,dpsxigp);
    mintele.MapToParent(dmxigp,dpmxigp);

    // evaluate the GP gap function derivatives
    map<int,double> dgapgp;

    // we need the participating slave and master nodes
    DRT::Node** snodes = sele.Nodes();
    DRT::Node** mnodes = mele.Nodes();
    vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
    vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());

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
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;
    map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
      map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

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
    map<int,double> dmap_nxsl_gp_unit;
    map<int,double> dmap_nysl_gp_unit;
    map<int,double> dmap_nzsl_gp_unit;

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
    else if (shapefcn_ == INPAR::MORTAR::shape_dual &&
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
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
					map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

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
            map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
							map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[sgid];

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
							map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

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
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
					map<int,double>& ddmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

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
    else if (shapefcn_ == INPAR::MORTAR::shape_dual &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
        int sgid = mymrtrnode->Id();

        // for dual shape functions ddmap_jj and dmmap_jk can be calculated together
        map<int,double>& ddmap_jj = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivD()[sgid];

        // integrate LinM and LinD
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivM()[mgid];

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
        map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

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
        map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();
  
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
    else if (shapefcn_ == INPAR::MORTAR::shape_dual &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      for (int j=0;j<nrow;++j)
      {
        MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[j]);
        if (!mymrtrnode) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
        double fac = 0.0;

        // get the corresponding map as a reference
        map<int,double>& dgmap = static_cast<CONTACT::CoNode*>(mymrtrnode)->CoData().GetDerivG();

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
 |  Compute penalty scaling factor kappa                      popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::IntegrateKappaPenalty(MORTAR::MortarElement& sele,
                                                  double* sxia, double* sxib,
                                                  RCP<Epetra_SerialDenseVector> gseg)
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
  typedef map<int,double>::const_iterator CI;

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
                                                  RCP<Epetra_SerialDenseVector> gseg,
                                                  INPAR::MORTAR::LagMultQuad3D lmtype)
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
  typedef map<int,double>::const_iterator CI;

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
                                    vector<map<int,double> >& derivxi,
                                    bool& startslave, bool& endslave)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
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
    //cout << "cmxib: " << cmxib << endl;

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
    //cout << "cmxia: " << cmxia << endl;

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
    //cout << "csxia: " << csxia << endl;
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
    //cout << "csxib: " << csxib << endl;
  }

  // prepare linearizations
  typedef map<int,double>::const_iterator CI;

  // *********************************************************************
  // finally compute Lin(XiAB_master)
  // *********************************************************************
  // build DerivXiBMaster if start node = slave node
  if (startslave==true)
  {
    map<int,double> dmap_mxib;
    map<int,double>& nxmap_b = static_cast<CONTACT::CoNode*>(smrtrnodes[0])->CoData().GetDerivN()[0];
    map<int,double>& nymap_b = static_cast<CONTACT::CoNode*>(smrtrnodes[0])->CoData().GetDerivN()[1];

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
    map<int,double> dmap_mxia;
    map<int,double>& nxmap_a = static_cast<CONTACT::CoNode*>(smrtrnodes[1])->CoData().GetDerivN()[0];
    map<int,double>& nymap_a = static_cast<CONTACT::CoNode*>(smrtrnodes[1])->CoData().GetDerivN()[1];

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
    map<int,double> dmap_sxia;

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
      map<int,double>& nxmap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      map<int,double>& nymap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

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
    map<int,double> dmap_sxib;

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
      map<int,double>& nxmap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
      map<int,double>& nymap_curr = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

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
                                    const map<int,double>& derivsxi,
                                    map<int,double>& derivmxi)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
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
  //cout << "cmxigp: " << cmxigp << endl;

  fac_xmsl_gp -= sgpx[0];
  fac_ymsl_gp -= sgpx[1];

  // prepare linearization
  typedef map<int,double>::const_iterator CI;

  // build directional derivative of slave GP coordinates
  map<int,double> dmap_xsl_gp;
  map<int,double> dmap_ysl_gp;

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
  map<int,double> dmap_nxsl_gp;
  map<int,double> dmap_nysl_gp;

  double sgpnmod[3] = {0.0,0.0,0.0};
  for (int i=0;i<3;++i) sgpnmod[i]=sgpn[i]*length;

  map<int,double> dmap_nxsl_gp_mod;
  map<int,double> dmap_nysl_gp_mod;

  for (int i=0;i<numsnode;++i)
  {
    map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
    map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];

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
                                      const vector<map<int,double> >& derivsxi,
                                      vector<map<int,double> >& derivmxi,
                                      double& alpha)
{
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  vector<MORTAR::MortarNode*> smrtrnodes(sele.NumNode());
  vector<MORTAR::MortarNode*> mmrtrnodes(mele.NumNode());
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
  typedef map<int,double>::const_iterator CI;
  map<int,double> dmap_nxsl_gp;
  map<int,double> dmap_nysl_gp;
  map<int,double> dmap_nzsl_gp;

  for (int i=0;i<numsnode;++i)
  {
    map<int,double>& dmap_nxsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[0];
    map<int,double>& dmap_nysl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[1];
    map<int,double>& dmap_nzsl_i = static_cast<CONTACT::CoNode*>(smrtrnodes[i])->CoData().GetDerivN()[2];

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
  typedef map<int,double>::const_iterator CI;
  cout << "\nLinearization of current master GP:" << endl;
  cout << "-> Coordinate 1:" << endl;
  for (CI p=derivmxi[0].begin();p!=derivmxi[0].end();++p)
    cout << p->first << " " << p->second << endl;
  cout << "-> Coordinate 2:" << endl;
  for (CI p=derivmxi[1].begin();p!=derivmxi[1].end();++p)
      cout << p->first << " " << p->second << endl;
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute deriv. of XiGP slave / master AuxPlane (3D)       popp 03/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegrator::DerivXiGP3DAuxPlane(MORTAR::MortarElement& ele,
                                        double* xigp, double* auxn,
                                        vector<map<int,double> >& derivxi, double& alpha,
                                        const vector<map<int,double> >& derivauxn,
                                        const vector<map<int,double> >& derivgp)
{
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating element nodes
  DRT::Node** nodes = ele.Nodes();
  vector<MORTAR::MortarNode*> mrtrnodes(ele.NumNode());
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
  typedef map<int,double>::const_iterator CI;

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
  typedef map<int,double>::const_iterator CI;
  cout << "\nLinearization of current slave / master GP:" << endl;
  cout << "-> Coordinate 1:" << endl;
  for (CI p=derivxi[0].begin();p!=derivxi[0].end();++p)
    cout << p->first << " " << p->second << endl;
  cout << "-> Coordinate 2:" << endl;
  for (CI p=derivxi[1].begin();p!=derivxi[1].end();++p)
    cout << p->first << " " << p->second << endl;
  cout << "-> Coordinate 3:" << endl;
  for (CI p=derivxi[2].begin();p!=derivxi[2].end();++p)
    cout << p->first << " " << p->second << endl;
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
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
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
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(snodes[master]);
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
    cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    map<int, double> nodemap0 = (snode->GetD())[0];
    map<int, double> nodemap1 = (snode->GetD())[1];
    typedef map<int,double>::const_iterator CI;

    cout << "Row dof id: " << sdofs[0] << endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      cout << p->first << '\t' << p->second << endl;

    cout << "Row dof id: " << sdofs[1] << endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      cout << p->first << '\t' << p->second << endl;
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
    MORTAR::MortarNode* sintnode = static_cast<MORTAR::MortarNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all slave nodes ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(snodes[master]);
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
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
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
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
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
    cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    map<int, double> nodemap0 = (snode->GetM())[0];
    map<int, double> nodemap1 = (snode->GetM())[1];
    typedef map<int,double>::const_iterator CI;

    cout << "Row dof id: " << sdofs[0] << endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      cout << p->first << '\t' << p->second << endl;

    cout << "Row dof id: " << sdofs[1] << endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      cout << p->first << '\t' << p->second << endl;
#endif // #ifdef DEBUG
     */
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
    MORTAR::MortarNode* sintnode = static_cast<MORTAR::MortarNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
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
    cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    cout << "Weighted gap: " << snode->Getg() << endl;
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
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
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
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(snodes[master]);
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
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(mnodes[slave]);
    int sndof = snode->NumDof();
    
    // FIXGIT: not working in parallel yet
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
      // loop over all master nodes again ("master nodes")
      for (int master=0;master<mele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
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

#endif //#ifdef CCADISCRET
