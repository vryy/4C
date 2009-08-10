/*!----------------------------------------------------------------------
\file drt_contact_integrator.cpp
\brief A class to perform integrations of Mortar matrices on the overlap
 			 of two CElements in 1D and 2D

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

#include "drt_contact_integrator.H"
#include "drt_celement.H"
#include "drt_contact_projector.H"
#include "contactdefines.H"
#include "../drt_fem_general/drt_utils_integration.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 08/08|
 *----------------------------------------------------------------------*/
CONTACT::Integrator::Integrator(DRT::Element::DiscretizationType eletype) :
shapefcn_(Interface::Undefined)
{
  InitializeGP(eletype);
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 07/09|
 *----------------------------------------------------------------------*/
CONTACT::Integrator::Integrator(const Interface::ShapeFcnType shapefcn,
                                DRT::Element::DiscretizationType eletype) :
shapefcn_(shapefcn)
{
  InitializeGP(eletype);
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::InitializeGP(DRT::Element::DiscretizationType eletype)
{
  //*********************************************************************
  // Create integration points according to eletype!
  // Note that our standard Gauss rules are:
  // 5 points: for integrals on 1D lines                 (1,2,3,4,5)
  // 7 points: for integrals on 2D triangles             (1,3,6,7,12,37)
  // 9 points: for integrals on 2D quadrilaterals        (1,4,9)
  //**********************************************************************
  switch(eletype)
  {
  case DRT::Element::line2:
  case DRT::Element::line3:
  {
    dim_=2;
    const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_5point);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=0.0;
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    dim_=3;
    const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_7point);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    dim_=3;
    const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_quad_9point);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  default:
    dserror("ERROR: Integrator: This contact element type is not implemented!");
  } // switch(eletype)
  
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
void CONTACT::Integrator::IntegrateDerivSlave2D3D(
      CONTACT::CElement& sele, double* sxia, double* sxib,
      RCP<Epetra_SerialDenseMatrix> dseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == Interface::Undefined)
    dserror("ERROR: IntegrateDerivSlave2D3D called without specific shape function defined!");
    
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: IntegrateDerivSlave2D3D called on a non-slave CElement!");
  if ((sxia[0]<-1.0) || (sxia[1]<-1.0) || (sxib[0]>1.0) || (sxib[1]>1.0))
    dserror("ERROR: IntegrateDerivSlave2D3D called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();
  int ndof = Dim();
  int ncol = nrow;

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,2,true);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);

  // get slave element nodes themselves for GP normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if (!mynodes) dserror("ERROR: IntegrateDerivSlave2D3D: Null pointer!");

  // map iterator
  typedef map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except line2 (1D) & tri3 (2D)
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if (shapefcn_ == Interface::DualFunctions)
  {
    if (sele.Shape()!=CElement::line2 && sele.Shape()!=CElement::tri3)
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

    // evaluate trace space and dual space shape functions
    sele.EvaluateShape(eta,val,deriv,nrow);
    if (shapefcn_ == Interface::DualFunctions)
      sele.EvaluateShapeDual(eta,dualval,dualderiv,nrow);

    // evaluate the Jacobian det
    double dxdsxi = sele.Jacobian(eta);

    // evaluate the Jacobian derivative
    map<int,double> derivjac;
    sele.DerivJacobian(eta,derivjac);

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
        double prod = 0.0;
        if (shapefcn_ == Interface::DualFunctions)
          prod = dualval[jindex]*val[kindex];
        else if (shapefcn_ == Interface::StandardFunctions)
          prod =val[jindex]*val[kindex];

        // isolate the dseg entries to be filled
        // (both the main diagonal and every other secondary diagonal)
        // and add current Gauss point's contribution to dseg
        if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          (*dseg)(j,k) += prod*dxdsxi*wgt;
      }
    }
    // compute element D matrix ******************************************

    // compute element D linearization ***********************************
    
    // here we have to adapt the alogrithm to support arbitrary shape functions
    // because we can't rely on the diagonality of D anymore
    // this was built analogous to building derivM from IntegrateDerivSegment2D
    
    for (int i=0;i<nrow;++i)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[i]);
      if (!mycnode) dserror("ERROR: IntegrateAndDerivSlave: Null pointer!");
      bool bound = mycnode->IsOnBound();

      //******************************************************************
      // standard case (node i reqiures no edge node modification)
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
          map<int,double>& ddmap_ik = mycnode->GetDerivD()[sgid];
          
          // multiply the corresponding two shape functions
          double prod = 0.0;
          if (shapefcn_ == Interface::DualFunctions)
            prod = wgt*val[k]*dualval[i];
          else if (shapefcn_ == Interface::StandardFunctions)
            prod = wgt*val[k]*val[i];

          // derivative of Jacobian
          for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
            ddmap_ik[p->first] += prod*(p->second);
          
          // derivative of dual shape function
          if (duallin)
          {
            for (int j=0;j<nrow;++j)
            {
              double fac = wgt*val[i]*val[j]*dxdsxi;
              for (CI p=dualmap[i][j].begin();p!=dualmap[i][j].end();++p)
                ddmap_ik[p->first] += fac*(p->second);
            }
          }
        } 
      }

      //******************************************************************
      // edge node case (node i reqiures edge node modification)
      //******************************************************************
      // all nodes on edges cause additional entries in derivM
      // this is a special case depending on the bound-flag of cnode
      // at the moment this is only possible for dual shape functions in 2D
      else
      {
        // check the shape function type
        if (shapefcn_ == Interface::StandardFunctions)
          dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for standard shape functions");
                
        //check for problem dimension
        if (Dim()==3)
          dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for 3D");

        // get gid of current boundary node
        int bgid = mycnode->Id();

        // loop over other nodes (interior nodes)
        for (int k=0;k<nrow;++k)
        {
          CONTACT::CNode* mycnode2 = static_cast<CONTACT::CNode*>(mynodes[k]);
          if (!mycnode2) dserror("ERROR: IntegrateAndDerivSlave: Null pointer!");
          bool bound2 = mycnode2->IsOnBound();
          if (bound2) continue;
          map<int,double>& nodemmap = mycnode2->GetDerivM()[bgid];

          // derivative of Jacobian
          double fac = wgt*val[i]*dualval[k];
          for (CI p=derivjac.begin();p!=derivjac.end();++p)
            nodemmap[p->first] -= fac*(p->second);
          
          // derivative of dual shape functions
          if (duallin)
          {
            for (int j=0;j<nrow;++j)
            {
              LINALG::SerialDenseVector vallin(nrow-1);
              LINALG::SerialDenseMatrix derivlin(nrow-1,1);
              if (i==0) sele.ShapeFunctions(CElement::dual1D_base_for_edge0,eta,vallin,derivlin);
              else if (i==1) sele.ShapeFunctions(CElement::dual1D_base_for_edge1,eta,vallin,derivlin);
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
void CONTACT::Integrator::IntegrateDerivSegment2D(
     CONTACT::CElement& sele, double& sxia, double& sxib,
     CONTACT::CElement& mele, double& mxia, double& mxib,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == Interface::Undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of CElement pair!");
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
  LINALG::SerialDenseMatrix ssecderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  mele.GetNodalCoords(mcoord);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

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
  if ((shapefcn_ == Interface::DualFunctions) && (sele.Shape()==CElement::line3))
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

    // coordinate transformation sxi->eta (slave CElement->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    CONTACT::Projector projector(2);
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

    // evaluate dual space shape functions (on slave element)
    if (shapefcn_ == Interface::DualFunctions)
      sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    // evaluate 2nd deriv of trace space shape functions (on slave element)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mycnode->n()[0];
      gpn[1]+=sval[i]*mycnode->n()[1];
      gpn[2]+=sval[i]*mycnode->n()[2];

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

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    // evaluate linearizations *******************************************
    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    sele.DJacDXi(djacdxi,sxi,ssecderiv);
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
    vector<CONTACT::CNode*> scnodes(sele.NumNode());
    vector<CONTACT::CNode*> mcnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
      if (!scnodes[i]) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
      if (!mcnodes[i]) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = scnodes[i]->GetDerivN()[0];
      map<int,double>& dmap_nysl_i = scnodes[i]->GetDerivN()[1];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
      {
        double valx =  sderiv(i,0)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*scnodes[i]->n()[1];
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
        dgapgp[scnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * scnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<2;++k)
      {
        dgapgp[mcnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mcnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************

    // reflect the CONTACTONEMORTARLOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef CONTACTONEMORTARLOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef CONTACTONEMORTARLOOP
    
    // decide whether boundary modification has to be considered or not
    // this is element-specific (is there a boundary node in this element?)
    bool bound = false;
#ifdef CONTACTBOUNDMOD
    for (int k=0;k<nrow;++k)
    {
      // check the shape function type (not really necessary because only dual shape functions arrive here)
      if (shapefcn_ == Interface::StandardFunctions)
        dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for standard shape functions");
      
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[k]);
      if (!mycnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
      bound += mycnode->IsOnBound();
    }
#endif // #ifdef CONTACTBOUNDMOD

    // compute segment D/M matrix ****************************************
    
    if (shapefcn_ == Interface::StandardFunctions)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for standard shape functions we use the same algorithm
        // for dseg as in IntegrateDerivSlave2D3D (but with modified integration area)
        // hence, mseg and dseg can not be combined into one loop

        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = sval[jindex]*mval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to mseg and dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
          }
        }

        if (dod)
        {
          for (int k=0; k<nrow*ndof; ++k)
          {
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);

            // multiply the two shape functions
            double prod = sval[jindex]*sval[kindex];

            // isolate the mseg and dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to mseg and dseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            {
              (*dseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
            }
          }
        }
      } // nrow*ndof loop
    }
    
    else if (shapefcn_ == Interface::DualFunctions)
    { 
      // loop over all mseg matrix entries
      // nrow represents the slave Lagrange multipliers !!!
      // ncol represents the master dofs !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0;j<nrow*ndof;++j)
      {
        // evaluate dseg entries seperately
        // (only for one mortar loop AND boundary modification)
        if (dod && bound)
        {
          CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[(int)(j/ndof)]);
          if (!mycnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
          bool j_boundnode = mycnode->IsOnBound();
          
          for (int k=0;k<nrow*ndof;++k)
          {
            CONTACT::CNode* mycnode2 = static_cast<CONTACT::CNode*>(mynodes[(int)(k/ndof)]);
            if (!mycnode2) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
            bool k_boundnode = mycnode2->IsOnBound();
            
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);
            
            // do not assemble off-diagonal terms if j,k are both non-boundary nodes
            if (!j_boundnode && !k_boundnode && (jindex!=kindex)) continue;
              
            // multiply the two shape functions 
            double prod = dualval[jindex]*sval[kindex];
    
            // isolate the dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to dseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
          }
        }
        
        // evaluate mseg entries
        // (and dseg entries for one mortar loop and NO boundary modification)
        for (int k=0;k<ncol*ndof;++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = dualval[jindex]*mval[kindex];
  
          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal for mseg)
          // (only the main diagonal for dseg)
          // and add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
            if(dod && !bound) (*dseg)(j,j) += prod*dxdsxi*dsxideta*wgt;
          }
        } 
      }
    } // shapefcn_ switch
    // compute segment D/M matrix ****************************************

    // compute segment gap vector ****************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!
    for (int j=0;j<nrow;++j)
    {
      double prod = 0.0;
#ifdef CONTACTPETROVGALERKIN
      if (shapefcn_ == Interface::StandardFunctions)
        dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
      prod = sval[j]*gap;
#else
      if (shapefcn_ == Interface::DualFunctions)
        prod = dualval[j]*gap;
      else if (shapefcn_ == Interface::StandardFunctions)
        prod = sval[j]*gap;
#endif // #ifdef CONTACTPETROVGALERKIN

      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*dxdsxi*dsxideta*wgt;
    }
    // compute segment gap vector ****************************************

    // compute segment D/M linearization *********************************
    // In the case with separate loops for D and M, this is done by the
    // method IntegrateDerivSlave2D3D(). But in the CONTACTONEMORTARLOOP
    // case we want to combine the linearization of D and M, just as we
    // combine the computation of D and M itself.
    
    // **************** edge modification ********************************
    if (dod && bound)
    {
      // check the shape function type (not really necessary because only dual shape functions arrive here)
      if (shapefcn_ == Interface::StandardFunctions)
        dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for standard shape functions");
      
      for (int j=0;j<nrow;++j)
      {
        CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
        if (!mycnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
        bool boundnode = mycnode->IsOnBound();
        int sgid = mycnode->Id();
        map<int,double>& nodemap = mycnode->GetDerivD()[sgid];
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
          fac = wgt*dualderiv(j,0)*sval[j]*dsxideta*dxdsxi;
          for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
            nodemap[p->first] += fac*(p->second);
  
          // (3) Lin(NSlave) - slave GP coordinates
          fac = wgt*dualval[j]*sderiv(j,0)*dsxideta*dxdsxi;
          for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
            nodemap[p->first] += fac*(p->second);
  
          // (4) Lin(dsxideta) - segment end coordinates
          fac = wgt*dualval[j]*sval[j]*dxdsxi;
          for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
            nodemap[p->first] -= 0.5*fac*(p->second);   
          for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
            nodemap[p->first] += 0.5*fac*(p->second);
  
          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*dualval[j]*sval[j]*dsxideta;
          for (CI p=derivjac.begin();p!=derivjac.end();++p)
            nodemap[p->first] += fac*(p->second);
  
          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*dualval[j]*sval[j]*dsxideta*dxdsxidsxi;
          for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
            nodemap[p->first] += fac*(p->second);
        }

        //******************************************************************
        // edge case (node j is a boundary node)
        //******************************************************************
        else
        {
          // get gid of current boundary node
          int bgid = mycnode->Id();

          // loop over other nodes (interior nodes)
          for (int k=0;k<nrow;++k)
          {
            CONTACT::CNode* mycnode2 = static_cast<CONTACT::CNode*>(mynodes[k]);
            if (!mycnode2) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
            bool boundnode2 = mycnode2->IsOnBound();
            if (boundnode2) continue;
            map<int,double>& nodemmap = mycnode2->GetDerivM()[bgid];

            // (1) Lin(Phi) - dual shape functions
            if (duallin)
              for (int m=0;m<nrow;++m)
              {
                LINALG::SerialDenseVector vallin(nrow-1);
                LINALG::SerialDenseMatrix derivlin(nrow-1,1);
                if (j==0) sele.ShapeFunctions(CElement::dual1D_base_for_edge0,sxi,vallin,derivlin);
                else if (j==1) sele.ShapeFunctions(CElement::dual1D_base_for_edge1,sxi,vallin,derivlin);
                double fac = wgt*sval[j]*vallin[m]*dsxideta*dxdsxi;
                for (CI p=dualmap[k][m].begin();p!=dualmap[k][m].end();++p)
                  nodemmap[p->first] -= fac*(p->second);
              }
    
            // (2) Lin(Phi) - slave GP coordinates
            fac = wgt*dualderiv(k,0)*sval[j]*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
              nodemmap[p->first] -= fac*(p->second);
    
            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*dualval[k]*sderiv(j,0)*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
              nodemmap[p->first] -= fac*(p->second);
    
            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*dualval[k]*sval[j]*dxdsxi;
            for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
              nodemmap[p->first] += 0.5*fac*(p->second);   
            for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
              nodemmap[p->first] -= 0.5*fac*(p->second);
    
            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*dualval[k]*sval[j]*dsxideta;
            for (CI p=derivjac.begin();p!=derivjac.end();++p)
              nodemmap[p->first] -= fac*(p->second);
    
            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*dualval[k]*sval[j]*dsxideta*dxdsxidsxi;
            for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
              nodemmap[p->first] -= fac*(p->second);
          }
        }
      }
    }
    
    // **************** no edge modification *****************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

      int sgid = mycnode->Id();

      // for standard shape functions we use the same algorithm
      // for ddmap_jk as in IntegrateDerivSlave2D3D 
      // this means that ddmap_jk and dmmap_jk have to be calculated separately
      if (shapefcn_ == Interface::StandardFunctions)
      {
        // loop over master nodes for building nodal derivM 
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions 

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*sderiv(j, 0)*mval[k]*dsxideta*dxdsxi;
          for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*sval(j, 0)*mderiv(k, 0)*dsxideta*dxdsxi;
          for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - segment end coordinates
          fac = wgt*sval[j]*mval[k]*dxdsxi;
          for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
            dmmap_jk[p->first] -= 0.5*fac*(p->second);
          for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
            dmmap_jk[p->first] += 0.5*fac*(p->second);

          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*sval[j]*mval[k]*dsxideta;
          for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*sval[j]*mval[k]*dsxideta*dxdsxidsxi;
          for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes

        if (dod)
        {
          // loop over slave nodes for building nodal derivD
          for (int k=0; k<nrow; ++k)
          {
            // global slave node ID
            int sgid = sele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            map<int,double>& ddmap_jk = mycnode->GetDerivD()[sgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*sderiv(j, 0)*sval[k]*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*sval(j, 0)*sderiv(k, 0)*dsxideta*dxdsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - segment end coordinates
            fac = wgt*sval[j]*sval[k]*dxdsxi;
            for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
              ddmap_jk[p->first] -= 0.5*fac*(p->second);
            for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
              ddmap_jk[p->first] += 0.5*fac*(p->second);

            // (5) Lin(dxdsxi) - slave GP Jacobian
            fac = wgt*sval[j]*sval[k]*dsxideta;
            for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*sval[j]*sval[k]*dsxideta*dxdsxidsxi;
            for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
          } // loop over slave nodes
        }
      }
      // for dual shape functions we can make use
      // of the row summing lemma: D_jj = Sum(k) M_jk which applies also for the linearizations
      // this enables us to compute ddmap_jk and dmmap_jk in parallel
      else if (shapefcn_ == Interface::DualFunctions)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jk = mycnode->GetDerivD()[sgid];

        // loop over master nodes
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          for (int m=0; m<nrow; ++m)
          {
            fac = wgt*sval[m]*mval[k]*dsxideta*dxdsxi;
            for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
            {
              dmmap_jk[p->first] += fac*(p->second);
              if (dod && !bound) ddmap_jk[p->first] += fac*(p->second);
            }
          }

          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*dualderiv(j, 0)*mval[k]*dsxideta*dxdsxi;

          for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod && !bound) ddmap_jk[p->first] += fac*(p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*dualval(j, 0)*mderiv(k, 0)*dsxideta*dxdsxi;

          for (CI p=dmxigp.begin(); p!=dmxigp.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod && !bound) ddmap_jk[p->first] += fac*(p->second);
          }

          // (4) Lin(dsxideta) - segment end coordinates
          fac = wgt*dualval[j]*mval[k]*dxdsxi;

          for (CI p=ximaps[0].begin(); p!=ximaps[0].end(); ++p)
          {
            dmmap_jk[p->first] -= 0.5*fac*(p->second);
            if (dod && !bound) ddmap_jk[p->first] -= 0.5*fac*(p->second);
          }
          for (CI p=ximaps[1].begin(); p!=ximaps[1].end(); ++p)
          {
            dmmap_jk[p->first] += 0.5*fac*(p->second);
            if (dod && !bound) ddmap_jk[p->first] += 0.5*fac*(p->second);
          }

          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*dualval[j]*mval[k]*dsxideta;

          for (CI p=derivjac.begin(); p!=derivjac.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod && !bound) ddmap_jk[p->first] += fac*(p->second);
          }

          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*dualval[j]*mval[k]*dsxideta*dxdsxidsxi;

          for (CI p=dsxigp.begin(); p!=dsxigp.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod && !bound) ddmap_jk[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    }
    // compute segment D/M linearization *********************************

    // compute segment gap linearization *********************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      map<int,double>& dgmap = mycnode->GetDerivG();

#ifdef CONTACTPETROVGALERKIN
      
      if (shapefcn_ == Interface::StandardFunctions)
        dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
      
      // (2) Lin(N) - slave GP coordinates
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
#else
      // (1) Lin(Phi) - dual shape functions
      if (shapefcn_ == Interface::DualFunctions)
      {
        for (int m=0;m<nrow;++m)
        {
          fac = wgt*sval[m]*gap*dsxideta*dxdsxi;
          for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
            dgmap[p->first] += fac*(p->second);
        }
      }

      // (2) Lin(Phi) - slave GP coordinates
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualderiv(j,0)*gap*dsxideta*dxdsxi;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sderiv(j,0)*gap*dsxideta*dxdsxi;
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (3) Lin(g) - gap function
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*dsxideta*dxdsxi;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*dsxideta*dxdsxi;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - segment end coordinates
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*gap*dxdsxi;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*gap*dxdsxi;
      for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
        dgmap[p->first] -= 0.5*fac*(p->second);
      for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
        dgmap[p->first] += 0.5*fac*(p->second);

      // (5) Lin(dxdsxi) - slave GP Jacobian
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*gap*dsxideta;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*gap*dsxideta;
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (6) Lin(dxdsxi) - slave GP coordinates
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*gap*dsxideta*dxdsxidsxi;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*gap*dsxideta*dxdsxidsxi;
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        dgmap[p->first] += fac*(p->second);
#endif  // #ifdef CONTACTPETROVGALERKIN

    }
    // compute segment gap linearization *********************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiAB (2D)               popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivXiAB2D(CONTACT::CElement& sele,
                                    double& sxia, double& sxib,
                                    CONTACT::CElement& mele,
                                    double& mxia, double& mxib,
                                    vector<map<int,double> >& derivxi,
                                    bool& startslave, bool& endslave)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  vector<CONTACT::CNode*> scnodes(sele.NumNode());
  vector<CONTACT::CNode*> mcnodes(mele.NumNode());
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  for (int i=0;i<numsnode;++i)
  {
    scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
    if (!scnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  for (int i=0;i<nummnode;++i)
  {
    mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
    if (!mcnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
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
      fac_dxm_b += derivmxib(i,0)*(mcnodes[i]->xspatial()[0]);
      fac_dym_b += derivmxib(i,0)*(mcnodes[i]->xspatial()[1]);
      fac_xmsl_b += valmxib[i]*(mcnodes[i]->xspatial()[0]);
      fac_ymsl_b += valmxib[i]*(mcnodes[i]->xspatial()[1]);
    }

    cmxib = -1/(fac_dxm_b*(scnodes[0]->n()[1])-fac_dym_b*(scnodes[0]->n()[0]));
    //cout << "cmxib: " << cmxib << endl;

    fac_xmsl_b -= scnodes[0]->xspatial()[0];
    fac_ymsl_b -= scnodes[0]->xspatial()[1];
  }

  // compute leading constant for DerivXiAMaster if end node = slave node
  if (endslave==true)
  {
    for (int i=0;i<nummnode;++i)
    {
      fac_dxm_a += derivmxia(i,0)*(mcnodes[i]->xspatial()[0]);
      fac_dym_a += derivmxia(i,0)*(mcnodes[i]->xspatial()[1]);
      fac_xmsl_a += valmxia[i]*(mcnodes[i]->xspatial()[0]);
      fac_ymsl_a += valmxia[i]*(mcnodes[i]->xspatial()[1]);
    }

    cmxia = -1/(fac_dxm_a*(scnodes[1]->n()[1])-fac_dym_a*(scnodes[1]->n()[0]));
    //cout << "cmxia: " << cmxia << endl;

    fac_xmsl_a -= scnodes[1]->xspatial()[0];
    fac_ymsl_a -= scnodes[1]->xspatial()[1];
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
      fac_dxsl_a += derivsxia(i,0)*(scnodes[i]->xspatial()[0]);
      fac_dysl_a += derivsxia(i,0)*(scnodes[i]->xspatial()[1]);
      fac_xslm_a += valsxia[i]*(scnodes[i]->xspatial()[0]);
      fac_yslm_a += valsxia[i]*(scnodes[i]->xspatial()[1]);
      fac_dnx_a  += derivsxia(i,0)*(scnodes[i]->n()[0]);
      fac_dny_a  += derivsxia(i,0)*(scnodes[i]->n()[1]);
      fac_nx_a   += valsxia[i]*(scnodes[i]->n()[0]);
      fac_ny_a   += valsxia[i]*(scnodes[i]->n()[1]);
    }

    fac_xslm_a -= mcnodes[1]->xspatial()[0];
    fac_yslm_a -= mcnodes[1]->xspatial()[1];

    csxia = -1/(fac_dxsl_a*fac_ny_a - fac_dysl_a*fac_nx_a + fac_xslm_a*fac_dny_a - fac_yslm_a*fac_dnx_a);
    //cout << "csxia: " << csxia << endl;
  }

  // compute leading constant for DerivXiBSlave if end node = master node
  if (endslave==false)
  {
    for (int i=0;i<numsnode;++i)
    {
      fac_dxsl_b  += derivsxib(i,0)*(scnodes[i]->xspatial()[0]);
      fac_dysl_b  += derivsxib(i,0)*(scnodes[i]->xspatial()[1]);
      fac_xslm_b += valsxib[i]*(scnodes[i]->xspatial()[0]);
      fac_yslm_b += valsxib[i]*(scnodes[i]->xspatial()[1]);
      fac_dnx_b  += derivsxib(i,0)*(scnodes[i]->n()[0]);
      fac_dny_b  += derivsxib(i,0)*(scnodes[i]->n()[1]);
      fac_nx_b   += valsxib[i]*(scnodes[i]->n()[0]);
      fac_ny_b   += valsxib[i]*(scnodes[i]->n()[1]);
    }

    fac_xslm_b -= mcnodes[0]->xspatial()[0];
    fac_yslm_b -= mcnodes[0]->xspatial()[1];

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
    map<int,double>& nxmap_b = scnodes[0]->GetDerivN()[0];
    map<int,double>& nymap_b = scnodes[0]->GetDerivN()[1];

    // add derivative of slave node coordinates
    dmap_mxib[scnodes[0]->Dofs()[0]] -= scnodes[0]->n()[1];
    dmap_mxib[scnodes[0]->Dofs()[1]] += scnodes[0]->n()[0];

    // add derivatives of master node coordinates
    for (int i=0;i<nummnode;++i)
    {
      dmap_mxib[mcnodes[i]->Dofs()[0]] += valmxib[i]*(scnodes[0]->n()[1]);
      dmap_mxib[mcnodes[i]->Dofs()[1]] -= valmxib[i]*(scnodes[0]->n()[0]);
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
    map<int,double>& nxmap_a = scnodes[1]->GetDerivN()[0];
    map<int,double>& nymap_a = scnodes[1]->GetDerivN()[1];

    // add derivative of slave node coordinates
    dmap_mxia[scnodes[1]->Dofs()[0]] -= scnodes[1]->n()[1];
    dmap_mxia[scnodes[1]->Dofs()[1]] += scnodes[1]->n()[0];

    // add derivatives of master node coordinates
    for (int i=0;i<nummnode;++i)
    {
      dmap_mxia[mcnodes[i]->Dofs()[0]] += valmxia[i]*(scnodes[1]->n()[1]);
      dmap_mxia[mcnodes[i]->Dofs()[1]] -= valmxia[i]*(scnodes[1]->n()[0]);
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
    dmap_sxia[mcnodes[1]->Dofs()[0]] -= fac_ny_a;
    dmap_sxia[mcnodes[1]->Dofs()[1]] += fac_nx_a;

    // add derivatives of slave node coordinates
    for (int i=0;i<numsnode;++i)
    {
      dmap_sxia[scnodes[i]->Dofs()[0]] += valsxia[i]*fac_ny_a;
      dmap_sxia[scnodes[i]->Dofs()[1]] -= valsxia[i]*fac_nx_a;
    }

    // add derivatives of slave node normals
    for (int i=0;i<numsnode;++i)
    {
      map<int,double>& nxmap_curr = scnodes[i]->GetDerivN()[0];
      map<int,double>& nymap_curr = scnodes[i]->GetDerivN()[1];

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
    dmap_sxib[mcnodes[0]->Dofs()[0]] -= fac_ny_b;
    dmap_sxib[mcnodes[0]->Dofs()[1]] += fac_nx_b;

    // add derivatives of slave node coordinates
    for (int i=0;i<numsnode;++i)
    {
      dmap_sxib[scnodes[i]->Dofs()[0]] += valsxib[i]*fac_ny_b;
      dmap_sxib[scnodes[i]->Dofs()[1]] -= valsxib[i]*fac_nx_b;
    }

    // add derivatives of slave node normals
    for (int i=0;i<numsnode;++i)
    {
      map<int,double>& nxmap_curr = scnodes[i]->GetDerivN()[0];
      map<int,double>& nymap_curr = scnodes[i]->GetDerivN()[1];

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
void CONTACT::Integrator::DerivXiGP2D(CONTACT::CElement& sele,
                                    CONTACT::CElement& mele,
                                    double& sxigp, double& mxigp,
                                    const map<int,double>& derivsxi,
                                    map<int,double>& derivmxi)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // we need the participating slave and master nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  vector<CONTACT::CNode*> scnodes(sele.NumNode());
  vector<CONTACT::CNode*> mcnodes(mele.NumNode());
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  for (int i=0;i<numsnode;++i)
  {
    scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
    if (!scnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
  }

  for (int i=0;i<nummnode;++i)
  {
    mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
    if (!mcnodes[i]) dserror("ERROR: DerivXiAB2D: Null pointer!");
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
    sgpn[0]+=valsxigp[i]*scnodes[i]->n()[0];
    sgpn[1]+=valsxigp[i]*scnodes[i]->n()[1];
    sgpn[2]+=valsxigp[i]*scnodes[i]->n()[2];

    sgpx[0]+=valsxigp[i]*scnodes[i]->xspatial()[0];
    sgpx[1]+=valsxigp[i]*scnodes[i]->xspatial()[1];
    sgpx[2]+=valsxigp[i]*scnodes[i]->xspatial()[2];
  }

  // FIXME: This does not have to be the UNIT normal (see 3D)!
  // The reason for this is that we linearize the Gauss point
  // projection from slave to master side here and this condition
  // only includes the Gauss point normal in a cross product.
  // When looking at Projector::ProjectGaussPoint, one can see
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
    fac_dxm_gp += derivmxigp(i,0)*(mcnodes[i]->xspatial()[0]);
    fac_dym_gp += derivmxigp(i,0)*(mcnodes[i]->xspatial()[1]);

    fac_xmsl_gp += valmxigp[i]*(mcnodes[i]->xspatial()[0]);
    fac_ymsl_gp += valmxigp[i]*(mcnodes[i]->xspatial()[1]);
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
    dmap_xsl_gp[scnodes[i]->Dofs()[0]] += valsxigp[i];
    dmap_ysl_gp[scnodes[i]->Dofs()[1]] += valsxigp[i];

    for (CI p=derivsxi.begin();p!=derivsxi.end();++p)
    {
      double facx = derivsxigp(i,0)*(scnodes[i]->xspatial()[0]);
      double facy = derivsxigp(i,0)*(scnodes[i]->xspatial()[1]);
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
    map<int,double>& dmap_nxsl_i = scnodes[i]->GetDerivN()[0];
    map<int,double>& dmap_nysl_i = scnodes[i]->GetDerivN()[1];

    for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
      dmap_nxsl_gp_mod[p->first] += valsxigp[i]*(p->second);
    for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
      dmap_nysl_gp_mod[p->first] += valsxigp[i]*(p->second);

    for (CI p=derivsxi.begin();p!=derivsxi.end();++p)
    {
      double valx =  derivsxigp(i,0)*scnodes[i]->n()[0];
      dmap_nxsl_gp_mod[p->first] += valx*(p->second);
      double valy =  derivsxigp(i,0)*scnodes[i]->n()[1];
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
    derivmxi[mcnodes[i]->Dofs()[0]] += valmxigp[i]*sgpn[1];
    derivmxi[mcnodes[i]->Dofs()[1]] -= valmxigp[i]*sgpn[0];
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
 |  Integrate Mmod on slave / master overlap (2D)             popp 01/08|
 |  This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local      |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an Epetra_SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateMmod2D(CONTACT::CElement& sele,
                                                                 double& sxia, double& sxib,
                                                                 CONTACT::CElement& mele,
                                                                 double& mxia, double& mxib)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateMmod2D called on a wrong type of CElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateMmod2D called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
      dserror("ERROR: IntegrateMmod2D called with infeasible master limits!");

  // create empty mmodseg object and wrap it with RCP
  int nrow  = sele.NumNode();
  int nrowdof = Dim();
  int ncol  = mele.NumNode();
  int ncoldof = Dim();
  RCP<Epetra_SerialDenseMatrix> mmodseg = rcp(new Epetra_SerialDenseMatrix(nrow*nrowdof,ncol*ncoldof));

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);

  // loop over all Gauss points for integration
  for (int gp=0;gp<nGP();++gp)
  {
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // coordinate transformation sxi->eta (slave CElement->Overlap)
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    CONTACT::Projector projector(2);
    projector.ProjectGaussPoint(sele,sxi,mele,mxi);

    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
    {
      cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
      cout << "Gauss point: " << sxi[0] << " " << sxi[1] << endl;
      cout << "Projection: " << mxi[0] << " " << mxi[1] << endl;
      dserror("ERROR: IntegrateMmod2D: Gauss point projection failed!");
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // build the delta function of slave side shape functions
    double deltasval = sval[0]-sval[1];

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    /* loop over all mmodseg matrix entries
       nrow represents the slave Lagrange multipliers !!!
       ncol represents the master dofs !!!
       (this DOES matter here for mmodseg, as it might
       sometimes be rectangular, not quadratic!)              */
    for (int j=0;j<nrow*nrowdof;++j)
    {
      for (int k=0;k<ncol*ncoldof;++k)
      {
        // multiply the two shape functions
        int mindex = (int)(k/ncoldof);
        double prod = 0.5*deltasval*mval[mindex];
        // add current Gauss point's contribution to mmodseg
        (*mmodseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
      }
    }
  } // for (int gp=0;gp<nGP();++gp)

  // prepare computation of purely geometric part of Mmod entries
  CNode* snode0 = static_cast<CNode*>(sele.Nodes()[0]);
  CNode* snode1 = static_cast<CNode*>(sele.Nodes()[1]);

  // normals
  double n[2][2];
  n[0][0] = snode0->n()[0];
  n[0][1] = snode0->n()[1];
  n[1][0] = snode1->n()[0];
  n[1][1] = snode1->n()[1];

  // tangents
  double t[2][2];
  t[0][0] = -n[0][1];
  t[0][1] =  n[0][0];
  t[1][0] = -n[1][1];
  t[1][1] =  n[1][0];

  // scalar product n1 * n2
  double n1n2 = 0.0;
  for (int i=0;i<3;++i)
    n1n2+=n[0][i]*n[1][i];

  // vector product n1 x n2
  double n1xn2 = n[0][0]*n[1][1] - n[0][1]*n[1][0];

  // // multiply geometric part onto Mmod
  for (int i=0;i<ncol;++i)
  {
    (*mmodseg)(0,0+i*ncoldof) *=  (1.0-n1n2);
    (*mmodseg)(1,0+i*ncoldof) *=  n1xn2;
    (*mmodseg)(0,1+i*ncoldof) *= -n1xn2;
    (*mmodseg)(1,1+i*ncoldof) *=  (1.0-n1n2);

    (*mmodseg)(2,0+i*ncoldof) *=  (n1n2-1.0);
    (*mmodseg)(3,0+i*ncoldof) *=  n1xn2;
    (*mmodseg)(2,1+i*ncoldof) *= -n1xn2;
    (*mmodseg)(3,1+i*ncoldof) *=  (n1n2-1.0);
  }

  return mmodseg;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 02/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::IntegrateDerivCell3D(
     CONTACT::CElement& sele, CONTACT::CElement& mele,
     RCP<CONTACT::Intcell> cell,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == Interface::Undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType dt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3D called on a wrong type of CElement pair!");
  if (cell==null)
    dserror("ERROR: IntegrateDerivCell3D called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);

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
  if ((shapefcn_ == Interface::DualFunctions) && (sele.Shape()!=CElement::tri3))
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
    CONTACT::Projector projector(3);
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

    // evaluate dual space shape functions (on slave element)
    if(shapefcn_ == Interface::DualFunctions)
      sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mycnode->n()[0];
      gpn[1]+=sval[i]*mycnode->n()[1];
      gpn[2]+=sval[i]*mycnode->n()[2];

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

    // evaluate the two Jacobians (int. cell and slave element)
    double jaccell = cell->Jacobian(eta);
    double jacslave = sele.Jacobian(sxi);

    // evaluate linearizations *******************************************
    // evaluate the derivative djacdxi = (Jac,xi , Jac,eta)
    double djacdxi[2] = {0.0, 0.0};
    sele.DJacDXi(djacdxi,sxi,ssecderiv);

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
    vector<CONTACT::CNode*> scnodes(sele.NumNode());
    vector<CONTACT::CNode*> mcnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
      if (!scnodes[i]) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
      if (!mcnodes[i]) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;
    map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = scnodes[i]->GetDerivN()[0];
      map<int,double>& dmap_nysl_i = scnodes[i]->GetDerivN()[1];
      map<int,double>& dmap_nzsl_i = scnodes[i]->GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*scnodes[i]->n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*scnodes[i]->n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*scnodes[i]->n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*scnodes[i]->n()[2];
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
        dgapgp[scnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * scnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,1) * scnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[mcnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mcnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,1) * mcnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************

    // reflect the CONTACTONEMORTARLOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef CONTACTONEMORTARLOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef CONTACTONEMORTARLOOP

    // compute cell D/M matrix *******************************************
    
    if (shapefcn_ == Interface::StandardFunctions)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for standard shape functions we use the same algorithm
        // for dseg as in IntegrateDerivSlave2D3D (but with modified integration area)
        // hence, mseg and dseg can not be combined into one loop
        
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = sval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg  
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jaccell*jacslave*wgt;
        }
        
        if (dod)
        {
          for (int k=0; k<nrow*ndof; ++k)
          {
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);
    
            // multiply the two shape functions
            double prod = sval[jindex]*sval[kindex];
    
            // isolate the mseg entries to be filled and
            // add current Gauss point's contribution to mseg  
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                (*dseg)(j, k) += prod*jaccell*jacslave*wgt;
          }
        }
      }
      
    }
    
    else if (shapefcn_ == Interface::DualFunctions)
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
        
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = dualval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg  
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jaccell*jacslave*wgt;
            if (dod) (*dseg)(j, j) += prod*jaccell*jacslave*wgt;
          }
        }
      }
    }
    // compute cell D/M matrix *******************************************

    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = 0.0;
#ifdef CONTACTPETROVGALERKIN
      if (shapefcn_ == Interface::StandardFunctions)
        dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
      prod = sval[j]*gap;
#else
      if (shapefcn_ == Interface::StandardFunctions)
        prod = sval[j]*gap;
      else if (shapefcn_ == Interface::DualFunctions)
        prod = dualval[j]*gap;
#endif // #ifdef CONTACTPETROVGALERKIN

      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jaccell*jacslave*wgt;
    }
    // compute cell gap vector *******************************************

    // compute cell D/M linearization ************************************
    // In the case with separate loops for D and M, this is done by the
    // method IntegrateDerivSlave2D3D(). But in the CONTACTONEMORTARLOOP
    // case we want to combine the linearization of D and M, just as we
    // combine the computation of D and M itself.
    
    // loop over all slave nodes
    for (int j=0; j<nrow; ++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode)
        dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

      int sgid = mycnode->Id();
      
      // for standard shape functions we use the same algorithm
      // for ddmap_jk as in IntegrateDerivSlave2D3D 
      // this means that ddmap_jk and dmmap_jk have to be calculated separately
      if (shapefcn_ == Interface::StandardFunctions)
      {   
        // loop over master nodes for building nodal derivM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;
  
          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];
  
          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions
  
          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*sderiv(j, 0)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sderiv(j, 1)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
  
          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*sval[j]*mderiv(k, 0)*jaccell*jacslave;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sval[j]*mderiv(k, 1)*jaccell*jacslave;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
  
          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*sval[j]*mval[k]*jacslave;
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
          fac = wgt*sval[j]*mval[k]*jaccell;
          for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
  
          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*sval[j]*mval[k]*jaccell*djacdxi[0];
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sval[j]*mval[k]*jaccell*djacdxi[1];
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes
        
        if (dod)
        {
          // loop over slave nodes for building nodal derivD
          for (int k=0; k<nrow; ++k)
          {
            // global master node ID
            int sgid = mele.Nodes()[k]->Id();
            double fac = 0.0;
    
            // get the correct map as a reference
            map<int,double>& ddmap_jk = mycnode->GetDerivD()[sgid];
    
            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions
    
            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*sderiv(j, 0)*sval[k]*jaccell*jacslave;
            for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sderiv(j, 1)*sval[k]*jaccell*jacslave;
            for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
    
            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*sval[j]*sderiv(k, 0)*jaccell*jacslave;
            for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sval[j]*sderiv(k, 1)*jaccell*jacslave;
            for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
    
            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt*sval[j]*sval[k]*jacslave;
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
            fac = wgt*sval[j]*sval[k]*jaccell;
            for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
    
            // (6) Lin(dxdsxi) - slave GP coordinates
            fac = wgt*sval[j]*sval[k]*jaccell*djacdxi[0];
            for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sval[j]*sval[k]*jaccell*djacdxi[1];
            for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
          } // loop over slave nodes
        }
      }
      
      else if (shapefcn_ == Interface::DualFunctions)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jj = mycnode->GetDerivD()[sgid];
  
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;
  
          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];
  
          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0; m<nrow; ++m)
            {
              fac = wgt*sval[m]*mval[k]*jaccell*jacslave;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                if (dod) ddmap_jj[p->first] += fac*(p->second);
              }
            }
  
          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*dualderiv(j, 0)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*dualderiv(j, 1)*mval[k]*jaccell*jacslave;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
  
          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*dualval[j]*mderiv(k, 0)*jaccell*jacslave;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*dualval[j]*mderiv(k, 1)*jaccell*jacslave;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
  
          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*dualval[j]*mval[k]*jacslave;
          for (int m=0; m<(int)jacintcellvec.size(); ++m)
          {
            int v = m/2; // which vertex?
            int dof = m%2; // which dof?
            for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            {
              dmmap_jk[p->first] += fac * jacintcellvec[m] * (p->second);
              if (dod) ddmap_jj[p->first] += fac * jacintcellvec[m] *(p->second);
            }
          }
  
          // (5) Lin(dxdsxi) - slave GP Jacobian
          fac = wgt*dualval[j]*mval[k]*jaccell;
          for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
  
          // (6) Lin(dxdsxi) - slave GP coordinates
          fac = wgt*dualval[j]*mval[k]*jaccell*djacdxi[0];
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*dualval[j]*mval[k]*jaccell*djacdxi[1];
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    } // loop over slave nodes
    // compute cell D/M linearization ************************************

    // compute cell gap linearization ************************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      map<int,double>& dgmap = mycnode->GetDerivG();

#ifdef CONTACTPETROVGALERKIN
      
      if( shapefcn_ == Interface::StandardFunctions )
        dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
      
      // (2) Lin(N) - slave GP coordinates
      fac = wgt*sderiv(j,0)*gap*jaccell*jacslave;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      fac = wgt*sderiv(j,1)*gap*jaccell*jacslave;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);

      // (3) Lin(g) - gap function
      fac = wgt*sval[j]*jaccell*jacslave;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt*sval[j]*gap*jacslave;
      for (int m=0;m<(int)jacintcellvec.size();++m)
      {
        int v = m/2;   // which vertex?
        int dof = m%2; // which dof?
        for (CI p=(cell->GetDerivVertex(v))[dof].begin();p!=(cell->GetDerivVertex(v))[dof].end();++p)
          dgmap[p->first] += fac * jacintcellvec[m] * (p->second);
      }

      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt*sval[j]*gap*jaccell;
      for (CI p=jacslavemap.begin();p!=jacslavemap.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (6) Lin(dxdsxi) - slave GP coordinates
      fac = wgt*sval[j]*gap*jaccell*djacdxi[0];
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      fac = wgt*sval[j]*gap*jaccell*djacdxi[1];
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);
#else
      if( shapefcn_ == Interface::StandardFunctions )
      {
        // (1) Lin(Phi) - dual shape functions
        // this vanishes here since there are no deformation-dependent dual functions
   
         // (2) Lin(NSlave) - slave GP coordinates
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
      else if( shapefcn_ == Interface::DualFunctions )
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
        fac = wgt*dualderiv(j, 0)*gap*jaccell*jacslave;
        for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          dgmap[p->first] += fac*(p->second);
        fac = wgt*dualderiv(j, 1)*gap*jaccell*jacslave;
        for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          dgmap[p->first] += fac*(p->second);
  
        // (3) Lin(g) - gap function
        fac = wgt*dualval[j]*jaccell*jacslave;
        for (CI p=dgapgp.begin(); p!=dgapgp.end(); ++p)
          dgmap[p->first] += fac*(p->second);
  
        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*dualval[j]*gap*jacslave;
        for (int m=0; m<(int)jacintcellvec.size(); ++m)
        {
          int v = m/2; // which vertex?
          int dof = m%2; // which dof?
          for (CI p=(cell->GetDerivVertex(v))[dof].begin(); p!=(cell->GetDerivVertex(v))[dof].end(); ++p)
            dgmap[p->first] += fac * jacintcellvec[m] * (p->second);
        }
  
        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*dualval[j]*gap*jaccell;
        for (CI p=jacslavemap.begin(); p!=jacslavemap.end(); ++p)
          dgmap[p->first] += fac*(p->second);
  
        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*dualval[j]*gap*jaccell*djacdxi[0];
        for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          dgmap[p->first] += fac*(p->second);
        fac = wgt*dualval[j]*gap*jaccell*djacdxi[1];
        for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          dgmap[p->first] += fac*(p->second);
      }
#endif // #ifdef CONTACTPETROVGALERKIN
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
void CONTACT::Integrator::IntegrateDerivCell3DAuxPlane(
     CONTACT::CElement& sele, CONTACT::CElement& mele,
     RCP<CONTACT::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == Interface::Undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of CElement pair!");
  if (cell==null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord(3,sele.NumNode());
  sele.GetNodalCoords(scoord);
  LINALG::SerialDenseMatrix mcoord(3,mele.NumNode());
  mele.GetNodalCoords(mcoord);

  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

  // map iterator
  typedef map<int,double>::const_iterator CI;

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool duallin = false;
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if ( (shapefcn_ == Interface::DualFunctions) && (sele.Shape()!=CElement::tri3))
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
    CONTACT::Projector projector(3);
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

    // evaluate dual space shape functions (on slave element)
    if (shapefcn_ == Interface::DualFunctions)
      sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mycnode->n()[0];
      gpn[1]+=sval[i]*mycnode->n()[1];
      gpn[2]+=sval[i]*mycnode->n()[2];

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

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

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
    vector<CONTACT::CNode*> scnodes(sele.NumNode());
    vector<CONTACT::CNode*> mcnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
      if (!scnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
      if (!mcnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;
    map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = scnodes[i]->GetDerivN()[0];
      map<int,double>& dmap_nysl_i = scnodes[i]->GetDerivN()[1];
      map<int,double>& dmap_nzsl_i = scnodes[i]->GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*scnodes[i]->n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*scnodes[i]->n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*scnodes[i]->n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*scnodes[i]->n()[2];
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
        dgapgp[scnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * scnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,1) * scnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[mcnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mcnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,1) * mcnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************

    // reflect the CONTACTONEMORTARLOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef CONTACTONEMORTARLOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef CONTACTONEMORTARLOOP

    // compute cell D/M matrix *******************************************
    if (shapefcn_ == Interface::StandardFunctions)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for standard shape functions we use the same algorithm
        // for dseg as in IntegrateDerivSlave2D3D (but with modified integration area)
        // hence, mseg and dseg can not be combined into one loop

        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = sval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jac*wgt;
        }

        if (dod)
        {
          for (int k=0; k<nrow*ndof; ++k)
          {
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);

            // multiply the two shape functions
            double prod = sval[jindex]*sval[kindex];

            // isolate the mseg entries to be filled and
            // add current Gauss point's contribution to mseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j, k) += prod*jac*wgt;
          }
        }
      }
    }
    
    else if (shapefcn_ == Interface::DualFunctions)
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

        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = dualval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jac*wgt;
            if (dod) (*dseg)(j, j) += prod*jac*wgt;
          }
        }
      } // nrow*ndof loop
    } // shapefcn_ switch
    // compute cell D/M matrix *******************************************

    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = 0.0;
#ifdef CONTACTPETROVGALERKIN
      if( shapefcn_ == Interface::StandardFunctions )
        dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
      prod = sval[j]*gap;
#else
      if (shapefcn_ == Interface::StandardFunctions)
        prod = sval[j]*gap;
      else if (shapefcn_ == Interface::DualFunctions)
        prod = dualval[j]*gap;
#endif // #ifdef CONTACTPETROVGALERKIN
      
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jac*wgt;
    }
    // compute cell gap vector *******************************************

    // compute cell D/M linearization ************************************
    // In the case with separate loops for D and M, this is done by the
    // method IntegrateDerivSlave2D3D(). But in the CONTACTONEMORTARLOOP
    // case we want to combine the linearization of D and M, just as we
    // combine the computation of D and M itself.
    
    // loop over slave nodes
    for (int j=0; j<nrow; ++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      int sgid = mycnode->Id();

      // for standard shape functions we use the same algorithm
      // for ddmap_jk as in IntegrateDerivSlave2D3D 
      // this means that ddmap_jk and dmmap_jk have to be calculated separately
      if (shapefcn_ == Interface::StandardFunctions)
      {
        // loop over master nodes for building nodal derivM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*sderiv(j, 0)*mval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sderiv(j, 1)*mval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*sval[j]*mderiv(k, 0)*jac;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sval[j]*mderiv(k, 1)*jac;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*sval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes

        if (dod)
        {
          // loop over slave nodes for building nodal derivD
          for (int k=0; k<nrow; ++k)
          {
            // global slave node ID
            int sgid = sele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            map<int,double>& ddmap_jk = mycnode->GetDerivD()[sgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*sderiv(j, 0)*sval[k]*jac;
            for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sderiv(j, 1)*sval[k]*jac;
            for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*sval[j]*sderiv(k, 0)*jac;
            for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sval[j]*sderiv(k, 1)*jac;
            for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt*sval[j]*sval[k];
            for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
          } // loop over slave nodes          
        }
      }
      
      else if (shapefcn_ == Interface::DualFunctions)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jk = mycnode->GetDerivD()[sgid];

        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0; m<nrow; ++m)
            {
              fac = wgt*sval[m]*mval[k]*jac;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                if (dod) ddmap_jk[p->first] += fac*(p->second);
              }
            }

          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*dualderiv(j, 0)*mval[k]*jac;
          for (CI p=dsxigp[0].begin(); p!=dsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jk[p->first] += fac*(p->second);
          }
          fac = wgt*dualderiv(j, 1)*mval[k]*jac;
          for (CI p=dsxigp[1].begin(); p!=dsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jk[p->first] += fac*(p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*dualval[j]*mderiv(k, 0)*jac;
          for (CI p=dmxigp[0].begin(); p!=dmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jk[p->first] += fac*(p->second);
          }
          fac = wgt*dualval[j]*mderiv(k, 1)*jac;
          for (CI p=dmxigp[1].begin(); p!=dmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jk[p->first] += fac*(p->second);
          }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*dualval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jk[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    } // loop over slave nodes
    // compute cell D/M linearization ************************************

    // compute cell gap linearization ************************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      map<int,double>& dgmap = mycnode->GetDerivG();

#ifdef CONTACTPETROVGALERKIN
      if( shapefcn_ == Interface::StandardFunctions )
        dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
      
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
#else
      // (1) Lin(Phi) - dual shape functions
      if (duallin)
        for (int m=0;m<nrow;++m)
        {
          fac = wgt*sval[m]*gap*jac;
          for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
            dgmap[p->first] += fac*(p->second);
        }

      // (2) Lin(Phi) - slave GP coordinates
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualderiv(j,0)*gap*jac;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sderiv(j,0)*gap*jac;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualderiv(j,1)*gap*jac;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sderiv(j,1)*gap*jac;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);

      // (3) Lin(g) - gap function
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*jac;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*jac;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*gap;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*gap;
      for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
        dgmap[p->first] += fac*(p->second);
#endif // #ifdef CONTACTPETROVGALERKIN
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
void CONTACT::Integrator::IntegrateDerivCell3DAuxPlaneQuad(
     CONTACT::CElement& sele, CONTACT::CElement& mele,
     CONTACT::IntElement& sintele, CONTACT::IntElement& mintele,
     RCP<CONTACT::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  
  // explicitely defined shapefunction type needed
  if (shapefcn_ == Interface::Undefined)
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
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of CElement pair!");
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
  LINALG::SerialDenseMatrix ssecderiv(nrow,3);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);
  LINALG::SerialDenseVector sintval(nintrow);
  LINALG::SerialDenseMatrix sintderiv(nintrow,2,true);
  LINALG::SerialDenseVector dualintval(nintrow);
  LINALG::SerialDenseMatrix dualintderiv(nintrow,2,true);

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
  if ( (shapefcn_ == Interface::DualFunctions) && (sele.Shape()!=CElement::tri3))
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }

  // prepare directional derivative of dual shape functions
  // this is necessary for all slave element types except tri3
  bool dualintlin = false;
  vector<vector<map<int,double> > > dualintmap(nintrow,vector<map<int,double> >(nintrow));
  if ( (shapefcn_ == Interface::DualFunctions) && (sintele.Shape()!=CElement::tri3))
  {
    dualintlin = true;
    sintele.DerivShapeDual(dualintmap);
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
    CONTACT::Projector projector(3);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,sintele,sxi,sprojalpha);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,mintele,mxi,mprojalpha);

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

    // evaluate dual space shape functions (on slave element)
    if (shapefcn_ == Interface::DualFunctions)
    {
      sele.EvaluateShapeDual(psxi,dualval,dualderiv,nrow);
      sintele.EvaluateShapeDual(sxi,dualintval,dualintderiv,nintrow);
    }
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(psxi,sval,sderiv,nrow);
    sele.Evaluate2ndDerivShape(psxi,ssecderiv,nrow);
    mele.EvaluateShape(pmxi,mval,mderiv,ncol);
    sintele.EvaluateShape(sxi,sintval,sintderiv,nintrow);

    // build interpolation of slave GP normal and coordinates
    double gpn[3] = {0.0,0.0,0.0};
    double sgpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nrow;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      gpn[0]+=sval[i]*mycnode->n()[0];
      gpn[1]+=sval[i]*mycnode->n()[1];
      gpn[2]+=sval[i]*mycnode->n()[2];

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

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

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
    vector<CONTACT::CNode*> scnodes(sele.NumNode());
    vector<CONTACT::CNode*> mcnodes(mele.NumNode());

    for (int i=0;i<nrow;++i)
    {
      scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
      if (!scnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    for (int i=0;i<ncol;++i)
    {
      mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
      if (!mcnodes[i]) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
    }

    // build directional derivative of slave GP normal (non-unit)
    map<int,double> dmap_nxsl_gp;
    map<int,double> dmap_nysl_gp;
    map<int,double> dmap_nzsl_gp;

    for (int i=0;i<nrow;++i)
    {
      map<int,double>& dmap_nxsl_i = scnodes[i]->GetDerivN()[0];
      map<int,double>& dmap_nysl_i = scnodes[i]->GetDerivN()[1];
      map<int,double>& dmap_nzsl_i = scnodes[i]->GetDerivN()[2];

      for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
        dmap_nxsl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
        dmap_nysl_gp[p->first] += sval[i]*(p->second);
      for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
        dmap_nzsl_gp[p->first] += sval[i]*(p->second);

      for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
      {
        double valx =  sderiv(i,0)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,0)*scnodes[i]->n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,0)*scnodes[i]->n()[2];
        dmap_nzsl_gp[p->first] += valz*(p->second);
      }

      for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
      {
        double valx =  sderiv(i,1)*scnodes[i]->n()[0];
        dmap_nxsl_gp[p->first] += valx*(p->second);
        double valy =  sderiv(i,1)*scnodes[i]->n()[1];
        dmap_nysl_gp[p->first] += valy*(p->second);
        double valz =  sderiv(i,1)*scnodes[i]->n()[2];
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
        dgapgp[scnodes[z]->Dofs()[k]] -= sval[z] * gpn[k];

        for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,0) * scnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
          dgapgp[p->first] -= gpn[k] * sderiv(z,1) * scnodes[z]->xspatial()[k] * (p->second);
      }
    }

    for (int z=0;z<ncol;++z)
    {
      for (int k=0;k<3;++k)
      {
        dgapgp[mcnodes[z]->Dofs()[k]] += mval[z] * gpn[k];

        for (CI p=dpmxigp[0].begin();p!=dpmxigp[0].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,0) * mcnodes[z]->xspatial()[k] * (p->second);

        for (CI p=dpmxigp[1].begin();p!=dpmxigp[1].end();++p)
          dgapgp[p->first] += gpn[k] * mderiv(z,1) * mcnodes[z]->xspatial()[k] * (p->second);
      }
    }
    // evaluate linearizations *******************************************

    // reflect the CONTACTONEMORTARLOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef CONTACTONEMORTARLOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef CONTACTONEMORTARLOOP

    // compute cell D/M matrix *******************************************
    if (shapefcn_ == Interface::StandardFunctions)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for standard shape functions we use the same algorithm
        // for dseg as in IntegrateDerivSlave2D3D (but with modified integration area)
        // hence, mseg and dseg can not be combined into one loop
        
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = sval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jac*wgt;
        }
  
        if (dod)
        {
          for (int k=0; k<nrow*ndof; ++k)
          {
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);
  
            // multiply the two shape functions
            double prod = sval[jindex]*sval[kindex];
  
            // isolate the mseg entries to be filled and
            // add current Gauss point's contribution to mseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j, k) += prod*jac*wgt;
          }
        }
      }
    }
    
    else if (shapefcn_ == Interface::DualFunctions)
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

        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = dualval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jac*wgt;
            if (dod) (*dseg)(j, j) += prod*jac*wgt;
          }
        }
      } // nrow*ndof loop
    } // shapefcn_ switch
    // compute cell D/M matrix *******************************************

    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */

    // WATCH OUT!!!
    // Applying a PetrovGalerkin scheme in the 3D quadratic case, means
    // not only using standard instead of dual shape functions for the
    // weighted gap interpolation, but also reducing polynomial order by one!
#ifdef CONTACTPETROVGALERKIN  
    if( shapefcn_ == Interface::StandardFunctions )
      dserror("Integrator::IntegrateDerivAuxPlaneQuad : CONTACTPETROVGALERKIN flag invalid for standard shape functions");
    
    for (int j=0;j<nintrow;++j)
    {
      double prod = dualintval[j]*gap;
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jac*wgt;
    }
#else
    for (int j=0;j<nrow;++j)
    {
      double prod = 0;
      if (shapefcn_ == Interface::DualFunctions)
        prod = dualval[j]*gap;
      else if (shapefcn_ == Interface::StandardFunctions)
        prod = sval[j]*gap;
      
      // add current Gauss point's contribution to gseg
      (*gseg)(j) += prod*jac*wgt;
    }
#endif // #ifdef CONTACTPETROVGALERKIN

    // compute cell gap vector *******************************************

    // compute cell D/M linearization ************************************
    // In the case with separate loops for D and M, this is done by the
    // method IntegrateDerivSlave2D3D(). But in the CONTACTONEMORTARLOOP
    // case we want to combine the linearization of D and M, just as we
    // combine the computation of D and M itself.
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      int sgid = mycnode->Id();

      // for standard shape functions we use the same algorithm
      // for ddmap_jk as in IntegrateDerivSlave2D3D 
      // this means that ddmap_jk and dmmap_jk have to be calculated separately
      if (shapefcn_ == Interface::StandardFunctions)
      {
        // loop over master nodes for building nodal derivM
        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          // this vanishes here since there are no deformation-dependent dual functions

          // (2) Lin(NSlave) - slave GP coordinates
          fac = wgt*sderiv(j, 0)*mval[k]*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sderiv(j, 1)*mval[k]*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*sval[j]*mderiv(k, 0)*jac;
          for (CI p=dpmxigp[0].begin(); p!=dpmxigp[0].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          fac = wgt*sval[j]*mderiv(k, 1)*jac;
          for (CI p=dpmxigp[1].begin(); p!=dpmxigp[1].end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*sval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
            dmmap_jk[p->first] += fac*(p->second);
        } // loop over master nodes

        if (dod)
        {
          // loop over slave nodes for building nodal derivD
          for (int k=0; k<nrow; ++k)
          {
            // global master node ID
            int sgid = sele.Nodes()[k]->Id();
            double fac = 0.0;

            // get the correct map as a reference
            map<int,double>& ddmap_jk = mycnode->GetDerivD()[sgid];

            // (1) Lin(Phi) - dual shape functions
            // this vanishes here since there are no deformation-dependent dual functions

            // (2) Lin(NSlave) - slave GP coordinates
            fac = wgt*sderiv(j, 0)*sval[k]*jac;
            for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sderiv(j, 1)*sval[k]*jac;
            for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (3) Lin(NSlave) - slave GP coordinates
            fac = wgt*sval[j]*sderiv(k, 0)*jac;
            for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            fac = wgt*sval[j]*sderiv(k, 1)*jac;
            for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);

            // (4) Lin(dsxideta) - intcell GP Jacobian
            fac = wgt*sval[j]*sval[k];
            for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
              ddmap_jk[p->first] += fac*(p->second);
          } // loop over slave nodes        
        }
      }
      
      else if (shapefcn_ == Interface::DualFunctions)
      {
        // get the D-map as a reference
        map<int,double>& ddmap_jj = mycnode->GetDerivD()[sgid];

        for (int k=0; k<ncol; ++k)
        {
          // global master node ID
          int mgid = mele.Nodes()[k]->Id();
          double fac = 0.0;

          // get the correct map as a reference
          map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];

          // (1) Lin(Phi) - dual shape functions
          if (duallin)
            for (int m=0; m<nrow; ++m)
            {
              fac = wgt*sval[m]*mval[k]*jac;
              for (CI p=dualmap[j][m].begin(); p!=dualmap[j][m].end(); ++p)
              {
                dmmap_jk[p->first] += fac*(p->second);
                if (dod) ddmap_jj[p->first] += fac*(p->second);
              }
            }

          // (2) Lin(Phi) - slave GP coordinates
          fac = wgt*dualderiv(j, 0)*mval[k]*jac;
          for (CI p=dpsxigp[0].begin(); p!=dpsxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*dualderiv(j, 1)*mval[k]*jac;
          for (CI p=dpsxigp[1].begin(); p!=dpsxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }

          // (3) Lin(NMaster) - master GP coordinates
          fac = wgt*dualval[j]*mderiv(k, 0)*jac;
          for (CI p=dpmxigp[0].begin(); p!=dpmxigp[0].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
          fac = wgt*dualval[j]*mderiv(k, 1)*jac;
          for (CI p=dpmxigp[1].begin(); p!=dpmxigp[1].end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }

          // (4) Lin(dsxideta) - intcell GP Jacobian
          fac = wgt*dualval[j]*mval[k];
          for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
          {
            dmmap_jk[p->first] += fac*(p->second);
            if (dod) ddmap_jj[p->first] += fac*(p->second);
          }
        } // loop over master nodes
      } // shapefcn_ switch
    } // loop over slave nodes
    // compute cell D/M linearization ************************************

    // compute cell gap linearization ************************************
#ifdef CONTACTPETROVGALERKIN
    
    if (shapefcn_ == Interface::StandardFunctions)
      dserror("CONTACTPETROVGALERKIN flag invalid for standard shape functions");
    
    for (int j=0;j<nintrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(myintnodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      map<int,double>& dgmap = mycnode->GetDerivG();

      // (1) Lin(Phi) - dual shape functions
      if (dualintlin)
        for (int m=0;m<nintrow;++m)
        {
          fac = wgt*sintval[m]*gap*jac;
          for (CI p=dualintmap[j][m].begin();p!=dualintmap[j][m].end();++p)
            dgmap[p->first] += fac*(p->second);
        }

      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt*dualintderiv(j,0)*gap*jac;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      fac = wgt*dualintderiv(j,1)*gap*jac;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);

      // (3) Lin(g) - gap function
      fac = wgt*dualintval[j]*jac;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt*dualintval[j]*gap;
      for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
        dgmap[p->first] += fac*(p->second);
    }
#else
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");

      double fac = 0.0;

      // get the corresponding map as a reference
      map<int,double>& dgmap = mycnode->GetDerivG();

      // (1) Lin(Phi) - dual shape functions
      if (duallin)
        for (int m=0;m<nrow;++m)
        {
          fac = wgt*sval[m]*gap*jac;
          for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
            dgmap[p->first] += fac*(p->second);
        }

      // (2) Lin(Phi) - slave GP coordinates
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualderiv(j,0)*gap*jac;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sderiv(j,0)*gap*jac;
      for (CI p=dpsxigp[0].begin();p!=dpsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualderiv(j,1)*gap*jac;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sderiv(j,1)*gap*jac;
      for (CI p=dpsxigp[1].begin();p!=dpsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);

      // (3) Lin(g) - gap function
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*jac;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*jac;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);

      // (4) Lin(dsxideta) - intcell GP Jacobian
      if (shapefcn_ == Interface::DualFunctions)
        fac = wgt*dualval[j]*gap;
      else if (shapefcn_ == Interface::StandardFunctions)
        fac = wgt*sval[j]*gap;
      for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
        dgmap[p->first] += fac*(p->second);
    }
#endif // #ifdef CONTACTPETROVGALERKIN
    // compute cell gap linearization ************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiGP master (3D)        popp 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivXiGP3D(CONTACT::CElement& sele,
                                      CONTACT::CElement& mele,
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
  vector<CONTACT::CNode*> scnodes(sele.NumNode());
  vector<CONTACT::CNode*> mcnodes(mele.NumNode());
  int numsnode = sele.NumNode();
  int nummnode = mele.NumNode();

  for (int i=0;i<numsnode;++i)
  {
    scnodes[i] = static_cast<CONTACT::CNode*>(snodes[i]);
    if (!scnodes[i]) dserror("ERROR: DerivXiGP3D: Null pointer!");
  }

  for (int i=0;i<nummnode;++i)
  {
    mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
    if (!mcnodes[i]) dserror("ERROR: DerivXiGP3D: Null pointer!");
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
      sgpn[k]+=valsxigp[i]*scnodes[i]->n()[k];
      sgpx[k]+=valsxigp[i]*scnodes[i]->xspatial()[k];
    }

  // build 3x3 factor matrix L
  LINALG::Matrix<3,3> lmatrix(true);
  for (int k=0;k<3;++k) lmatrix(k,2) = -sgpn[k];
  for (int z=0;z<nummnode;++z)
    for (int k=0;k<3;++k)
    {
      lmatrix(k,0) += derivmxigp(z,0) * mcnodes[z]->xspatial()[k];
      lmatrix(k,1) += derivmxigp(z,1) * mcnodes[z]->xspatial()[k];
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
    map<int,double>& dmap_nxsl_i = scnodes[i]->GetDerivN()[0];
    map<int,double>& dmap_nysl_i = scnodes[i]->GetDerivN()[1];
    map<int,double>& dmap_nzsl_i = scnodes[i]->GetDerivN()[2];

    for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
      dmap_nxsl_gp[p->first] += valsxigp[i]*(p->second);
    for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
      dmap_nysl_gp[p->first] += valsxigp[i]*(p->second);
    for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
      dmap_nzsl_gp[p->first] += valsxigp[i]*(p->second);

    for (CI p=derivsxi[0].begin();p!=derivsxi[0].end();++p)
    {
      double valx =  derivsxigp(i,0)*scnodes[i]->n()[0];
      dmap_nxsl_gp[p->first] += valx*(p->second);
      double valy =  derivsxigp(i,0)*scnodes[i]->n()[1];
      dmap_nysl_gp[p->first] += valy*(p->second);
      double valz =  derivsxigp(i,0)*scnodes[i]->n()[2];
      dmap_nzsl_gp[p->first] += valz*(p->second);
    }

    for (CI p=derivsxi[1].begin();p!=derivsxi[1].end();++p)
    {
      double valx =  derivsxigp(i,1)*scnodes[i]->n()[0];
      dmap_nxsl_gp[p->first] += valx*(p->second);
      double valy =  derivsxigp(i,1)*scnodes[i]->n()[1];
      dmap_nysl_gp[p->first] += valy*(p->second);
      double valz =  derivsxigp(i,1)*scnodes[i]->n()[2];
      dmap_nzsl_gp[p->first] += valz*(p->second);
    }
  }

  // start to fill linearization maps for master GP
  // (1) all master nodes coordinates part
  for (int z=0;z<nummnode;++z)
  {
    for (int k=0;k<3;++k)
    {
      derivmxi[0][mcnodes[z]->Dofs()[k]] -= valmxigp[z] * lmatrix(0,k);
      derivmxi[1][mcnodes[z]->Dofs()[k]] -= valmxigp[z] * lmatrix(1,k);
    }
  }

  // (2) slave Gauss point coordinates part
  for (int z=0;z<numsnode;++z)
  {
    for (int k=0;k<3;++k)
    {
      derivmxi[0][scnodes[z]->Dofs()[k]] += valsxigp[z] * lmatrix(0,k);
      derivmxi[1][scnodes[z]->Dofs()[k]] += valsxigp[z] * lmatrix(1,k);

      for (CI p=derivsxi[0].begin();p!=derivsxi[0].end();++p)
      {
        derivmxi[0][p->first] += derivsxigp(z,0) * scnodes[z]->xspatial()[k] * lmatrix(0,k) * (p->second);
        derivmxi[1][p->first] += derivsxigp(z,0) * scnodes[z]->xspatial()[k] * lmatrix(1,k) * (p->second);
      }

      for (CI p=derivsxi[1].begin();p!=derivsxi[1].end();++p)
      {
        derivmxi[0][p->first] += derivsxigp(z,1) * scnodes[z]->xspatial()[k] *lmatrix(0,k) * (p->second);
        derivmxi[1][p->first] += derivsxigp(z,1) * scnodes[z]->xspatial()[k] *lmatrix(1,k) * (p->second);
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
void CONTACT::Integrator::DerivXiGP3DAuxPlane(CONTACT::CElement& ele,
                                        double* xigp, double* auxn,
                                        vector<map<int,double> >& derivxi, double& alpha,
                                        const vector<map<int,double> >& derivauxn,
                                        const vector<map<int,double> >& derivgp)
{
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // we need the participating element nodes
  DRT::Node** nodes = ele.Nodes();
  vector<CONTACT::CNode*> cnodes(ele.NumNode());
  int numnode = ele.NumNode();

  for (int i=0;i<numnode;++i)
  {
    cnodes[i] = static_cast<CONTACT::CNode*>(nodes[i]);
    if (!cnodes[i]) dserror("ERROR: DerivXiGP3DAuxPlane: Null pointer!");
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
      lmatrix(k,0) += derivxigp(z,0) * cnodes[z]->xspatial()[k];
      lmatrix(k,1) += derivxigp(z,1) * cnodes[z]->xspatial()[k];
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
      derivxi[0][cnodes[z]->Dofs()[k]] -= valxigp[z] * lmatrix(0,k);
      derivxi[1][cnodes[z]->Dofs()[k]] -= valxigp[z] * lmatrix(1,k);
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
bool CONTACT::Integrator::AssembleD(const Epetra_Comm& comm,
                                    CONTACT::CElement& sele,
                                    Epetra_SerialDenseMatrix& dseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleD: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
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
        CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(snodes[master]);
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
            snode->AddMValue(sdof,col,minusval);
          }
          else
            snode->AddDValue(sdof,col,val);
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
 |  Assemble M contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave and master  |
 |  overlap pair to the M map of the adjacent slave nodes.              |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleM(const Epetra_Comm& comm,
                                    CONTACT::CElement& sele,
                                    CONTACT::CElement& mele,
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
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
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
        CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sndof+sdof,master*mndof+mdof);
          snode->AddMValue(sdof,col,val);
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
 |  Assemble Mmod contribution (2D)                            popp 01/08|
 |  This method assembles the contribution of a 1D slave / master        |
 |  overlap pair to the Mmod map of the adjacent slave nodes.            |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleMmod(const Epetra_Comm& comm,
                                       CONTACT::CElement& sele,
                                       CONTACT::CElement& mele,
                                       Epetra_SerialDenseMatrix& mmodseg)
{
  // get adjacent slave nodes and master nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleMmod: Null pointer for snodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleMmod: Null pointer for mnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
    //const int* sdofs = snode->Dofs();
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
        CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mmodseg(slave*sndof+sdof,master*mndof+mdof);
          snode->AddMmodValue(sdof,col,val);
        }
      }
    }
    /*
#ifdef DEBUG
    cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
    map<int, double> nodemap0 = (snode->GetMmod())[0];
    map<int, double> nodemap1 = (snode->GetMmod())[1];
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
 |  Assemble g~ contribution (2D / 3D)                        popp 01/08|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleG(const Epetra_Comm& comm,
                                    CONTACT::CElement& sele,
                                    Epetra_SerialDenseVector& gseg)
{

#ifdef CONTACTPETROVGALERKIN
  CONTACT::IntElement& sintref = static_cast<CONTACT::IntElement&>(sele);
#else
  CONTACT::CElement& sintref = sele;
#endif // #ifdef CONTACTPETROVGALERKIN

  // get adjacent slave to assemble to
  DRT::Node** snodes = sintref.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sintref.NumNode();++slave)
  {
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

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

#endif //#ifdef CCADISCRET
