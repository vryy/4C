/*!----------------------------------------------------------------------
\file drt_contact_integrator.cpp
\brief A class to perform integrations of Mortar matrices on the overlap
 			 of two CElements in 1D and 2D

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
CONTACT::Integrator::Integrator(DRT::Element::DiscretizationType eletype)
{
  //*********************************************************************
  // Create integration points according to eletype!
  // Note that our standard Gauss rules are:
  // 5 points: for integrals on 1D lines                 (1,2,3,4,5)
  // 6 points: for integrals on 2D triangles             (1,3,6,7,12,37)
  // 9 points: for integrals on 2D quadrilaterals        (1,4,9)
  //**********************************************************************
  switch(eletype)
  {
  case DRT::Element::line2:
  {
    dim_=2;
    const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_5point);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),1);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::line3:
  {
    dim_=2;
    const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_5point);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),1);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i];
      weights_[i]=intpoints.qwgt[i];
    }
    break;    
  }
  case DRT::Element::tri3:
  {
    dim_=3;
    const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
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
}


/*----------------------------------------------------------------------*
 |  Integrate a 1D slave element overlap                      popp 01/08|
 |  This method integrates 2 functions on the same (slave) CEelement    |
 |  from given local coordinates sxia to sxib                           |
 |  Output is an Epetra_SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateD(CONTACT::CElement& sele,
                                                              double sxia, double sxib)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::IntegrateD not yet implemented for 3D");
  
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: IntegrateD called on a non-slave CElement!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateD called with infeasible slave limits!");
  
  // create empty dseg object and wrap it with RCP
  int nrow = sele.NumNode();
  int ndof = Dim();
  int ncol = nrow;
  
  RCP<Epetra_SerialDenseMatrix> dtemp = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
  
  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,1);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);

  // loop over all Gauss points for integration
  for (int gp=0;gp<nGP();++gp)
  {
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);
    
    // coordinate transformation sxi->eta (CElement->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;
    
    // evaluate trace space and dual space shape functions
    sele.EvaluateShape(sxi,val,deriv,nrow);
    sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);
    
    // evaluate the two Jacobians
    double dxdsxi = sele.Jacobian(eta);
    double dsxideta = -0.5*sxia + 0.5*sxib;
    
    /* loop over all dtemp matrix entries
       nrow represents the Lagrange multipliers !!!
       ncol represents the dofs !!!
       (although this does not really matter here for dseg,
       as it will turn out to be diagonal anyway)             */
    for (int j=0;j<nrow;++j)
    {
      for (int k=0;k<ncol;++k)
      {
        // multiply the two shape functions
        double prod = dualval[j]*val[k];
        // add current Gauss point's contribution to dtemp  
        (*dtemp)(j,k) += prod*dxdsxi*dsxideta*wgt; 
      }
    }  
  } // for (int gp=0;gp<nGP();++gp)

  // fill dseg matrix with dtemp matrix entries
  // (each dtemp value is multiplied with a (dof)-unit-matrix)
  for (int j=0;j<nrow*ndof;++j)
  {
    for (int k=0;k<ncol*ndof;++k)
    {
      int jindex = (int)(j/ndof);
      int kindex = (int)(k/ndof);
      // isolate the dseg entries to be filled
      if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
        (*dseg)(j,k) = (*dtemp)(jindex,kindex);
    }
  }
  
  return dseg;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of D                       popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivD(CONTACT::CElement& sele,
                                 double sxia, double sxib)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::DerivD not yet implemented for 3D");
    
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: DerivD called on a non-slave CElement!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: DerivD called with infeasible slave limits!");
  
  int nrow = sele.NumNode();
  
  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,1);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);
  
  // get nodal coords for Jacobian evaluation
  LINALG::SerialDenseMatrix coord = sele.GetNodalCoords();
  
  // prepare directional derivative of dual shape functions
  // this is only necessary for qudratic shape functions in 2D
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if (sele.Shape()==CElement::line3)
    sele.DerivShapeDual(dualmap);
  
  // loop over all Gauss points for integration
  for (int gp=0;gp<nGP();++gp)
  {
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);
    typedef map<int,double>::const_iterator CI;
    
    // coordinate transformation sxi->eta (CElement->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;
    
    // evaluate trace space and dual space shape functions
    sele.EvaluateShape(sxi,val,deriv,nrow);
    sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);
    
    // evaluate the two Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the Jacobian derivative
    map<int,double> testmap;
    sele.DerivJacobian(val,deriv,coord,testmap);
    
    // compute contribution oj J to nodal D-derivative-maps
    DRT::Node** mynodes = sele.Nodes();
    if (!mynodes) dserror("ERROR: IntegrateD: Null pointer!");
    
    for (int i=0;i<nrow;++i)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[i]);
      if (!mycnode) dserror("ERROR: Integrate1D: Null pointer!");
      bool bound = mycnode->IsOnBound();
      map<int,double>& nodemap = mycnode->GetDerivD();
      
      //******************************************************************
      // standard case (no edge node modification)
      //******************************************************************
      if (!bound)
      {
        // contribution of current element / current GP
        double fac = wgt*val[i]*dualval[i]*dsxideta;
        for (CI p=testmap.begin();p!=testmap.end();++p)
          nodemap[p->first] += fac*(p->second);
      }
      
      //******************************************************************
      // edge node modification case
      //******************************************************************
      else
      {
        // get gid of current boundary node
        int bgid = mycnode->Id();
        
        // loop over other nodes (interior nodes)
        for (int k=0;k<nrow;++k)
        {
          CONTACT::CNode* mycnode2 = static_cast<CONTACT::CNode*>(mynodes[k]);
          if (!mycnode2) dserror("ERROR: Integrate1D: Null pointer!");
          bool bound2 = mycnode2->IsOnBound();
          if (bound2) continue;
          map<int,double>& nodemmap = mycnode2->GetDerivM()[bgid];
          
          // contribution to DerivM of current element / current GP
          double fac = wgt*val[i]*dualval[k]*dsxideta;
          for (CI p=testmap.begin();p!=testmap.end();++p)
            nodemmap[p->first] -= fac*(p->second);
        }
      }
    }
    
    // compute contribution of dual shape fct. to nodal D-derivative-maps
    if (sele.Shape()!=CElement::line3) continue;
    for (int i=0;i<nrow;++i)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[i]);
      if (!mycnode) dserror("ERROR: Integrate1D: Null pointer!");
      bool bound = mycnode->IsOnBound();
      map<int,double>& nodemap = mycnode->GetDerivD();
      
      //******************************************************************
      // standard case (no edge node modification)
      //******************************************************************
      if (!bound)
      {
        // contribution of current element / current GP
        for (int j=0;j<nrow;++j)
        {
          double fac = wgt*val[i]*val[j]*dsxideta*dxdsxi;
          for (CI p=dualmap[i][j].begin();p!=dualmap[i][j].end();++p)
            nodemap[p->first] += fac*(p->second);
        }
      }
      
      //******************************************************************
      // edge node modification case
      //******************************************************************
      else
      {
        dserror("ERROR: edge node modification + quad elements + full lin not yet working");
        // get gid of current boundary node
        int bgid = mycnode->Id();
        
        // loop over other nodes (interior nodes)
        for (int k=0;k<nrow;++k)
        {
          CONTACT::CNode* mycnode2 = static_cast<CONTACT::CNode*>(mynodes[k]);
          if (!mycnode2) dserror("ERROR: Integrate1D: Null pointer!");
          bool bound2 = mycnode2->IsOnBound();
          if (bound2) continue;
          map<int,double>& nodemmap = mycnode2->GetDerivM()[bgid];
          
          // contribution of current element / current GP
          for (int j=0;j<nrow;++j)
          {
            LINALG::SerialDenseVector vallin(nrow-1);
            LINALG::SerialDenseMatrix derivlin(nrow-1,1);
            if (i==0) sele.ShapeFunctions(CElement::dual1D_base_for_edge0,sxi,vallin,derivlin);
            else if (i==1) sele.ShapeFunctions(CElement::dual1D_base_for_edge1,sxi,vallin,derivlin);
            double fac = wgt*val[i]*vallin[j]*dsxideta*dxdsxi;
            for (CI p=dualmap[k][j].begin();p!=dualmap[k][j].end();++p)
              nodemmap[p->first] -= fac*(p->second);
          }
        }
      }
    }
    
  } // for (int gp=0;gp<nGP();++gp)
  
  return;
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave / master overlap                     popp 01/08|
 |  This method integrates a slave side function (dual shape fct.)      |
 |  and a master side function (standard shape fct.) from given local   |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an Epetra_SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateM(CONTACT::CElement& sele,
                                                              double sxia, double sxib,
                                                              CONTACT::CElement& mele,
                                                              double mxia, double mxib)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::IntegrateM not yet implemented for 3D");
    
  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateM called on a wrong type of CElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateM called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
      dserror("ERROR: IntegrateM called with infeasible master limits!");
  
  // create empty mseg object and wrap it with RCP
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();
  
  RCP<Epetra_SerialDenseMatrix> mtemp = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
  RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
  
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
    
    // simple version (no Gauss point projection)
    // double test[2] = {0.0, 0.0};
    // test[0] = 0.5*(1-eta[0])*mxib + 0.5*(1+eta[0])*mxia;
    // mxi[0]=test[0];
    
    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
    {
      cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
      cout << "Slave nodes: " << sele.NodeIds()[0] << " " << sele.NodeIds()[1] << endl;
      cout << "Master nodes: " << mele.NodeIds()[0] << " " << mele.NodeIds()[1] << endl;
      cout << "sxia: " << sxia << " sxib: " << sxib << endl;
      cout << "mxia: " << mxia << " mxib: " << mxib << endl;
      dserror("ERROR: IntegrateM: Gauss point projection failed! mxi=%d",mxi[0]);
    }
    // evaluate dual space shape functions (on slave element)
    sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);
    
    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;
    
    /* loop over all mseg matrix entries
       nrow represents the slave Lagrange multipliers !!!
       ncol represents the master dofs !!!
       (this DOES matter here for mseg, as it might
       sometimes be rectangular, not quadratic!)              */
    for (int j=0;j<nrow;++j)
    {
      for (int k=0;k<ncol;++k)
      {
        // multiply the two shape functions
        double prod = dualval[j]*mval[k];
        // add current Gauss point's contribution to mseg  
        (*mtemp)(j,k) += prod*dxdsxi*dsxideta*wgt; 
      }
    }  
  } // for (int gp=0;gp<nGP();++gp)
  
  // fill mseg matrix with mtemp matrix entries
  // (each mtemp value is multiplied with a (dof)-unit-matrix)
  for (int j=0;j<nrow*ndof;++j)
  {
    for (int k=0;k<ncol*ndof;++k)
    {
      int jindex = (int)(j/ndof);
      int kindex = (int)(k/ndof);
      // isolate the mseg entries to be filled
      if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
        (*mseg)(j,k) = (*mtemp)(jindex,kindex);
    }
  }

  return mseg;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of M                       popp 05/08|
 |  IMPORTANT NOTE:                                                     |
 |  If CONTACTONEMORTARLOOP is defined then this method also computes   |
 |  the contribution of a 1D slave element part to the D-derivative     |
 |  map of the adjacent slave nodes                                     |
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivM(CONTACT::CElement& sele,
                                 double sxia, double sxib,
                                 CONTACT::CElement& mele,
                                 double mxia, double mxib)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::DerivM not yet implemented for 3D");
    
  // *********************************************************************
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
  typedef map<int,double>::const_iterator CI;
  
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
  DerivXiAB(sele,sxia,sxib,mele,mxia,mxib,ximaps,startslave,endslave);
  
  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: DerivM called on a wrong type of CElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: DerivM called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: DerivM called with infeasible master limits!");
  
  // create empty mseg object and wrap it with RCP
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  
  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseMatrix ssecderiv(nrow,1); //only 2D so far
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);

  // get slave nodal coords for Jacobian evaluation
  LINALG::SerialDenseMatrix scoord = sele.GetNodalCoords();
  
  // prepare directional derivative of dual shape functions
  // this is only necessary for qudratic shape functions in 2D
  vector<vector<map<int,double> > > dualmap(nrow,vector<map<int,double> >(nrow));
  if (sele.Shape()==CElement::line3)
    sele.DerivShapeDual(dualmap);
    
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
    
    // simple version (no Gauss point projection)
    // double test[2] = {0.0, 0.0};
    // test[0] = 0.5*(1-eta[0])*mxib + 0.5*(1+eta[0])*mxia;
    // mxi[0]=test[0];
    
    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
    {
      cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
      cout << "Slave nodes: " << sele.NodeIds()[0] << " " << sele.NodeIds()[1] << endl;
      cout << "Master nodes: " << mele.NodeIds()[0] << " " << mele.NodeIds()[1] << endl;
      cout << "sxia: " << sxia << " sxib: " << sxib << endl;
      cout << "mxia: " << mxia << " mxib: " << mxib << endl;
      dserror("ERROR: IntegrateM: Gauss point projection failed! mxi=%d",mxi[0]);
    }
    // evaluate dual space shape functions (on slave element)
    sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    sele.Evaluate2ndDerivShape(sxi,ssecderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);
    
    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the derivative dxdsxidsxi = Jac,xi
    double djacdxi[2] = {0.0, 0.0};
    sele.DJacDXi(djacdxi,sval,sderiv,ssecderiv,scoord);
    double dxdsxidsxi=djacdxi[0]; // only 2D so far
    
    // evalute the GP slave coordinate derivatives
    map<int,double> dsxigp;
    for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
      dsxigp[p->first] += 0.5*(1-eta[0])*(p->second);
    for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
      dsxigp[p->first] += 0.5*(1+eta[0])*(p->second);
    
    // evalute the GP master coordinate derivatives
    map<int,double> dmxigp;
    DerivXiGP(sele,sxia,sxib,mele,mxia,mxib,sxi[0],mxi[0],dsxigp,dmxigp);
    
    // simple version (no Gauss point projection)
    //for (CI p=ximaps[2].begin();p!=ximaps[2].end();++p)
    //  dmxigp[p->first] += 0.5*(1+eta[0])*(p->second);
    //for (CI p=ximaps[3].begin();p!=ximaps[3].end();++p)
    //  dmxigp[p->first] += 0.5*(1-eta[0])*(p->second);
   
    // evaluate the Jacobian derivative
    map<int,double> testmap;
    sele.DerivJacobian(sval,sderiv,scoord,testmap);
    
    // contributions to DerivM_jk
    DRT::Node** mynodes = sele.Nodes();
    if (!mynodes) dserror("ERROR: DerivM: Null pointer!");
        
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: DerivM: Null pointer!");
            
      for (int k=0;k<ncol;++k)
      {
        // global master node ID
        int mgid = mele.Nodes()[k]->Id();
        double fac = 0.0;
        
        // get the correct map as a reference
        map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];
        
        // (1) Lin(Phi) - dual shape functions
        for (int m=0;m<nrow;++m)
        {
          fac = wgt*sval[m]*mval[k]*dsxideta*dxdsxi;
          for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
            dmmap_jk[p->first] += fac*(p->second);
        }
        
        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*dualderiv(j,0)*mval[k]*dsxideta*dxdsxi;
        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt*dualval(j,0)*mderiv(k,0)*dsxideta*dxdsxi;
        for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (4) Lin(dsxideta) - segment end coordinates
        fac = wgt*dualval[j]*mval[k]*dxdsxi;
        for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
          dmmap_jk[p->first] -= 0.5*fac*(p->second);
        for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
          dmmap_jk[p->first] += 0.5*fac*(p->second);
        
        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*dualval[j]*mval[k]*dsxideta;
        for (CI p=testmap.begin();p!=testmap.end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*dualval[j]*mval[k]*dsxideta*dxdsxidsxi;
        for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
          dmmap_jk[p->first] += fac*(p->second);
      }
    }
    
    //********************************************************************
    // Compute directional derivative of D
    //********************************************************************
    // In the case with separate loops for D and M, this is done by the
    // method DerivD() above. But in the CONTACTONEMORTARLOOP case we
    // want to combine the linearization of D and M, just as we combine
    // the computation of D and M itself.
    // Of course, this measn we get additional terms for the linearization
    // of D because the slave segment end coordinates can now vary. In
    // the case with separate D and M loops we simply could integrate on
    // the full slave element from xi=-1 to xi=1, now we have to consider
    // the coordinates sxia,sxib and their linearizations!
    //********************************************************************
#ifdef CONTACTONEMORTARLOOP
    // contributions to DerivD_jj
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: DerivM: Null pointer!");
      double fac = 0.0;
      
      // get the D-map as a reference
      map<int,double>& ddmap_jj = mycnode->GetDerivD();
      
      // (1) Lin(Phi) - dual shape functions
      for (int k=0;k<nrow;++k)
      {
        fac = wgt*sval[j]*sval[k]*dsxideta*dxdsxi;
        for (CI p=dualmap[j][k].begin();p!=dualmap[j][k].end();++p)
          ddmap_jj[p->first] += fac*(p->second);
      }
      
      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt*dualderiv(j,0)*sval[j]*dsxideta*dxdsxi;
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        ddmap_jj[p->first] += fac*(p->second);
      
      // (3) Lin(NSlave) - slave GP coordinates
      fac = wgt*dualval[j]*sderiv(j,0)*dsxideta*dxdsxi;
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        ddmap_jj[p->first] += fac*(p->second);
      
      // (4) Lin(dsxideta) - segment end coordinates
      fac = wgt*dualval[j]*sval[j]*dxdsxi;
      for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
        ddmap_jj[p->first] -= 0.5*fac*(p->second);
      for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
        ddmap_jj[p->first] += 0.5*fac*(p->second);
      
      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt*dualval[j]*sval[j]*dsxideta;
      for (CI p=testmap.begin();p!=testmap.end();++p)
        ddmap_jj[p->first] += fac*(p->second);
      
      // (6) Lin(dxdsxi) - slave GP coordinates
      fac = wgt*dualval[j]*sval[j]*dsxideta*dxdsxidsxi;
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        ddmap_jj[p->first] += fac*(p->second);
    }
#endif // #ifdef CONTACTONEMORTARLOOP
    
  } // for (int gp=0;gp<nGP();++gp)
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of XiAB                    popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivXiAB(CONTACT::CElement& sele,
                                    double sxia, double sxib,
                                    CONTACT::CElement& mele,
                                    double mxia, double mxib,
                                    vector<map<int,double> >& derivxi,
                                    bool startslave, bool endslave)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::DerivXiAB not yet implemented for 3D");
    
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
    if (!scnodes[i]) dserror("ERROR: DerivXiAB: Null pointer!");
  }
  
  for (int i=0;i<nummnode;++i)
  {
    mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
    if (!mcnodes[i]) dserror("ERROR: DerivXiAB: Null pointer!");
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
 |  Compute directional derivative of XiGP master             popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::DerivXiGP(CONTACT::CElement& sele,
                                    double sxia, double sxib,
                                    CONTACT::CElement& mele,
                                    double mxia, double mxib,
                                    double sxigp, double mxigp,
                                    const map<int,double>& derivsxi,
                                    map<int,double>& derivmxi)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::DerivXiGP not yet implemented for 3D");
    
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
    if (!scnodes[i]) dserror("ERROR: DerivXiAB: Null pointer!");
  }
  
  for (int i=0;i<nummnode;++i)
  {
    mcnodes[i] = static_cast<CONTACT::CNode*>(mnodes[i]);
    if (!mcnodes[i]) dserror("ERROR: DerivXiAB: Null pointer!");
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
  
  // normalize interpolated GP normal back to length 1.0 !!!
  double length = sqrt(sgpn[0]*sgpn[0]+sgpn[1]*sgpn[1]+sgpn[2]*sgpn[2]);
  if (length<1.0e-12) dserror("ERROR: DerivXiGP: Divide by zero!");
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
 |  Integrate a 1D slave / master overlap                     popp 01/08|
 |  This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local      |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an Epetra_SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateMmod(CONTACT::CElement& sele,
                                                                 double sxia, double sxib,
                                                                 CONTACT::CElement& mele,
                                                                 double mxia, double mxib)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::IntegrateMmod not yet implemented for 3D");
    
  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateMmod called on a wrong type of CElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateMmod called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
      dserror("ERROR: IntegrateMmod called with infeasible master limits!");
  
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
      dserror("ERROR: IntegrateMmod: Gauss point projection failed!");
    
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
 |  Integrate gap on a 1D slave / master overlap              popp 01/08|
 |  This method integrates a slave side function (dual shape fct.)      |
 |  and the gap function g = ( ( sx - mx ) * n )  from given local      |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an Epetra_SerialDenseVector holding the int. values       |
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseVector> CONTACT::Integrator::IntegrateG(CONTACT::CElement& sele,
                                                              double sxia, double sxib,
                                                              CONTACT::CElement& mele,
                                                              double mxia, double mxib)
{
  //check for problem dimension
  if (Dim()==3) dserror("ERROR: Integrator::IntegrateG not yet implemented for 3D");
    
  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateG called on a wrong type of CElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateG called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
      dserror("ERROR: IntegrateG called with infeasible master limits!");
  
  // create empty gseg object and wrap it with RCP
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
  
  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);

  // get slave and master nodal coords for Jacobian / GP evaluation
  LINALG::SerialDenseMatrix scoord = sele.GetNodalCoords();
  LINALG::SerialDenseMatrix mcoord = mele.GetNodalCoords();
  
  // get slave element nodes themselves for normal evaluation
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateG: Null pointer!");
  
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
    
    // simple version (no Gauss point projection)
    // double test[2] = {0.0, 0.0};
    // test[0] = 0.5*(1-eta[0])*mxib + 0.5*(1+eta[0])*mxia;
    // mxi[0]=test[0];
    
    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
      dserror("ERROR: IntegrateG: Gauss point projection failed!");
    
    // evaluate dual space shape functions (on slave element)
    sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
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
    if (length<1.0e-12) dserror("ERROR: IntegrateG: Divide by zero!");
    
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
    
    // build normal gap at current GP
    double gap = 0.0;
    for (int i=0;i<3;++i)
      gap+=(mgpx[i]-sgpx[i])*gpn[i];
    
#ifdef DEBUG
    //cout << "GP gap: " << gap << endl;
#endif // #ifdef DEBUG
    
    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;
    
    /* loop over all gseg vector entries
       nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = dualval[j]*gap;
      // add current Gauss point's contribution to gseg  
      (*gseg)(j) += prod*dxdsxi*dsxideta*wgt; 
    }
    
  } // for (int gp=0;gp<nGP();++gp)

  return gseg;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution                                   popp 01/08|
 |  This method assembles the contrubution of a 1D slave element        |
 |  to the D map of the adjacent slave nodes.                           |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleD(CONTACT::Interface& inter,
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
    if (snode->Owner() != inter.Comm().MyPID())
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
            snode->AddMValue(sdof,col,-val);
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
 |  Assemble M contribution                                   popp 01/08|
 |  This method assembles the contrubution of a 1D slave / master       |
 |  overlap pair to the M map of the adjacent slave nodes.              |
 |  IMPORTANT NOTE:                                                     |
 |  If CONTACTONEMORTARLOOP is defined then this method also assembles  |
 |  the contribution of a 1D slave element part to the D map of the     |
 |  adjacent slave nodes via the connection D = sum (M)                 |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleM(CONTACT::Interface& inter,
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
#ifdef CONTACTONEMORTARLOOP
    const int* sdofs = snode->Dofs();
#endif // #ifdef CONTACTONEMORTARLOOP
    
    int sndof = snode->NumDof();
    
    // only process slave node rows that belong to this proc
    if (snode->Owner() != inter.Comm().MyPID())
      continue;
    
    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;
    
    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // initialize sum(M) of current row
      double msum=0.0;
      
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
          msum+=val;
          snode->AddMValue(sdof,col,val);
        }
      }
      
#ifdef CONTACTONEMORTARLOOP
      // we know that D(snode) = sum_mnode (M(mnode,snode))
      snode->AddDValue(sdof,sdofs[sdof],msum);
  #ifdef CONTACTBOUNDMOD
      dserror("ERROR: Combination 1 mortar loop <-> boundary modification not yet impl.");
  #endif // #ifdef CONTACTBOUNDMOD
#endif // #ifdef CONTACTONEMORTARLOOP
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
 |  Assemble Mmod contribution                                popp 01/08|
 |  This method assembles the contribution of a 1D slave / master        |
 |  overlap pair to the Mmod map of the adjacent slave nodes.            |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleMmod(CONTACT::Interface& inter,
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
    if (snode->Owner() != inter.Comm().MyPID())
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
 |  Assemble g~ contribution                                  popp 01/08|
 |  This method assembles the contribution of a 1D slave / master        |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleG(CONTACT::Interface& inter,
                                    CONTACT::CElement& sele,
                                    Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleG: Null pointer for snodes!");
  
  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
    
    // only process slave node rows that belong to this proc
    if (snode->Owner() != inter.Comm().MyPID())
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
