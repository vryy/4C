/*!----------------------------------------------------------------------
\file drt_contact_integrator_tools.cpp
\brief

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
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 02/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::IntegrateDerivCell3D(
     CONTACT::CElement& sele, CONTACT::CElement& mele,
     RCP<CONTACT::Intcell> cell,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg,
     vector<vector<double> >& testgps, vector<vector<double> >& testgpm,
     vector<vector<double> >& testjs, vector<vector<double> >& testji,
     bool printderiv)
{
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
  if (sele.Shape()!=CElement::tri3)
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }
  
  // finite difference checks
  vector<double> testgpscurr(12);
  vector<double> testgpmcurr(12);
  vector<double> testjscurr(6);
  vector<double> testjicurr(6);
  static int gpcounter = 0;
    
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
    
    // store results in testgps and testgpm for FD checks
    // store results in testjs and testji for FD checks
    testgpscurr[2*gp] = sxi[0];
    testgpscurr[2*gp+1] = sxi[1];
    testgpmcurr[2*gp] = mxi[0];
    testgpmcurr[2*gp+1] = mxi[1];
    testjscurr[gp] = jacslave;
    testjicurr[gp] = jaccell;
        
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
    LINALG::SerialDenseVector svalcell(nrow);
    LINALG::SerialDenseMatrix sderivcell(nrow,2,true);
    cell->EvaluateShape(eta,svalcell,sderivcell);
    
    for (int v=0;v<cell->NumVertices();++v)
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
    
    // finite difference checks (print to screen)
    if (printderiv)
    {
      // prepare Jacobian derivatives
      map<int,double> djsgp;
      map<int,double> djigp;
      
      for (CI p=jacslavemap.begin();p!=jacslavemap.end();++p)
        djsgp[p->first] += (p->second);
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        djsgp[p->first] += djacdxi[0] * (p->second);
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        djsgp[p->first] += djacdxi[1] * (p->second);
      
      for (int m=0;m<(int)jacintcellvec.size();++m)
      {
        int v = m/2;   // which vertex?
        int dof = m%2; // which dof?
        for (CI p=(cell->GetDerivVertex(v))[dof].begin();p!=(cell->GetDerivVertex(v))[dof].end();++p)
          djigp[p->first] += jacintcellvec[m] * (p->second);
      }
      
      typedef map<int,double>::const_iterator CI;

      cout << "\nAnalytical derivative for Intcell " << gpcounter << " (x-component). Slave GP " << gp << endl;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      cout << "Analytical derivative for Intcell " << gpcounter << " (y-component). Slave GP " << gp << endl;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;

      cout << "\nAnalytical derivative for Intcell " << gpcounter << " (x-component). Master GP " << gp << endl;
      for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      cout << "Analytical derivative for Intcell " << gpcounter << " (y-component). Master GP " << gp << endl;
      for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      
      cout << "\nAnalytical derivative for Slave Jacobian Intcell " << gpcounter << " GP " << gp << endl;
      for (CI p=djsgp.begin();p!=djsgp.end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      
      cout << "\nAnalytical derivative for Intcell Jacobian Intcell " << gpcounter << " GP " << gp << endl;
      for (CI p=djigp.begin();p!=djigp.end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
    }
        
    // compute cell M matrix *********************************************
    // loop over all mseg matrix entries
    // nrow represents the slave Lagrange multipliers !!!
    // ncol represents the master dofs !!!
    // (this DOES matter here for mseg, as it might
    // sometimes be rectangular, not quadratic!)
    for (int j=0;j<nrow*ndof;++j)
    {
      for (int k=0;k<ncol*ndof;++k)
      {
        int jindex = (int)(j/ndof);
        int kindex = (int)(k/ndof);
              
        // multiply the two shape functions
        double prod = dualval[jindex]*mval[kindex];
        
        // isolate the mseg entries to be filled and
        // add current Gauss point's contribution to mseg  
        if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          (*mseg)(j,k) += prod*jaccell*jacslave*wgt; 
      }
    }  
    // compute cell M matrix *********************************************
    
    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = dualval[j]*gap;
      // add current Gauss point's contribution to gseg  
      (*gseg)(j) += prod*jaccell*jacslave*wgt; 
    }
    // compute cell gap vector *******************************************
    
    // compute cell M linearization **************************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
            
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
          fac = wgt*sval[m]*mval[k]*jaccell*jacslave;
          for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
            dmmap_jk[p->first] += fac*(p->second);
        }
        
        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*dualderiv(j,0)*mval[k]*jaccell*jacslave;
        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        fac = wgt*dualderiv(j,1)*mval[k]*jaccell*jacslave;
        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt*dualval[j]*mderiv(k,0)*jaccell*jacslave;
        for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        fac = wgt*dualval[j]*mderiv(k,1)*jaccell*jacslave;
        for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*dualval[j]*mval[k]*jacslave;
        for (int m=0;m<(int)jacintcellvec.size();++m)
        {
          int v = m/2;   // which vertex?
          int dof = m%2; // which dof?
          for (CI p=(cell->GetDerivVertex(v))[dof].begin();p!=(cell->GetDerivVertex(v))[dof].end();++p)
            dmmap_jk[p->first] += fac * jacintcellvec[m] * (p->second);
        }
        
        // (5) Lin(dxdsxi) - slave GP Jacobian
        fac = wgt*dualval[j]*mval[k]*jaccell;
        for (CI p=jacslavemap.begin();p!=jacslavemap.end();++p)
          dmmap_jk[p->first] += fac*(p->second);
                
        // (6) Lin(dxdsxi) - slave GP coordinates
        fac = wgt*dualval[j]*mval[k]*jaccell*djacdxi[0];
        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        fac = wgt*dualval[j]*mval[k]*jaccell*djacdxi[1];
        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
      }
    }
    // compute cell M linearization **************************************
    
    // compute cell gap linearization ************************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3D: Null pointer!");
            
      double fac = 0.0;
      
      // get the corresponding map as a reference
      map<int,double>& dgmap = mycnode->GetDerivG();
      
      // (1) Lin(Phi) - dual shape functions
      for (int m=0;m<nrow;++m)
      {
        fac = wgt*sval[m]*gap*jaccell*jacslave;
        for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
          dgmap[p->first] += fac*(p->second);
      }
      
      // (2) Lin(Phi) - slave GP coordinates
      fac = wgt*dualderiv(j,0)*gap*jaccell*jacslave;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      fac = wgt*dualderiv(j,1)*gap*jaccell*jacslave;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);
      
      // (3) Lin(g) - gap function
      fac = wgt*dualval[j]*jaccell*jacslave;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);
      
      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt*dualval[j]*gap*jacslave;
      for (int m=0;m<(int)jacintcellvec.size();++m)
      {
        int v = m/2;   // which vertex?
        int dof = m%2; // which dof?
        for (CI p=(cell->GetDerivVertex(v))[dof].begin();p!=(cell->GetDerivVertex(v))[dof].end();++p)
          dgmap[p->first] += fac * jacintcellvec[m] * (p->second);
      }
      
      // (5) Lin(dxdsxi) - slave GP Jacobian
      fac = wgt*dualval[j]*gap*jaccell;
      for (CI p=jacslavemap.begin();p!=jacslavemap.end();++p)
        dgmap[p->first] += fac*(p->second);
              
      // (6) Lin(dxdsxi) - slave GP coordinates
      fac = wgt*dualval[j]*gap*jaccell*djacdxi[0];
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      fac = wgt*dualval[j]*gap*jaccell*djacdxi[1];
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);
    }
    // compute cell gap linearization ************************************ 
  }
  //**********************************************************************
  
  // finite difference checks
  testgps.push_back(testgpscurr);
  testgpm.push_back(testgpmcurr);
  testjs.push_back(testjscurr);
  testji.push_back(testjicurr);
  if (printderiv) gpcounter = gpcounter +1;
  
#ifdef CONTACTONEMORTARLOOP
    dserror("ERROR: IntegrateDerivCell3D: One mortar loop case not yet impl. for 3D!");
#endif // #ifdef CONTACTONEMORTARLOOP
    
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
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
void CONTACT::Integrator::IntegrateDerivCell3DAuxPlane(
     CONTACT::CElement& sele, CONTACT::CElement& mele,
     RCP<CONTACT::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg,
     vector<vector<double> >& testgps, vector<vector<double> >& testgpm,
     vector<vector<double> >& testjs, vector<vector<double> >& testji,
     bool printderiv)
{
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
  if (sele.Shape()!=CElement::tri3)
  {
    duallin = true;
    sele.DerivShapeDual(dualmap);
  }
  
  // finite difference checks
  vector<double> testgpscurr(12);
  vector<double> testgpmcurr(12);
  vector<double> testjscurr(6);
  vector<double> testjicurr(6);
  static int gpcounter = 0;
    
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

    // store results in testgps and testgpm for FD checks
    // store results in testjs and testji for FD checks
    testgpscurr[2*gp] = sxi[0];
    testgpscurr[2*gp+1] = sxi[1];
    testgpmcurr[2*gp] = mxi[0];
    testgpmcurr[2*gp+1] = mxi[1];
    testjscurr[gp] = 0.0;
    testjicurr[gp] = jac;
        
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
    
    // finite difference checks (print to screen)
    if (printderiv)
    {
      typedef map<int,double>::const_iterator CI;

      cout << "\nAnalytical derivative for Intcell " << gpcounter << " (x-component). Slave GP " << gp << endl;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      cout << "Analytical derivative for Intcell " << gpcounter << " (y-component). Slave GP " << gp << endl;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;

      cout << "\nAnalytical derivative for Intcell " << gpcounter << " (x-component). Master GP " << gp << endl;
      for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      cout << "Analytical derivative for Intcell " << gpcounter << " (y-component). Master GP " << gp << endl;
      for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      
      cout << "\nAnalytical derivative for Intcell Jacobian Intcell " << gpcounter << " GP " << gp << endl;
      for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
    }
        
    // compute cell M matrix *********************************************
    // loop over all mseg matrix entries
    // nrow represents the slave Lagrange multipliers !!!
    // ncol represents the master dofs !!!
    // (this DOES matter here for mseg, as it might
    // sometimes be rectangular, not quadratic!)
    for (int j=0;j<nrow*ndof;++j)
    {
      for (int k=0;k<ncol*ndof;++k)
      {
        int jindex = (int)(j/ndof);
        int kindex = (int)(k/ndof);
              
        // multiply the two shape functions
        double prod = dualval[jindex]*mval[kindex];
       
        // isolate the mseg entries to be filled and
        // add current Gauss point's contribution to mseg
        if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          (*mseg)(j,k) += prod*jac*wgt; 
      }
    }   
    // compute cell M matrix *********************************************
    
    // compute cell gap vector *******************************************
    // loop over all gseg vector entries
    // nrow represents the slave side dofs !!!  */
    for (int j=0;j<nrow;++j)
    {
      double prod = dualval[j]*gap;
      // add current Gauss point's contribution to gseg  
      (*gseg)(j) += prod*jac*wgt; 
    }
    // compute cell gap vector *******************************************
    
    // compute cell M linearization **************************************
    for (int j=0;j<nrow;++j)
    {
      CONTACT::CNode* mycnode = static_cast<CONTACT::CNode*>(mynodes[j]);
      if (!mycnode) dserror("ERROR: IntegrateDerivCell3DAuxPlane: Null pointer!");
            
      for (int k=0;k<ncol;++k)
      {
        // global master node ID
        int mgid = mele.Nodes()[k]->Id();
        double fac = 0.0;
        
        // get the correct map as a reference
        map<int,double>& dmmap_jk = mycnode->GetDerivM()[mgid];
        
        // (1) Lin(Phi) - dual shape functions
        if (duallin)
          for (int m=0;m<nrow;++m)
          {
            fac = wgt*sval[m]*mval[k]*jac;
            for (CI p=dualmap[j][m].begin();p!=dualmap[j][m].end();++p)
              dmmap_jk[p->first] += fac*(p->second);
          }
        
        // (2) Lin(Phi) - slave GP coordinates
        fac = wgt*dualderiv(j,0)*mval[k]*jac;
        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        fac = wgt*dualderiv(j,1)*mval[k]*jac;
        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (3) Lin(NMaster) - master GP coordinates
        fac = wgt*dualval[j]*mderiv(k,0)*jac;
        for (CI p=dmxigp[0].begin();p!=dmxigp[0].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        fac = wgt*dualval[j]*mderiv(k,1)*jac;
        for (CI p=dmxigp[1].begin();p!=dmxigp[1].end();++p)
          dmmap_jk[p->first] += fac*(p->second);
        
        // (4) Lin(dsxideta) - intcell GP Jacobian
        fac = wgt*dualval[j]*mval[k];
        for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
          dmmap_jk[p->first] += fac*(p->second);
      }
    }
    // compute cell M linearization **************************************
    
    // compute cell gap linearization ************************************
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
      fac = wgt*dualderiv(j,0)*gap*jac;
      for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        dgmap[p->first] += fac*(p->second);
      fac = wgt*dualderiv(j,1)*gap*jac;
      for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        dgmap[p->first] += fac*(p->second);
      
      // (3) Lin(g) - gap function
      fac = wgt*dualval[j]*jac;
      for (CI p=dgapgp.begin();p!=dgapgp.end();++p)
        dgmap[p->first] += fac*(p->second);
      
      // (4) Lin(dsxideta) - intcell GP Jacobian
      fac = wgt*dualval[j]*gap;
      for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
        dgmap[p->first] += fac*(p->second);
    }
    // compute cell gap linearization ************************************ 
  }
  //**********************************************************************
  
  // finite difference checks
  testgps.push_back(testgpscurr);
  testgpm.push_back(testgpmcurr);
  testjs.push_back(testjscurr);
  testji.push_back(testjicurr);
  if (printderiv) gpcounter = gpcounter +1;
    
#ifdef CONTACTONEMORTARLOOP
    dserror("ERROR: IntegrateDerivCell3DAuxPlane: One mortar loop case not yet impl. for 3D!");
#endif // #ifdef CONTACTONEMORTARLOOP
    
  return;
}

#endif //#ifdef CCADISCRET
