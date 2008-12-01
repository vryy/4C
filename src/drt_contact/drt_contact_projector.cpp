/*!----------------------------------------------------------------------
\file drt_contact_projector.cpp
\brief A class to perform projections of nodes onto opposing CElements

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

#include "drt_contact_projector.H"
#include "drt_contact_interface.H"
#include "drt_celement.H"
#include "drt_cnode.H"
#include "contactdefines.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 01/08|
 *----------------------------------------------------------------------*/
CONTACT::Projector::Projector(int dim) :
dim_(dim)
{
  if (Dim()!=2 && Dim()!=3)
    dserror("ERROR: Contact problem must be 2D or 3D");
}

/*----------------------------------------------------------------------*
 |  Project a node along its nodal normal (public)            popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::ProjectNodalNormal(CONTACT::CNode& node,
                                            CONTACT::CElement& ele,
                                            double* xi)
{
  bool ok = true;
  if (Dim()==2)
  {
    // define variable to check if projection is outward w.r.t to slave
    bool outward = true;
    
    // define variable to check if inward projection is far off or nearby
    // (if it is nearby, this is a penetration case, which is feasible)
    double gap = 1.0e12;
    
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double deta = 0.0;
    double f = 0.0;
    double df = 0.0;
    int k=0;
    
    for (k=0;k<CONTACTMAXITER;++k)
    {
      f=EvaluateFNodalNormal(node,ele,eta,outward,gap);
      if (abs(f) < CONTACTCONVTOL) break;
      df=EvaluateGradFNodalNormal(node,ele,eta);
      deta=(-f)/df;
      eta[0]+=deta;
    }
    
    // get the result
    xi[0]=eta[0];
    
    // Newton iteration unconverged
    if (abs(f) > CONTACTCONVTOL)
    {
      ok = false;
      xi[0] = 9999.99;
      
      // Here (S->M projection) we only give a warning, no error!!!
      // This iteration sometimes diverges, when the projection is far off.
      // These cases are harmless, as these nodes then do not participate in
      // the overlap detection anyway!
      //cout << "***WARNING*** ProjectNodalNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and CElementID " << ele.Id() << endl;  
    }
    
    // no outward projection w.r.t to slave found    
    else if (!outward)
    {
      // we have to exclude the penetration case!
      // (there we also have an inward proj, but it is feasible of course!)
      if (gap>CONTACTCRITDIST)
      {
        ok = false;
        xi[0] = 9999.99;
      
        // At the moment we give a warning here, just to check!!!
        //cout << "***WARNING*** ProjectNodalNormal:" << " Inward projection for NodeID "
        //     << node.Id() << " and CElementID " << ele.Id() << endl;  
      }
    }
/*    
#ifdef DEBUG      
    // Newton iteration converged
    else
    {
      cout << "Newton iteration converged in " << k << " step(s)!" << endl
           << "The result is: " << xi[0] << endl;
    }
#endif // #ifdef DEBUG
*/    
  } // if (Dim()==2)
  
  else
    dserror("ERROR: ProjectNodalNormal: Called 2D version for 3D problem!");
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Project a node along element's normal field (public)      popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::ProjectElementNormal(CONTACT::CNode& node,
                                              CONTACT::CElement& ele,
                                              double* xi)
{
  bool ok = true;
  
  if (Dim()==2)
  {    
    // define variable to check if projection is outward w.r.t to slave
    bool outward = true;
    
    // define variable to check if inward projection is far off or nearby
    // (if it is nearby, this is a penetration case, which is feasible)
    double gap = 1.0e12;
        
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double deta = 0.0;
    double f = 0.0;
    double df = 0.0;
    int k=0;
    
    for (k=0;k<CONTACTMAXITER;++k)
    {
      f=EvaluateFElementNormal(node,ele,eta,outward,gap);
      if (abs(f) < CONTACTCONVTOL) break;
      df=EvaluateGradFElementNormal(node,ele,eta);
      deta=(-f)/df;
      eta[0]+=deta;
    }
    
    // get the result
    xi[0]=eta[0];
    
    // Newton iteration unconverged
    if (abs(f) > CONTACTCONVTOL)
    {
      ok = false;
      xi[0] = 9999.99;
      
      // Here (M->S projection) we only give a warning, no error!!!
      // This iteration sometimes diverges, when the projection is far off.
      // These cases are harmless, as these nodes then do not participate in
      // the overlap detection anyway!
      //cout << "***WARNING*** ProjectElementNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and CElementID " << ele.Id() << endl;
    }
    
    // no outward projection w.r.t to slave found
    else if (!outward)
    {
      // we have to exclude the penetration case!
      // (there we also have an inward proj, but it is feasible of course!)
      if (gap>CONTACTCRITDIST)
      {
        ok = false;
        xi[0] = 9999.99;
      
        // At the moment we give a warning here, just to check!!!
        //cout << "***WARNING*** ProjectElementNormal:" << " Inward projection for NodeID "
        //      << node.Id() << " and CElementID " << ele.Id() << endl;  
      }
    }
/*    
#ifdef DEBUG      
    // Newton iteration converged
    else
    {
      cout << "Newton iteration converged in " << k << " step(s)!" << endl
           << "The result is: " << xi[0] << endl;
    }
#endif // #ifdef DEBUG
*/    
  } // if (Dim()==2)
    
  else
    dserror("ERROR: ProjectElementNormal: Called 2D version for 3D problem!");
  
  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a node along element's normal field (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::ProjectElementNormal3D(CONTACT::CNode& node,
                                                CONTACT::CElement& ele,
                                                double* xi)
{
  // start in the element center
  DRT::Element::DiscretizationType dt = ele.Shape();
  double eta[2] = {0.0, 0.0};
  if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
  {
    eta[0] = 1.0/3;
    eta[1] = 1.0/3;
  }
  
  // coordinate increment
  double deta[2] = {0.0, 0.0};
  
  // auxiliary variable
  double alpha = 0.0;
  double dalpha = 0.0;
  
  // function f (vector-valued)
  double f[3] = {0.0, 0.0, 0.0};
  
  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  Epetra_SerialDenseMatrix df(3,3);
  
  // start iteration
  int k=0;
  double conv = 0.0;
  
  for (k=0;k<CONTACTMAXITER;++k)
  {
    EvaluateFElementNormal3D(f,node,ele,eta,alpha);
    conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    //cout << "Iteration " << k << ": -> |f|=" << conv << endl;
    if (conv <= CONTACTCONVTOL) break;
    EvaluateGradFElementNormal3D(df,node,ele,eta,alpha);
    
    // solve deta = - inv(df) * f
    LINALG::NonSymmetricInverse(df,3);
    deta[0] = -df(0,0)*f[0] - df(0,1)*f[1] - df(0,2)*f[2];
    deta[1] = -df(1,0)*f[0] - df(1,1)*f[1] - df(1,2)*f[2];
    dalpha  = -df(2,0)*f[0] - df(2,1)*f[1] - df(2,2)*f[2];
    
    // update eta and alpha
    eta[0] += deta[0];
    eta[1] += deta[1];
    alpha  += dalpha;
  }
      
  // Newton iteration unconverged
  if (conv > CONTACTCONVTOL)
    dserror("ERROR: ProjectElementNormal3D: Newton unconverged for NodeID %i "
            "and CElementID %i", node.Id(), ele.Id());

  // Newton iteration converged
  xi[0]=eta[0];
  xi[1]=eta[1];
  //cout << "Newton iteration converged in " << k << " steps!" << endl;
  
  return true;  
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (public)           popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::ProjectGaussPoint(CONTACT::CElement& gpele,
                                           const double* gpeta,
                                           CONTACT::CElement& ele,
                                           double* xi)
{
  bool ok = true;
  if (Dim()==2)
  {
    // collect necessary data (slave side, for GP)
    int nnodes = gpele.NumNode();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
    LINALG::SerialDenseMatrix coord(3,nnodes);
    DRT::Node** mynodes = gpele.Nodes();
    if(!mynodes) dserror("ERROR: ProjectGaussPoint: Null pointer!");
        
    // get shape function values and derivatives at gpeta
    ele.EvaluateShape(gpeta, val, deriv, nnodes);

    // get interpolated GP normal and GP coordinates
    double gpn[3] = {0.0, 0.0, 0.0};
    double gpx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nnodes;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      
      gpn[0]+=val[i]*mycnode->n()[0];
      gpn[1]+=val[i]*mycnode->n()[1];
      gpn[2]+=val[i]*mycnode->n()[2];
      
      coord(0,i) = mycnode->xspatial()[0];
      coord(1,i) = mycnode->xspatial()[1];
      coord(2,i) = mycnode->xspatial()[2];
      
      gpx[0]+=val[i]*coord(0,i);
      gpx[1]+=val[i]*coord(1,i);
      gpx[2]+=val[i]*coord(2,i);
    }
    
    // define variable to check if projection is outward w.r.t to slave
    bool outward = true;
    
    // define variable to check if inward projection is far off or nearby
    // (if it is nearby, this is a penetration case, which is feasible)
    double gap = 1.0e12;
        
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double deta = 0.0;
    double f = 0.0;
    double df = 0.0;
    int k=0;
    
    for (k=0;k<CONTACTMAXITER;++k)
    {
      f=EvaluateFGaussPoint(gpx,gpn,ele,eta,outward,gap);
      if (abs(f) < CONTACTCONVTOL) break;
      df=EvaluateGradFGaussPoint(gpn,ele,eta);
      deta=(-f)/df;
      eta[0]+=deta;
    }
        
    // Newton iteration unconverged
    if (abs(f) > CONTACTCONVTOL)
    {
      ok = false;
      dserror("ERROR: ProjectGaussPoint: Newton unconverged for GP at xi=%d"
              " from CElementID %i", gpeta[0], gpele.Id());
    }
    
    // no outward projection w.r.t to slave found    
    else if (!outward)
    {
      // we have to exclude the penetration case!
      // (there we also have an inward proj, but it is feasible of course!)
      if (gap>CONTACTCRITDIST)
      {
      ok = false;
      dserror("ERROR: ProjectGaussPoint: Inward projection for GP at xi=%d"
              " from CElementID %i", gpeta[0], gpele.Id());
      }
    }
    
    // Newton iteration converged
    xi[0]=eta[0];
/*
#ifdef DEBUG
    cout << "GP Newton iteration converged in " << k << " step(s)!" << endl
         << "The result is: " << xi[0] << endl;
#endif // #ifdef DEBUG
*/
  } // if (Dim()==2)
  
  else
    dserror("ERROR: ProjectGaussPoint: Called 2D version for 3D problem!");
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (3D)               popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::ProjectGaussPoint3D(CONTACT::CElement& gpele,
                                             const double* gpeta,
                                             CONTACT::CElement& ele,
                                             double* xi)
{
  // collect necessary data (slave side, for GP)
  int nnodes = gpele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
  LINALG::SerialDenseMatrix coord(3,nnodes);
  DRT::Node** mynodes = gpele.Nodes();
  if(!mynodes) dserror("ERROR: ProjectGaussPoint3D: Null pointer!");
      
  // get shape function values and derivatives at gpeta
  gpele.EvaluateShape(gpeta, val, deriv, nnodes);

  // get interpolated GP normal and GP coordinates
  double gpn[3] = {0.0, 0.0, 0.0};
  double gpx[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    
    gpn[0]+=val[i]*mycnode->n()[0];
    gpn[1]+=val[i]*mycnode->n()[1];
    gpn[2]+=val[i]*mycnode->n()[2];
    
    coord(0,i) = mycnode->xspatial()[0];
    coord(1,i) = mycnode->xspatial()[1];
    coord(2,i) = mycnode->xspatial()[2];
    
    gpx[0]+=val[i]*coord(0,i);
    gpx[1]+=val[i]*coord(1,i);
    gpx[2]+=val[i]*coord(2,i);
  }
    
  // start in the element center
  DRT::Element::DiscretizationType dt = ele.Shape();
  double eta[2] = {0.0, 0.0};
  if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
  {
    eta[0] = 1.0/3;
    eta[1] = 1.0/3;
  }
  
  // coordinate increment
  double deta[2] = {0.0, 0.0};
  
  // auxiliary variable
  double alpha = 0.0;
  double dalpha = 0.0;
  
  // function f (vector-valued)
  double f[3] = {0.0, 0.0, 0.0};
  
  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  Epetra_SerialDenseMatrix df(3,3);
  
  // start iteration
  int k=0;
  double conv = 0.0;
  
  for (k=0;k<CONTACTMAXITER;++k)
  {
    EvaluateFGaussPoint3D(f,gpx,gpn,ele,eta,alpha);
    conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    //cout << "Iteration " << k << ": -> |f|=" << conv << endl;
    if (conv <= CONTACTCONVTOL) break;
    EvaluateGradFGaussPoint3D(df,gpx,gpn,ele,eta,alpha);
    
    // solve deta = - inv(df) * f
    LINALG::NonSymmetricInverse(df,3);
    deta[0] = -df(0,0)*f[0] - df(0,1)*f[1] - df(0,2)*f[2];
    deta[1] = -df(1,0)*f[0] - df(1,1)*f[1] - df(1,2)*f[2];
    dalpha  = -df(2,0)*f[0] - df(2,1)*f[1] - df(2,2)*f[2];
    
    // update eta and alpha
    eta[0] += deta[0];
    eta[1] += deta[1];
    alpha  += dalpha;
  }
      
  // Newton iteration unconverged
  if (conv > CONTACTCONVTOL)
    dserror("ERROR: ProjectGaussPoint3D: Newton unconverged for GP at xi = %d,"
                  " %d from CElementID %i", gpeta[0], gpeta[1], gpele.Id());

  // Newton iteration converged
  xi[0]=eta[0];
  xi[1]=eta[1];
  
  //cout << "Newton iteration converged in " << k << " steps!" << endl;
  
  return true;
}
/*----------------------------------------------------------------------*
 |  Project a Gauss point along AuxPlane normal (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::ProjectGaussPointAuxn3D(const double* globgp,
                                                 const double* auxn,
                                                 CONTACT::CElement& ele,
                                                 double* xi)
{
  // start in the element center
  DRT::Element::DiscretizationType dt = ele.Shape();
  double eta[2] = {0.0, 0.0};
  if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
  {
    eta[0] = 1.0/3;
    eta[1] = 1.0/3;
  }
  
  // coordinate increment
  double deta[2] = {0.0, 0.0};
  
  // auxiliary variable
  double alpha = 0.0;
  double dalpha = 0.0;
  
  // function f (vector-valued)
  double f[3] = {0.0, 0.0, 0.0};
  
  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  Epetra_SerialDenseMatrix df(3,3);
  
  // start iteration
  int k=0;
  double conv = 0.0;
  
  for (k=0;k<CONTACTMAXITER;++k)
  {
    EvaluateFGaussPointAuxn3D(f,globgp,auxn,ele,eta,alpha);
    conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    //cout << "Iteration " << k << ": -> |f|=" << conv << endl;
    if (conv <= CONTACTCONVTOL) break;
    EvaluateGradFGaussPointAuxn3D(df,globgp,auxn,ele,eta,alpha);
    
    // solve deta = - inv(df) * f
    LINALG::NonSymmetricInverse(df,3);
    deta[0] = -df(0,0)*f[0] - df(0,1)*f[1] - df(0,2)*f[2];
    deta[1] = -df(1,0)*f[0] - df(1,1)*f[1] - df(1,2)*f[2];
    dalpha  = -df(2,0)*f[0] - df(2,1)*f[1] - df(2,2)*f[2];
    
    // update eta and alpha
    eta[0] += deta[0];
    eta[1] += deta[1];
    alpha  += dalpha;
  }
      
  // Newton iteration unconverged
  if (conv > CONTACTCONVTOL)
    dserror("ERROR: ProjectGaussPointAuxn3D: Newton unconverged for GP"
            "at xi = (%f,%f,%f) onto CElementID %i", globgp[0],globgp[1],globgp[2],ele.Id());

  // Newton iteration converged
  xi[0]=eta[0];
  xi[1]=eta[1];
  
  //cout << "Newton iteration converged in " << k << " steps!" << endl;
  
  return true;  
}

/*----------------------------------------------------------------------*
 |  Evaluate F for nodal normal case (public)                 popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::EvaluateFNodalNormal(CONTACT::CNode& node,
                                                CONTACT::CElement& ele,
                                                const double* eta,
                                                bool& outward,
                                                double& gap)
{
  /* Evaluate the function F(eta) = ( Ni * xim - xs ) x ns,
     or to be more precise the third component of this vector function!
    
       Ni  shape functions of element
       xim coords of element nodes (master side)
       xs  coords of node to be projected (slave side)
       ns   outward normal of node to be projected (slave side)          */
  double fval = 0.0;
  
  if (Dim()==2)
  {
    // collect necessary data (master side)
    int nnodes = ele.NumNode();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
    
    // build interpolation of master node coordinates for current eta
    double nx[3] = {0.0, 0.0, 0.0};
    ele.LocalToGlobal(eta,nx,0);
    
    // subtract slave node coordinates
    nx[0]-=node.xspatial()[0];
    nx[1]-=node.xspatial()[1];
    nx[2]-=node.xspatial()[2];
  
    // update boolean variable outward and gap
    gap = nx[0]*node.n()[0]+nx[1]*node.n()[1]+nx[2]*node.n()[2];
    if (gap<0.0)
      outward = false;
    else
      outward = true;
    gap=abs(gap);
    
    //calculate F
    fval = nx[0]*node.n()[1]-nx[1]*node.n()[0];
  }
  
  else
    dserror("ERROR: EvaluateFNodalNormal: 3D version not yet implemented!");
  
  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for nodal normal case (public)             popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::EvaluateGradFNodalNormal(CONTACT::CNode& node,
                                                    CONTACT::CElement& ele,
                                                    const double* eta)
{
  /* Evaluate the function GradF(eta)
     = Ni,eta * xim * nys - Ni,eta * yim * nxs,
    
       Ni,eta    shape function derivatives of element
       xim, yim  coords of element nodes (master side)
       nxs, nys   outward normal of node to be projected (slave side)   */
  
  double fgrad = 0.0;
    
  if (Dim()==2)
  {
    // collect necessary data (master side)
    int nnodes = ele.NumNode();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
        
    // build interpolation of master node coordinates for current eta
    // use shape function derivatives for interpolation (hence "1")
    double nxeta[3] = {0.0, 0.0, 0.0};
    ele.LocalToGlobal(eta,nxeta,1);
    
    // calculate GradF
    fgrad =  nxeta[0]*node.n()[1]-nxeta[1]*node.n()[0];
  }
  
  else
    dserror("ERROR: EvaluateGradFNodalNormal: 3D version not yet implemented!");
    
  return fgrad;
    
}

/*----------------------------------------------------------------------*
 |  Evaluate F for element normal case (public)               popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::EvaluateFElementNormal(CONTACT::CNode& node,
                                                  CONTACT::CElement& ele,
                                                  const double* eta,
                                                  bool& outward,
                                                  double& gap)
{
  /* Evaluate the function F(eta) = ( Ni * xis - xm ) x ( Nj * njs),
     or to be more precise the third component of this vector function!
      
       Ni  shape functions of element
       xis coords of element nodes (slave side)
       xm  coords of node to be projected (master side)
       nis outward normals of element nodes (slave side)                */
  
  double fval = 0.0;
    
  if (Dim()==2)
  {
    // collect necessary data (slave side)
    int nnodes = ele.NumNode();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
    LINALG::SerialDenseMatrix coord(3,nnodes);
    DRT::Node** mynodes = ele.Nodes();
    if(!mynodes) dserror("ERROR: EvaluateFElementNormal: Null pointer!");
        
    // get shape function values and derivatives at eta
    ele.EvaluateShape(eta, val, deriv, nnodes);
  
    // get interpolated normal and proj. coordinates for current eta
    double nn[3] = {0.0, 0.0, 0.0};
    double nx[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nnodes;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      nn[0]+=val[i]*mycnode->n()[0];
      nn[1]+=val[i]*mycnode->n()[1];
      nn[2]+=val[i]*mycnode->n()[2];
      
      coord(0,i) = mycnode->xspatial()[0];
      coord(1,i) = mycnode->xspatial()[1];
      coord(2,i) = mycnode->xspatial()[2];
      
      nx[0]+=val[i]*coord(0,i);
      nx[1]+=val[i]*coord(1,i);
      nx[2]+=val[i]*coord(2,i);
    }
      
    // subtract master node coordinates
    nx[0]-=node.xspatial()[0];
    nx[1]-=node.xspatial()[1];
    nx[2]-=node.xspatial()[2];
    
    // update boolean variable outward
    gap = -nx[0]*nn[0]-nx[1]*nn[1]-nx[2]*nn[2];
    if (gap<0.0)
      outward = false;
    else
      outward = true;
    gap=abs(gap);
    
    // calculate F
    fval = nx[0]*nn[1]-nx[1]*nn[0];
  }
    
  else
    dserror("ERROR: EvaluateFElementNormal: 3D version not yet implemented!");
      
  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for element normal case (public)           popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::EvaluateGradFElementNormal(CONTACT::CNode& node,
                                                      CONTACT::CElement& ele,
                                                      const double* eta)
{
  /* Evaluate the function GradF(eta)
      = ( Ni,eta * xis ) * ( Nj * nyjs )
      + ( Ni * xis - xm ) * ( Nj,eta * nyjs )
      - ( Ni,eta * yis ) * ( Nj * nxjs )
      - ( Ni * yis - ym ) * ( Nj,eta * nxjs )
        
       Ni,eta      shape function derivatives of element
       xis, yis   coords of element nodes (slave side)
       xm, ym     coords of node to be projected (master side)
       nxjs, nyjs outward normals of element nodes (slave side)         */
  
  double fgrad = 0.0;
      
  if (Dim()==2)
  {
    // collect necessary data (slave side)
    int nnodes = ele.NumNode();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
    LINALG::SerialDenseMatrix coord(3,nnodes);
    DRT::Node** mynodes = ele.Nodes();
    if(!mynodes) dserror("ERROR: EvaluateGradFElementNormal: Null pointer!");
          
    // get shape function values and derivatives at eta
    ele.EvaluateShape(eta, val, deriv, nnodes);
  
    // get interpolated normal and proj. coordinates for current eta
    double nn[3] = {0.0, 0.0, 0.0};
    double nneta[3] = {0.0, 0.0, 0.0};
    double nx[3] = {0.0, 0.0, 0.0};
    double nxeta[3] = {0.0, 0.0, 0.0};
    for (int i=0;i<nnodes;++i)
    {
      CNode* mycnode = static_cast<CNode*> (mynodes[i]);
      
      nn[0]+=val[i]*mycnode->n()[0];
      nn[1]+=val[i]*mycnode->n()[1];
      nn[2]+=val[i]*mycnode->n()[2];
      
      nneta[0]+=deriv(i,0)*mycnode->n()[0];
      nneta[1]+=deriv(i,0)*mycnode->n()[1];
      nneta[2]+=deriv(i,0)*mycnode->n()[2];
        
      coord(0,i) = mycnode->xspatial()[0];
      coord(1,i) = mycnode->xspatial()[1];
      coord(2,i) = mycnode->xspatial()[2];
      
      nx[0]+=val[i]*coord(0,i);
      nx[1]+=val[i]*coord(1,i);
      nx[2]+=val[i]*coord(2,i);
          
      nxeta[0]+=deriv(i,0)*coord(0,i);
      nxeta[1]+=deriv(i,0)*coord(1,i);
      nxeta[2]+=deriv(i,0)*coord(2,i);
    }
        
    // subtract master node coordinates
    nx[0]-=node.xspatial()[0];
    nx[1]-=node.xspatial()[1];
    nx[2]-=node.xspatial()[2];
      
    // calculate GradF
    fgrad =   nxeta[0]*nn[1] + nx[0]*nneta[1]
            - nxeta[1]*nn[0] - nx[1]*nneta[0];
  }
  
  else
      dserror("ERROR: EvaluateGradFElementNormal: 3D version not yet implemented!");
    
  return fgrad;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for element normal case (3D)                   popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::EvaluateFElementNormal3D(
                              double* f,
                              CONTACT::CNode& node,
                              CONTACT::CElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi + alpha * (Ni * ni) -xm
     which is a vector-valued function with 3 components!
    
       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       ni      outward normal of nodes of element to project on
       xm      node to project (master)                                 */

  // collect necessary data
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
  LINALG::SerialDenseMatrix coord(3,nnodes);
  DRT::Node** mynodes = ele.Nodes();
  if(!mynodes) dserror("ERROR: EvaluateFElementNormal3D: Null pointer!");
      
  // get shape function values and derivatives at eta
  ele.EvaluateShape(eta, val, deriv, nnodes);

  // get interpolated normal and proj. coordinates for current eta
  double nn[3] = {0.0, 0.0, 0.0};
  double nx[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    nn[0]+=val[i]*mycnode->n()[0];
    nn[1]+=val[i]*mycnode->n()[1];
    nn[2]+=val[i]*mycnode->n()[2];
    
    coord(0,i) = mycnode->xspatial()[0];
    coord(1,i) = mycnode->xspatial()[1];
    coord(2,i) = mycnode->xspatial()[2];
    
    nx[0]+=val[i]*coord(0,i);
    nx[1]+=val[i]*coord(1,i);
    nx[2]+=val[i]*coord(2,i);
  }
  
  // get coords of node to be projected
  double xm[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<3;++i)
    xm[i] = node.xspatial()[i];
  
  // evaluate function f
  for (int i=0;i<3;++i)
    f[i] = nx[i] + alpha * nn[i] - xm[i];
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for element normal case (3D)               popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::EvaluateGradFElementNormal3D(
                              Epetra_SerialDenseMatrix& fgrad,
                              CONTACT::CNode& node,
                              CONTACT::CElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
     + alpha * (Ni * ni) - xm, which is a (3x3)-matrix!
    
       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       ni      outward normal of nodes of element to project on
       xm      node to project (master)                                 */
  
  // collect necessary data
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
  LINALG::SerialDenseMatrix coord(3,nnodes);
  DRT::Node** mynodes = ele.Nodes();
  if(!mynodes) dserror("ERROR: EvaluateGradFElementNormal3D: Null pointer!");
      
  // get shape function values and derivatives at eta
  ele.EvaluateShape(eta, val, deriv, nnodes);

  // get interpolated normal + deriv and proj. coordinates deriv for current eta
  double nn[3] = {0.0, 0.0, 0.0};
  double nneta1[3] = {0.0, 0.0, 0.0};
  double nneta2[3] = {0.0, 0.0, 0.0};
  double nxeta1[3] = {0.0, 0.0, 0.0};
  double nxeta2[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    nn[0]+=val[i]*mycnode->n()[0];
    nn[1]+=val[i]*mycnode->n()[1];
    nn[2]+=val[i]*mycnode->n()[2];
        
    nneta1[0]+=deriv(i,0)*mycnode->n()[0];
    nneta1[1]+=deriv(i,0)*mycnode->n()[1];
    nneta1[2]+=deriv(i,0)*mycnode->n()[2];
    
    nneta2[0]+=deriv(i,1)*mycnode->n()[0];
    nneta2[1]+=deriv(i,1)*mycnode->n()[1];
    nneta2[2]+=deriv(i,1)*mycnode->n()[2];
    
    coord(0,i) = mycnode->xspatial()[0];
    coord(1,i) = mycnode->xspatial()[1];
    coord(2,i) = mycnode->xspatial()[2];
    
    nxeta1[0]+=deriv(i,0)*coord(0,i);
    nxeta1[1]+=deriv(i,0)*coord(1,i);
    nxeta1[2]+=deriv(i,0)*coord(2,i);
    
    nxeta2[0]+=deriv(i,1)*coord(0,i);
    nxeta2[1]+=deriv(i,1)*coord(1,i);
    nxeta2[2]+=deriv(i,1)*coord(2,i);
  }
  
  //evaluate function f gradient
  for (int i=0;i<3;++i)
  {
    fgrad(i,0) = nxeta1[i] + alpha*nneta1[i];
    fgrad(i,1) = nxeta2[i] + alpha*nneta2[i];
    fgrad(i,2) = nn[i];
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (public)                  popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::EvaluateFGaussPoint(const double* gpx,
                                               const double* gpn,
                                               CONTACT::CElement& ele,
                                               const double* eta,
                                               bool& outward,
                                               double& gap)
{
  /* Evaluate the function F(eta) = ( Ni * xim - gpx ) x gpn,
     or to be more precise the third component of this vector function!
    
       Ni  shape functions of element (master side)
       xim coords of element nodes (master side)
       gpx coords of GP to be projected (slave side)
       gpn outward normal of GP to be projected (slave side)          */
  
  double fval = 0.0;
      
  if (Dim()==2)
  {
    // collect necessary data (master side)
    int nnodes = ele.NumNode();
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,1);
    
    // build interpolation of master node coordinates for current eta
    double nx[3] = {0.0, 0.0, 0.0};
    ele.LocalToGlobal(eta,nx,0);
    
    // subtract GP coordinates
    nx[0]-=gpx[0];
    nx[1]-=gpx[1];
    nx[2]-=gpx[2];
  
    // update boolean variable outward
    gap = nx[0]*gpn[0]+nx[1]*gpn[1]+nx[2]*gpn[2];
    if (gap<0.0)
      outward = false;
    else
      outward = true;
    gap=abs(gap);
    
    // calculate F
    fval = nx[0]*gpn[1]-nx[1]*gpn[0];
  }
  
  else
    dserror("ERROR: EvaluateFGaussPoint: 3D version not yet implemented!");
        
  return fval;
  
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (public)              popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::EvaluateGradFGaussPoint(const double* gpn,
                                                      CONTACT::CElement& ele,
                                                      const double* eta)
{
  /* Evaluate the function GradF(eta)
      = Ni,eta * xim * gpny - Ni,eta * yim * gpnx,
    
       Ni,eta     shape function derivatives of element (master side)
       xim, yim   coords of element nodes (master side)
       gpnx, gpny outward normal of GP to be projected (slave side)   */
  
  double fgrad = 0.0;
        
  if (Dim()==2)
  {
  // collect necessary data (master side)
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,1);
  
  // build interpolation of master node coordinates for current eta
  // use shape function derivatives for interpolation (hence "1")
  double nxeta[3] = {0.0, 0.0, 0.0};
  ele.LocalToGlobal(eta,nxeta,1);

  // calculate GradF
  fgrad = nxeta[0]*gpn[1]-nxeta[1]*gpn[0];
  }
  
  else
    dserror("ERROR: EvaluateGradFGaussPoint: 3D version not yet implemented!");
      
  return fgrad;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (3D)                      popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::EvaluateFGaussPoint3D(
                              double* f,
                              const double* gpx,
                              const double* gpn,
                              CONTACT::CElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi - alpha * gpn - gpx
     which is a vector-valued function with 3 components!
    
       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       gpx     coords of GP to be projected
       gpn     normal of GP along which to project                  */

  // collect necessary data
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
    
  // build interpolation of ele node coordinates for current eta
  double nx[3] = {0.0, 0.0, 0.0};
  ele.LocalToGlobal(eta,nx,0);

  // evaluate function f
  for (int i=0; i<3; ++i)
    f[i] = nx[i] - alpha * gpn[i] - gpx[i];
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (3D)                  popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::EvaluateGradFGaussPoint3D(
                              Epetra_SerialDenseMatrix& fgrad,
                              const double* gpx,
                              const double* gpn,
                              CONTACT::CElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
     - alpha * gpn - gpx, which is a (3x3)-matrix!
    
       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       gpx     coords of GP to be projected
       gpn     normal of GP along which to project                  */
  
  // collect necessary data
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
    
  // build interpolation of ele node coordinates for current eta
  double nxeta1[3] = {0.0, 0.0, 0.0};
  double nxeta2[3] = {0.0, 0.0, 0.0};
  ele.LocalToGlobal(eta,nxeta1,1);
  ele.LocalToGlobal(eta,nxeta2,2);
  
  //evaluate function f gradient
  for (int i=0;i<3;++i)
  {
    fgrad(i,0) = nxeta1[i];
    fgrad(i,1) = nxeta2[i];
    fgrad(i,2) = -gpn[i];
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for AuxPlane Gauss point case (3D)             popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::EvaluateFGaussPointAuxn3D(
                              double* f,
                              const double* globgp,
                              const double* auxn,
                              CONTACT::CElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi - alpha * auxn - globgp
     which is a vector-valued function with 3 components!
    
       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       globgp  coords of AuxPlaneGP to be projected
       auxn    normal of AuxPlane along which to project            */

  // collect necessary data
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
    
  // build interpolation of ele node coordinates for current eta
  double nx[3] = {0.0, 0.0, 0.0};
  ele.LocalToGlobal(eta,nx,0);

  // evaluate function f
  for (int i=0; i<3; ++i)
    f[i] = nx[i] - alpha * auxn[i] - globgp[i];
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for AuxPlane Gauss point case (3D)         popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::EvaluateGradFGaussPointAuxn3D(
                              Epetra_SerialDenseMatrix& fgrad,
                              const double* globgp,
                              const double* auxn,
                              CONTACT::CElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
     - alpha * auxn - globgp, which is a (3x3)-matrix!
    
       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       globgp  coords of AuxPlaneGP to be projected
       auxn    normal of AuxPlane along which to project            */
  
  // collect necessary data
  int nnodes = ele.NumNode();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
    
  // build interpolation of ele node coordinates for current eta
  double nxeta1[3] = {0.0, 0.0, 0.0};
  double nxeta2[3] = {0.0, 0.0, 0.0};
  ele.LocalToGlobal(eta,nxeta1,1);
  ele.LocalToGlobal(eta,nxeta2,2);
  
  //evaluate function f gradient
  for (int i=0;i<3;++i)
  {
    fgrad(i,0) = nxeta1[i];
    fgrad(i,1) = nxeta2[i];
    fgrad(i,2) = -auxn[i];
  }

  return true;
}

#endif //#ifdef CCADISCRET
