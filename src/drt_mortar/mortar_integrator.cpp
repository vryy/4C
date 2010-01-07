/*!----------------------------------------------------------------------
\file mortar_integrator.cpp
\brief A class to perform integrations of Mortar matrices on the overlap
 			 of two MortarElements in 1D and 2D

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

#include "mortar_integrator.H"
#include "mortar_element.H"
#include "mortar_projector.H"
#include "mortar_defines.H"
#include "../drt_fem_general/drt_utils_integration.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 08/08|
 *----------------------------------------------------------------------*/
MORTAR::Integrator::Integrator(DRT::Element::DiscretizationType eletype) :
shapefcn_(MortarInterface::Undefined)
{
  InitializeGP(eletype);
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 07/09|
 *----------------------------------------------------------------------*/
MORTAR::Integrator::Integrator(const MortarInterface::ShapeFcnType shapefcn,
                                DRT::Element::DiscretizationType eletype) :
shapefcn_(shapefcn)
{
  InitializeGP(eletype);
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
void MORTAR::Integrator::InitializeGP(DRT::Element::DiscretizationType eletype)
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
void MORTAR::Integrator::IntegrateDerivSlave2D3D(
      MORTAR::MortarElement& sele, double* sxia, double* sxib,
      RCP<Epetra_SerialDenseMatrix> dseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == MortarInterface::Undefined)
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
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);

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
    if (shapefcn_ == MortarInterface::DualFunctions)
      sele.EvaluateShapeDual(eta,dualval,dualderiv,nrow);

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
        double prod = 0.0;
        if (shapefcn_ == MortarInterface::DualFunctions)
          prod = dualval[jindex]*val[kindex];
        else if (shapefcn_ == MortarInterface::StandardFunctions)
          prod =val[jindex]*val[kindex];

        // isolate the dseg entries to be filled
        // (both the main diagonal and every other secondary diagonal)
        // and add current Gauss point's contribution to dseg
        if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          (*dseg)(j,k) += prod*dxdsxi*wgt;
      }
    }
    // compute element D matrix ******************************************
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
void MORTAR::Integrator::IntegrateDerivSegment2D(
     MORTAR::MortarElement& sele, double& sxia, double& sxib,
     MORTAR::MortarElement& mele, double& mxia, double& mxib,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == MortarInterface::Undefined)
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
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,1);
  
  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

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
    MORTAR::Projector projector(2);
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
    if (shapefcn_ == MortarInterface::DualFunctions)
      sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    // reflect the MORTARONELOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef MORTARONELOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef MORTARONELOOP
    
    // decide whether boundary modification has to be considered or not
    // this is element-specific (is there a boundary node in this element?)
    bool bound = false;
#ifdef MORTARBOUNDMOD
    for (int k=0;k<nrow;++k)
    {
      // check the shape function type (not really necessary because only dual shape functions arrive here)
      if (shapefcn_ == MortarInterface::StandardFunctions)
        dserror("ERROR: IntegrateAndDerivSlave: Edge node mod. called for standard shape functions");
      
      MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
      if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
      bound += mymrtrnode->IsOnBound();
    }
#endif // #ifdef MORTARBOUNDMOD

    // compute segment D/M matrix ****************************************
    
    if (shapefcn_ == MortarInterface::StandardFunctions)
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
    
    else if (shapefcn_ == MortarInterface::DualFunctions)
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
  }
  //**********************************************************************

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
RCP<Epetra_SerialDenseMatrix> MORTAR::Integrator::IntegrateMmod2D(MORTAR::MortarElement& sele,
                                                                 double& sxia, double& sxib,
                                                                 MORTAR::MortarElement& mele,
                                                                 double& mxia, double& mxib)
{
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateMmod2D called on a wrong type of MortarElement pair!");
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

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    MORTAR::Projector projector(2);
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
  MortarNode* snode0 = static_cast<MortarNode*>(sele.Nodes()[0]);
  MortarNode* snode1 = static_cast<MortarNode*>(sele.Nodes()[1]);

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
void MORTAR::Integrator::IntegrateDerivCell3D(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     RCP<MORTAR::Intcell> cell,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == MortarInterface::Undefined)
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
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);

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
    MORTAR::Projector projector(3);
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
    if(shapefcn_ == MortarInterface::DualFunctions)
      sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the two Jacobians (int. cell and slave element)
    double jaccell = cell->Jacobian(eta);
    double jacslave = sele.Jacobian(sxi);

    // reflect the MORTARONELOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef MORTARONELOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef MORTARONELOOP

    // compute cell D/M matrix *******************************************
    
    if (shapefcn_ == MortarInterface::StandardFunctions)
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
    
    else if (shapefcn_ == MortarInterface::DualFunctions)
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
void MORTAR::Integrator::IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     RCP<MORTAR::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == MortarInterface::Undefined)
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

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = Dim();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);

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
    MORTAR::Projector projector(3);
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
    if (shapefcn_ == MortarInterface::DualFunctions)
      sele.EvaluateShapeDual(sxi,dualval,dualderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // reflect the MORTARONELOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef MORTARONELOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef MORTARONELOOP

    // compute cell D/M matrix *******************************************
    if (shapefcn_ == MortarInterface::StandardFunctions)
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
    
    else if (shapefcn_ == MortarInterface::DualFunctions)
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
void MORTAR::Integrator::IntegrateDerivCell3DAuxPlaneQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     RCP<MORTAR::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg)
{
  
  // explicitely defined shapefunction type needed
  if (shapefcn_ == MortarInterface::Undefined)
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
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector dualval(nrow);
  LINALG::SerialDenseMatrix dualderiv(nrow,2,true);
  LINALG::SerialDenseVector sintval(nintrow);
  LINALG::SerialDenseMatrix sintderiv(nintrow,2,true);
  LINALG::SerialDenseVector dualintval(nintrow);
  LINALG::SerialDenseMatrix dualintderiv(nintrow,2,true);

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
    MORTAR::Projector projector(3);
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

    // evaluate dual space shape functions (on slave element)
    if (shapefcn_ == MortarInterface::DualFunctions)
    {
      sele.EvaluateShapeDual(psxi,dualval,dualderiv,nrow);
      sintele.EvaluateShapeDual(sxi,dualintval,dualintderiv,nintrow);
    }
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(psxi,sval,sderiv,nrow);
    mele.EvaluateShape(pmxi,mval,mderiv,ncol);
    sintele.EvaluateShape(sxi,sintval,sintderiv,nintrow);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // reflect the MORTARONELOOP flag:
    // decide whether D and LinD have to be integrated or not
    
    // this is the standard case
    // integration of D and LinD is done separately in IntegrateDerivSlave2D3D
    // (dissimilar ways of computing the entries is employed)
    bool dod = false;
    
#ifdef MORTARONELOOP
    // this is the special case
    // integration of M, LinM and D, LinD is done here in combination
    // (evaluation of occuring terms is handled similarly)
    dod = true;
#endif // #ifdef MORTARONELOOP

    // compute cell D/M matrix *******************************************
    if (shapefcn_ == MortarInterface::StandardFunctions)
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
    
    else if (shapefcn_ == MortarInterface::DualFunctions)
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
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave             |
 |  element to the D map of the adjacent slave nodes.                   |
 *----------------------------------------------------------------------*/
bool MORTAR::Integrator::AssembleD(const Epetra_Comm& comm,
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
            snode->AddMValue(sdof,col,minusval);
          }
          else
          {	
            snode->AddDValue(sdof,col,val);
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
 |  Assemble M contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave and master  |
 |  overlap pair to the M map of the adjacent slave nodes.              |
 *----------------------------------------------------------------------*/
bool MORTAR::Integrator::AssembleM(const Epetra_Comm& comm,
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
bool MORTAR::Integrator::AssembleMmod(const Epetra_Comm& comm,
                                       MORTAR::MortarElement& sele,
                                       MORTAR::MortarElement& mele,
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
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
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
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
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

#endif //#ifdef CCADISCRET
