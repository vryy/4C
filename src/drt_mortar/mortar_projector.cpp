/*!----------------------------------------------------------------------
\file mortar_projector.cpp
\brief A class to perform projections of nodes onto opposing MortarElements

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

#include "mortar_projector.H"
#include "mortar_interface.H"
#include "mortar_element.H"
#include "mortar_node.H"
#include "mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "mortar_calc_utils.H"
/*----------------------------------------------------------------------*
 |  impl...                                                  farah 01/14|
 *----------------------------------------------------------------------*/
MORTAR::MortarProjector* MORTAR::MortarProjector::Impl(MortarElement& ele)
{
  switch (ele.Shape())
  {
  case DRT::Element::quad4:
  {
    return MortarProjectorCalc<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return MortarProjectorCalc<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return MortarProjectorCalc<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return MortarProjectorCalc<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return MortarProjectorCalc<DRT::Element::tri6>::Instance();
  }
  case DRT::Element::line2:
  {
    return MortarProjectorCalc<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return MortarProjectorCalc<DRT::Element::line3>::Instance();
  }
  //==================================================
  //                     NURBS
  //==================================================
  case DRT::Element::nurbs2:
  {
    return MortarProjectorCalc<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:
  {
    return MortarProjectorCalc<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:
  {
    return MortarProjectorCalc<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs8:
  {
    return MortarProjectorCalc<DRT::Element::nurbs8>::Instance();
  }
  case DRT::Element::nurbs9:
  {
    return MortarProjectorCalc<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele.Shape(), ele.NumNode()); break;
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
MORTAR::MortarProjectorCalc<distype>::MortarProjectorCalc()
{
  //nothing
}

/*----------------------------------------------------------------------*
 |  Instance (public)                                        farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
MORTAR::MortarProjectorCalc<distype> * MORTAR::MortarProjectorCalc<distype>::Instance(bool create)
{
  static MortarProjectorCalc<distype> * instance;
  if (create)
  {
    if (instance==NULL)
      instance = new MortarProjectorCalc<distype>();
  }
  else
  {
    if (instance!=NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 |  Done (public)                                             farah 01/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void MORTAR::MortarProjectorCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}

/*----------------------------------------------------------------------*
 |  Project a node along its nodal normal (public)            popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectNodalNormal(MORTAR::MortarNode& node,
                                                 MORTAR::MortarElement& ele,
                                                 double* xi)
{
  bool ok = true;
  if (ndim_==2)
  {
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double f = 0.0;
    double df = 0.0;
    int k=0;

    for (k=0;k<MORTARMAXITER;++k)
    {
      f=EvaluateFNodalNormal(node,ele,eta);
      if (abs(f) < MORTARCONVTOL) break;
      df=EvaluateGradFNodalNormal(node,ele,eta);
      if (abs(df)<1.0e-12) dserror("ERROR: Singular Jacobian for projection");
      eta[0]+=(-f)/df;
    }

    // get the result
    xi[0]=eta[0];

    // Newton iteration unconverged
    if (abs(f) > MORTARCONVTOL)
    {
      ok = false;
      xi[0] = 1.0e12;

      // Here (S->M projection) we only give a warning, no error!!!
      // This iteration sometimes diverges, when the projection is far off.
      // These cases are harmless, as these nodes then do not participate in
      // the overlap detection anyway!
      //std::cout << "***WARNING*** ProjectNodalNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and MortarElementID " << ele.Id() << std::endl;
    }

    /*
    // Newton iteration converged
    else
    {
      std::cout << "Newton iteration converged in " << k << " step(s)!" << std::endl
           << "The result is: " << xi[0] << std::endl;
    }*/
  }

  else dserror("ERROR: ProjectNodalNormal: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a node along element's normal field (public)      popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectElementNormal(MORTAR::MortarNode& node,
                                                   MORTAR::MortarElement& ele,
                                                   double* xi)
{
  bool ok = true;
  if (ndim_==2)
  {
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double f = 0.0;
    double df = 0.0;
    int k=0;

    for (k=0;k<MORTARMAXITER;++k)
    {
      f=EvaluateFElementNormal(node,ele,eta);
      if (abs(f) < MORTARCONVTOL) break;
      df=EvaluateGradFElementNormal(node,ele,eta);
      if (abs(df)<1.0e-12) dserror("ERROR: Singular Jacobian for projection");
      eta[0]+=(-f)/df;
    }

    // get the result
    xi[0]=eta[0];

    // Newton iteration unconverged
    if (abs(f) > MORTARCONVTOL)
    {
      ok = false;
      xi[0] = 1.0e12;

      // Here (M->S projection) we only give a warning, no error!!!
      // This iteration sometimes diverges, when the projection is far off.
      // These cases are harmless, as these nodes then do not participate in
      // the overlap detection anyway!
      //std::cout << "***WARNING*** ProjectElementNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and MortarElementID " << ele.Id() << std::endl;
    }

    /*
    // Newton iteration converged
    else
    {
      std::cout << "Newton iteration converged in " << k << " step(s)!" << std::endl
           << "The result is: " << xi[0] << std::endl;
    }*/
  }

  else dserror("ERROR: ProjectElementNormal: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (public)           popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectGaussPoint(MORTAR::MortarElement& gpele,
                                                const double* gpeta,
                                                MORTAR::MortarElement& ele,
                                                double* xi)
{
  bool ok = true;
  if (ndim_==2)
  {
    LINALG::Matrix<n_,1>          val;
    LINALG::Matrix<ndim_,n_>   coord;

    DRT::Node** mynodes = gpele.Nodes();
    if(!mynodes) dserror("ERROR: ProjectGaussPoint: Null pointer!");

    // get shape function values and derivatives at gpeta
    if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
    {
      LINALG::SerialDenseVector auxval(n_);
      LINALG::SerialDenseMatrix deriv(n_,1);
      gpele.EvaluateShape(gpeta, auxval, deriv, gpele.NumNode());

      for(int i=0;i<n_;++i)
        val(i)=auxval(i);
    }
    else
      DRT::UTILS::shape_function_1D (val ,gpeta[0],distype);

    // get interpolated GP normal and GP coordinates
    double gpn[ndim_];
    double gpx[ndim_];
    for (int i=0;i<ndim_;++i)
    {
      gpn[i] = 0.0;
      gpx[i] = 0.0;
    }

    for (int i=0;i<n_;++i)
    {
      MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);

      for (int j=0;j<ndim_;++j)
      {
        gpn[j]     += val(i)*mymrtrnode->MoData().n()[j];

        coord(j,i) =  mymrtrnode->xspatial()[j];

        gpx[j]     +=val(i)*coord(j,i);
      }
    }

    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double f      = 0.0;
    double df     = 0.0;
    int k         = 0;

    for (k=0;k<MORTARMAXITER;++k)
    {
      f=EvaluateFGaussPoint(gpx,gpn,ele,eta);
      if (abs(f) < MORTARCONVTOL) break;
      df=EvaluateGradFGaussPoint(gpn,ele,eta);
      if (abs(df)<1.0e-12) dserror("ERROR: Singular Jacobian for projection");
      eta[0]+=(-f)/df;
    }

    // get the result
    xi[0]=eta[0];

    // Newton iteration unconverged
    if (abs(f) > MORTARCONVTOL)
    {
      ok = false;
      xi[0] = 1.0e12;

      dserror("ERROR: ProjectGaussPoint: Newton unconverged for GP at xi=%d"
              " from MortarElementID %i", gpeta[0], gpele.Id());
    }
  }

  else dserror("ERROR: ProjectGaussPoint: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 | Check projection for warped elements quad4 elements       farah 01/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::CheckProjection4AUXPLANE(MORTAR::MortarElement& ele,
                                             double* ngp, double* globgp )
{
  if (ele.Shape()==DRT::Element::tri3)
    dserror("ELEMENT SHAPE TRI3 -- NO WARPING");

  if (ele.Shape()!=DRT::Element::quad4)
  {
    return true;
  }
  
  int nnode = ele.NumNode();
  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("ERROR: Project: Null pointer!");

  //compute base-vectors
  std::vector<double> t0(3);
  std::vector<double> t1(3);
  std::vector<double> auxn(3);
  std::vector<double> auxc(3);
  std::vector<double> proj_gp(3);
  LINALG::Matrix<3,3> P;
  LINALG::Matrix<3,3> T;
  double n[3]={0.0 , 0.0 , 0.0};
  double length_t=0.0;
  double length_n=0.0;
  double a1=0.0;
  bool all_negative=true;
  MortarNode* mycnode_0=0;
  MortarNode* mycnode_1=0;
  MortarNode* mycnode_2=0;

  //project gp onto auxn
  for (int i=0;i<nnode;++i) //loop over edges
  {
    if(i==0)
    {
      mycnode_1 = static_cast<MortarNode*> (mynodes[0]);
      if (!mycnode_1) dserror("ERROR: Project: Null pointer!");

      mycnode_0 = static_cast<MortarNode*> (mynodes[3]);
      if (!mycnode_0) dserror("ERROR: Project: Null pointer!");

      mycnode_2 = static_cast<MortarNode*> (mynodes[1]);
      if (!mycnode_2) dserror("ERROR: Project: Null pointer!");
    }
    if(i==3)
    {
      mycnode_1 = static_cast<MortarNode*> (mynodes[3]);
      if (!mycnode_1) dserror("ERROR: Project: Null pointer!");

      mycnode_0 = static_cast<MortarNode*> (mynodes[2]);
      if (!mycnode_0) dserror("ERROR: Project: Null pointer!");

      mycnode_2 = static_cast<MortarNode*> (mynodes[0]);
      if (!mycnode_2) dserror("ERROR: Project: Null pointer!");
    }
    if (i==1 || i==2)
    {
      mycnode_1 = static_cast<MortarNode*> (mynodes[i]);
      if (!mycnode_1) dserror("ERROR: Project: Null pointer!");

      mycnode_0 = static_cast<MortarNode*> (mynodes[i-1]);
      if (!mycnode_0) dserror("ERROR: Project: Null pointer!");

      mycnode_2 = static_cast<MortarNode*> (mynodes[i+1]);
      if (!mycnode_2) dserror("ERROR: Project: Null pointer!");
    }

    //span vectors -- edges
    t0[0]=mycnode_0->xspatial()[0] - mycnode_1->xspatial()[0];
    t0[1]=mycnode_0->xspatial()[1] - mycnode_1->xspatial()[1];
    t0[2]=mycnode_0->xspatial()[2] - mycnode_1->xspatial()[2];

    t1[0]=mycnode_2->xspatial()[0] - mycnode_1->xspatial()[0];
    t1[1]=mycnode_2->xspatial()[1] - mycnode_1->xspatial()[1];
    t1[2]=mycnode_2->xspatial()[2] - mycnode_1->xspatial()[2];

    auxc[0]=mycnode_1->xspatial()[0];
    auxc[1]=mycnode_1->xspatial()[1];
    auxc[2]=mycnode_1->xspatial()[2];

    auxn[0]=t1[1]*t0[2] - t1[2]*t0[1];
    auxn[1]=t1[2]*t0[0] - t1[0]*t0[2];
    auxn[2]=t1[0]*t0[1] - t1[1]*t0[0];

    //fill Projection matrix P
    for (int j=0;j<3;++j)
      P(j,2)=-ngp[j];

    for (int j=0;j<3;++j)
      P(j,0)=t0[j];

    for (int j=0;j<3;++j)
      P(j,1)=t1[j];

    P.Invert();
    double lambda1= P(0,0)*(globgp[0]-auxc[0]) + P(0,1)*(globgp[1]-auxc[1]) + P(0,2)*(globgp[2]-auxc[2]);
    double lambda2= P(1,0)*(globgp[0]-auxc[0]) + P(1,1)*(globgp[1]-auxc[1]) + P(1,2)*(globgp[2]-auxc[2]);
    //double alph= P(2,0)*(globgp[0]-auxc[0]) + P(2,1)*(globgp[1]-auxc[1]) + P(2,2)*(globgp[2]-auxc[2]);

    proj_gp[0]=lambda1*t0[0] + lambda2*t1[0]+auxc[0];
    proj_gp[1]=lambda1*t0[1] + lambda2*t1[1]+auxc[1];
    proj_gp[2]=lambda1*t0[2] + lambda2*t1[2]+auxc[2];
    //std::cout << "proj_gp_AUX4PLANE= " << proj_gp[0] << " "  << proj_gp[1] << " " << proj_gp[2] << std::endl;

    //check
    //edge 1
    n[0]=-(t0[1]*auxn[2] - t0[2]*auxn[1]);
    n[1]=-(t0[2]*auxn[0] - t0[0]*auxn[2]);
    n[2]=-(t0[0]*auxn[1] - t0[1]*auxn[0]);

    length_t=sqrt(t0[0]*t0[0] + t0[1]*t0[1] +t0[2]*t0[2]);
    length_n=sqrt(n[0]*n[0] + n[1]*n[1] +n[2]*n[2]);
    //fill Projection matrix T
    for (int j=0;j<3;++j)
      T(j,0)=n[j]/length_n;

    for (int j=0;j<3;++j)
      T(j,1)=t0[j]/length_t;

    for (int j=0;j<3;++j)
      T(j,2)=auxc[j];

    T.Invert();
    a1=T(0,0)*proj_gp[0] + T(0,1)*proj_gp[1] + T(0,2)*proj_gp[2];

    if (a1>0.0)
      all_negative=false;

    //edge 2
    n[0]=(t1[1]*auxn[2] - t1[2]*auxn[1]);
    n[1]=(t1[2]*auxn[0] - t1[0]*auxn[2]);
    n[2]=(t1[0]*auxn[1] - t1[1]*auxn[0]);

    length_t=sqrt(t1[0]*t1[0] + t1[1]*t1[1] +t1[2]*t1[2]);
    length_n=sqrt(n[0]*n[0] + n[1]*n[1] +n[2]*n[2]);
    //fill Projection matrix T
    for (int j=0;j<3;++j)
      T(j,0)=n[j]/length_n;

    for (int j=0;j<3;++j)
      T(j,1)=t1[j]/length_t;

    for (int j=0;j<3;++j)
      T(j,2)=auxc[j];

    T.Invert();

    a1=T(0,0)*proj_gp[0] + T(0,1)*proj_gp[1] + T(0,2)*proj_gp[2];
    //a2=T(1,0)*proj_gp[0] + T(1,1)*proj_gp[1] + T(1,2)*proj_gp[2];

    if (a1>0.0)
      all_negative=false;
  }

  if(all_negative==false)
    return true;

  return false; // creates error in ProjectGaussPoint3D()
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (3D)               popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectGaussPoint3D(MORTAR::MortarElement& gpele,
                                                  const double* gpeta,
                                                  MORTAR::MortarElement& ele,
                                                  double* xi, double& par)
{
  if (ndim_==3)
  {
    LINALG::Matrix<n_,1>          val;
    LINALG::Matrix<ndim_,n_>      coord;
    coord.Clear();

    DRT::Node** mynodes = gpele.Nodes();
    if(!mynodes) dserror("ERROR: ProjectGaussPoint: Null pointer!");

    // get shape function values and derivatives at gpeta
    if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs8 ||
       distype==DRT::Element::nurbs9)
    {
      LINALG::SerialDenseVector auxval(n_);
      LINALG::SerialDenseMatrix deriv(n_,1);
      gpele.EvaluateShape(gpeta, auxval, deriv, gpele.NumNode());

      for(int i=0;i<n_;++i)
        val(i)=auxval(i);
    }
    else
      DRT::UTILS::shape_function_2D (val,gpeta[0],gpeta[1],distype);

    // get interpolated GP normal and GP coordinates
    double gpn[3];
    double gpx[3];
    for (int i=0;i<ndim_;++i)
    {
      gpn[i] = 0.0;
      gpx[i] = 0.0;
    }

    for (int i=0;i<n_;++i)
    {
      MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);

      for (int j=0;j<ndim_;++j)
      {
        gpn[j]     += val(i)*mymrtrnode->MoData().n()[j];

        coord(j,i) =  mymrtrnode->xspatial()[j];

        gpx[j]     +=val(i)*coord(j,i);
      }
    }

    // start in the element center
    DRT::Element::DiscretizationType dt = ele.Shape();
    double eta[2] = {0.0, 0.0};
    if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
    {
      eta[0] = 1.0/3;
      eta[1] = 1.0/3;
    }

    // auxiliary variable
    double alpha = 0.0;

    // function f (vector-valued)
    double f[3] = {0.0, 0.0, 0.0};

    // gradient of f (df/deta[0], df/deta[1], df/dalpha)
    LINALG::Matrix<3,3> df;

    // start iteration
    int k=0;
    double conv = 0.0;

    for (k=0;k<MORTARMAXITER;++k)
    {
      EvaluateFGaussPoint3D(f,gpx,gpn,ele,eta,alpha);
      conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
      if (conv <= MORTARCONVTOL) break;
      EvaluateGradFGaussPoint3D(df,gpx,gpn,ele,eta,alpha);

      // safety check: if projection normal is parallel to the master element --> det can be zero
      double det = df.Determinant();
      if (det > -1e-12 and det < 1e-12)
      {
        std::cout << "WARNING: GPProjection parallel to master element --> GP skipped for this master element!" << std::endl;
        // leave here
        return false;
      }

      // solve deta = - inv(df) * f
      double jacdet = df.Invert();
      if (abs(jacdet)<1.0e-12) dserror("ERROR: Singular Jacobian for projection");

      // update eta and alpha
      eta[0] += -df(0,0)*f[0] - df(0,1)*f[1] - df(0,2)*f[2];
      eta[1] += -df(1,0)*f[0] - df(1,1)*f[1] - df(1,2)*f[2];
      alpha  += -df(2,0)*f[0] - df(2,1)*f[1] - df(2,2)*f[2];
      
      //Projection Check
      if (k==MORTARMAXITER-1) 
      {
        bool check = CheckProjection4AUXPLANE(ele, gpn,gpx); 
        if (check==false)
          dserror("!!! STOP !!!   -->   Projection Error: Newton unconverged but GP on mele !!!");
      }
    }

    // Newton iteration unconverged
    if (conv > MORTARCONVTOL)
    {
      xi[0]=1e12;
      xi[1]=1e12;
    }
    else
    {
      xi[0]=eta[0];
      xi[1]=eta[1];
    }
    
    par=alpha;
  }

  else dserror("ERROR: ProjectGaussPoint: Called 3D version for 2D problem!");

  return true;
}
/*----------------------------------------------------------------------*
 |  Project a Gauss point along AuxPlane normal (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectGaussPointAuxn3D(const double* globgp,
                                                      const double* auxn,
                                                      MORTAR::MortarElement& ele,
                                                      double* xi, double& par)
{
  if (ndim_==3)
  {
    // start in the element center
    DRT::Element::DiscretizationType dt = ele.Shape();
    double eta[2] = {0.0, 0.0};
    if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
    {
      eta[0] = 1.0/3;
      eta[1] = 1.0/3;
    }

    // auxiliary variable
    double alpha = 0.0;

    // function f (vector-valued)
    double f[3] = {0.0, 0.0, 0.0};

    // gradient of f (df/deta[0], df/deta[1], df/dalpha)
    LINALG::Matrix<3,3> df;

    // start iteration
    int k=0;
    double conv = 0.0;

    for (k=0;k<MORTARMAXITER;++k)
    {
      EvaluateFGaussPointAuxn3D(f,globgp,auxn,ele,eta,alpha);
      conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
      //std::cout << "Iteration " << k << ": -> |f|=" << conv << std::endl;
      if (conv <= MORTARCONVTOL) break;
      EvaluateGradFGaussPointAuxn3D(df,globgp,auxn,ele,eta,alpha);

      // solve deta = - inv(df) * f
      double jacdet = df.Invert();
      if (abs(jacdet)<1.0e-12) dserror("ERROR: Singular Jacobian for projection");

      // update eta and alpha
      eta[0] += -df(0,0)*f[0] - df(0,1)*f[1] - df(0,2)*f[2];
      eta[1] += -df(1,0)*f[0] - df(1,1)*f[1] - df(1,2)*f[2];
      alpha  += -df(2,0)*f[0] - df(2,1)*f[1] - df(2,2)*f[2];
    }

    // Newton iteration unconverged
    if (conv > MORTARCONVTOL)
      dserror("ERROR: ProjectGaussPointAuxn3D: Newton unconverged for GP"
              "at xi = (%f,%f,%f) onto MortarElementID %i", globgp[0],globgp[1],globgp[2],ele.Id());

    // Newton iteration converged
    xi[0]=eta[0];
    xi[1]=eta[1];
    par = alpha;
    //std::cout << "Newton iteration converged in " << k << " steps!" << std::endl;
  }

  else dserror("ERROR: ProjectGaussPoint: Called 3D version for 2D problem!");

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for nodal normal case (public)                 popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateFNodalNormal(MORTAR::MortarNode& node,
                                                     MORTAR::MortarElement& ele,
                                                     const double* eta)
{
  /* Evaluate the function F(eta) = ( Ni * xim - xs ) x ns,
     or to be more precise the third component of this vector function!

       Ni  shape functions of element
       xim coords of element nodes (master side)
       xs  coords of node to be projected (slave side)
       ns   outward normal of node to be projected (slave side)          */
  double fval = 0.0;

  // build interpolation of master node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nx,0);

  // subtract slave node coordinates
  for (int i=0;i<ndim_;++i)
    nx[i]-=node.xspatial()[i];

  //calculate F
  fval = nx[0]*node.MoData().n()[1]-nx[1]*node.MoData().n()[0];

  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for nodal normal case (public)             popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateGradFNodalNormal(MORTAR::MortarNode& node,
                                                         MORTAR::MortarElement& ele,
                                                         const double* eta)
{
  /* Evaluate the function GradF(eta)
     = Ni,eta * xim * nys - Ni,eta * yim * nxs,

       Ni,eta    shape function derivatives of element
       xim, yim  coords of element nodes (master side)
       nxs, nys   outward normal of node to be projected (slave side)   */

  double fgrad = 0.0;

  // build interpolation of master node coordinates for current eta
  // use shape function derivatives for interpolation (hence "1")
  double nxeta[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nxeta,1);

  // calculate GradF
  fgrad =  nxeta[0]*node.MoData().n()[1]-nxeta[1]*node.MoData().n()[0];

  return fgrad;

}

/*----------------------------------------------------------------------*
 |  Evaluate F for element normal case (public)               popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateFElementNormal(MORTAR::MortarNode& node,
                                                       MORTAR::MortarElement& ele,
                                                       const double* eta)
{
  /* Evaluate the function F(eta) = ( Ni * xis - xm ) x ( Nj * njs),
     or to be more precise the third component of this vector function!

       Ni  shape functions of element
       xis coords of element nodes (slave side)
       xm  coords of node to be projected (master side)
       nis outward normals of element nodes (slave side)                */

  double fval = 0.0;

  // collect necessary data (slave side)
  DRT::Node** mynodes = ele.Nodes();
  if(!mynodes) dserror("ERROR: EvaluateFElementNormal: Null pointer!");

  LINALG::Matrix<n_,1>          val;
  LINALG::Matrix<ndim_,n_>   coord;

  // get shape function values and derivatives at gpeta
  if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
  {
    LINALG::SerialDenseVector auxval(n_);
    LINALG::SerialDenseMatrix deriv(n_,1);
    ele.EvaluateShape(eta, auxval, deriv, ele.NumNode());

    for(int i=0;i<n_;++i)
      val(i)=auxval(i);
  }
  else
    DRT::UTILS::shape_function_1D (val,eta[0],distype);

  // get interpolated normal and proj. coordinates for current eta
  double nn[ndim_];
  double nx[ndim_];
  for (int j=0;j<ndim_;++j)
  {
    nn[j]=0.0;
    nx[j]=0.0;
  }

  for (int i=0;i<n_;++i)
  {
    MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);

    for (int j=0;j<ndim_;++j)
    {
      nn[j]+=val(i)*mymrtrnode->MoData().n()[j];

      coord(j,i) = mymrtrnode->xspatial()[j];

      nx[j]+=val(i)*coord(j,i);
    }
  }

  // subtract master node coordinates
  for (int j=0;j<ndim_;++j)
    nx[j]-=node.xspatial()[j];

  // calculate F
  fval = nx[0]*nn[1]-nx[1]*nn[0];

  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for element normal case (public)           popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateGradFElementNormal(MORTAR::MortarNode& node,
                                                           MORTAR::MortarElement& ele,
                                                           const double* eta)
{
  if (ndim_==3)
    dserror("This Projector is only for 2D Problems!");

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

  // collect necessary data (slave side)
  LINALG::Matrix<n_,1>        val;
  LINALG::Matrix<ndim_-1,n_>  deriv;
  LINALG::Matrix<ndim_,n_>    coord;

  DRT::Node** mynodes = ele.Nodes();
  if(!mynodes) dserror("ERROR: EvaluateGradFElementNormal: Null pointer!");

  // get shape function values and derivatives at gpeta
  if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
  {
    LINALG::SerialDenseVector auxval(n_);
    LINALG::SerialDenseMatrix auxderiv(n_,1);
    ele.EvaluateShape(eta, auxval, auxderiv, ele.NumNode());

    for(int i=0;i<n_;++i)
    {
      val(i)=auxval(i);
      deriv(0,i)=auxderiv(i,0);
    }
  }
  else
  {
    DRT::UTILS::shape_function_1D        (val,eta[0],distype);
    DRT::UTILS::shape_function_1D_deriv1 (deriv,eta[0],distype);
  }

  // get interpolated normal and proj. coordinates for current eta
  double nn[ndim_];
  double nneta[ndim_];
  double nx[ndim_];
  double nxeta[ndim_];
  for (int j=0;j<ndim_;++j)
  {
    nn[j]=0.0;
    nneta[j]=0.0;
    nx[j]=0.0;
    nxeta[j]=0.0;
  }

  for (int i=0;i<n_;++i)
  {
    MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);

    for (int j=0;j<ndim_;++j)
    {
      nn[j]     += val(i)*mymrtrnode->MoData().n()[j];
      nneta[j]  += deriv(0,i)*mymrtrnode->MoData().n()[j];
      coord(j,i) = mymrtrnode->xspatial()[j];
      nx[j]     += val(i)*coord(j,i);
      nxeta[j]  += deriv(0,i)*coord(j,i);
    }
  }

  // subtract master node coordinates
  for(int j=0;j<ndim_;++j)
    nx[j]-=node.xspatial()[j];

  // calculate GradF
  fgrad =   nxeta[0]*nn[1] + nx[0]*nneta[1]
          - nxeta[1]*nn[0] - nx[1]*nneta[0];

  return fgrad;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (public)                  popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateFGaussPoint(const double* gpx,
                                                    const double* gpn,
                                                    MORTAR::MortarElement& ele,
                                                    const double* eta)
{
  /* Evaluate the function F(eta) = ( Ni * xim - gpx ) x gpn,
     or to be more precise the third component of this vector function!

       Ni  shape functions of element (master side)
       xim coords of element nodes (master side)
       gpx coords of GP to be projected (slave side)
       gpn outward normal of GP to be projected (slave side)          */

  double fval = 0.0;

  // build interpolation of master node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nx,0);

  // subtract GP coordinates
  nx[0]-=gpx[0];
  nx[1]-=gpx[1];

  // calculate F
  fval = nx[0]*gpn[1]-nx[1]*gpn[0];

  return fval;

}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (public)              popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateGradFGaussPoint(const double* gpn,
                                                        MORTAR::MortarElement& ele,
                                                        const double* eta)
{
  /* Evaluate the function GradF(eta)
      = Ni,eta * xim * gpny - Ni,eta * yim * gpnx,

       Ni,eta     shape function derivatives of element (master side)
       xim, yim   coords of element nodes (master side)
       gpnx, gpny outward normal of GP to be projected (slave side)   */

  double fgrad = 0.0;

  // build interpolation of master node coordinates for current eta
  // use shape function derivatives for interpolation (hence "1")
  double nxeta[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nxeta,1);

  // calculate GradF
  fgrad = nxeta[0]*gpn[1]-nxeta[1]*gpn[0];

  return fgrad;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (3D)                      popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::EvaluateFGaussPoint3D(
                              double* f,
                              const double* gpx,
                              const double* gpn,
                              MORTAR::MortarElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi - alpha * gpn - gpx
     which is a vector-valued function with 3 components!

       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       gpx     coords of GP to be projected
       gpn     normal of GP along which to project                  */

  // build interpolation of ele node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nx,0);

  // evaluate function f
  for (int i=0; i<ndim_; ++i)
    f[i] = nx[i] - alpha * gpn[i] - gpx[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (3D)                  popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::EvaluateGradFGaussPoint3D(
                              LINALG::Matrix<3,3>& fgrad,
                              const double* gpx,
                              const double* gpn,
                              MORTAR::MortarElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
     - alpha * gpn - gpx, which is a (3x3)-matrix!

       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       gpx     coords of GP to be projected
       gpn     normal of GP along which to project                  */

  // build interpolation of ele node coordinates for current eta
//  double nxeta1[3] = {0.0, 0.0, 0.0};
//  double nxeta2[3] = {0.0, 0.0, 0.0};
//  ele.LocalToGlobal(eta,nxeta1,1);
//  ele.LocalToGlobal(eta,nxeta2,2);

  double nxeta1[ndim_];
  double nxeta2[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nxeta1,1);
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nxeta2,2);

  //evaluate function f gradient
  for (int i=0;i<ndim_;++i)
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
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::EvaluateFGaussPointAuxn3D(
                              double* f,
                              const double* globgp,
                              const double* auxn,
                              MORTAR::MortarElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi - alpha * auxn - globgp
     which is a vector-valued function with 3 components!

       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       globgp  coords of AuxPlaneGP to be projected
       auxn    normal of AuxPlane along which to project            */

  // build interpolation of ele node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nx,0);

  // evaluate function f
  for (int i=0; i<ndim_; ++i)
    f[i] = nx[i] - alpha * auxn[i] - globgp[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for AuxPlane Gauss point case (3D)         popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::EvaluateGradFGaussPointAuxn3D(
                              LINALG::Matrix<3,3>& fgrad,
                              const double* globgp,
                              const double* auxn,
                              MORTAR::MortarElement& ele,
                              const double* eta,
                              const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
     - alpha * auxn - globgp, which is a (3x3)-matrix!

       Ni      shape functions of element to project on
       xi      coords of nodes of element to project on
       globgp  coords of AuxPlaneGP to be projected
       auxn    normal of AuxPlane along which to project            */

  // build interpolation of ele node coordinates for current eta
  double nxeta1[ndim_];
  double nxeta2[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nxeta1,1);
  MORTAR::UTILS::LocalToGlobal<distype>(ele,eta,nxeta2,2);

  // evaluate function f gradient
  for (int i=0;i<ndim_;++i)
  {
    fgrad(i,0) = nxeta1[i];
    fgrad(i,1) = nxeta2[i];
    fgrad(i,2) = -auxn[i];
  }

  return true;
}

