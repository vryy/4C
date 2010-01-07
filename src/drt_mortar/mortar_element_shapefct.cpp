/*!----------------------------------------------------------------------
\file mortar_element_shapefct.cpp
\brief Shape function repository for mortar coupling element

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

#include "mortar_element.H"
#include "mortar_integrator.H"
#include "mortar_defines.H"
#include "../drt_lib/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  1D/2D shape function repository                           popp 04/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::ShapeFunctions(MortarElement::ShapeType shape,
                                       const double* xi,
                                       LINALG::SerialDenseVector& val,
                                       LINALG::SerialDenseMatrix& deriv)
{
  switch(shape)
  {
  // *********************************************************************
  // 1D standard linear shape functions (line2)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::lin1D:
  {
    val[0] = 0.5*(1-xi[0]);
    val[1] = 0.5*(1+xi[0]);
    deriv(0,0) = -0.5;
    deriv(1,0) =  0.5;
    break;
  }
  // *********************************************************************
  // 2D standard linear shape functions (tri3)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::lin2D:
  {
    val[0] = 1-xi[0]-xi[1];
    val[1] = xi[0];
    val[2] = xi[1];
    deriv(0,0) = -1.0; deriv(0,1) = -1.0;
    deriv(1,0) =  1.0; deriv(1,1) =  0.0;
    deriv(2,0) =  0.0; deriv(2,1) =  1.0;
    break;
  }
  // *********************************************************************
  // 2D standard blinear shape functions (quad4)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::bilin2D:
  {
    val[0] = 0.25*(1-xi[0])*(1-xi[1]);
    val[1] = 0.25*(1+xi[0])*(1-xi[1]);
    val[2] = 0.25*(1+xi[0])*(1+xi[1]);
    val[3] = 0.25*(1-xi[0])*(1+xi[1]);
    deriv(0,0) = -0.25*(1-xi[1]); deriv(0,1) = -0.25*(1-xi[0]);
    deriv(1,0) =  0.25*(1-xi[1]); deriv(1,1) = -0.25*(1+xi[0]);
    deriv(2,0) =  0.25*(1+xi[1]); deriv(2,1) =  0.25*(1+xi[0]);
    deriv(3,0) = -0.25*(1+xi[1]); deriv(3,1) =  0.25*(1-xi[0]);
    break;
  }
  // *********************************************************************
  // 1D standard quadratic shape functions (line3)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::quad1D:
  {
    val[0] = 0.5*xi[0]*(xi[0]-1);
    val[1] = 0.5*xi[0]*(xi[0]+1);
    val[2] = (1-xi[0])*(1+xi[0]);
    deriv(0,0) = xi[0]-0.5;
    deriv(1,0) = xi[0]+0.5;
    deriv(2,0) = -2*xi[0];
    break;
  }
  // *********************************************************************
  // 2D standard quadratic shape functions (tri6)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::quad2D:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double t1 = 1.0-r-s;
    const double t2 = r;
    const double t3 = s;

    val[0] = t1*(2.0*t1-1.0);
    val[1] = t2*(2.0*t2-1.0);
    val[2] = t3*(2.0*t3-1.0);
    val[3] = 4.0*t2*t1;
    val[4] = 4.0*t2*t3;
    val[5] = 4.0*t3*t1;

    deriv(0,0)= -3.0+4.0*(r+s);
    deriv(0,1)= -3.0+4.0*(r+s);
    deriv(1,0)= 4.0*r-1.0;
    deriv(1,1)= 0.0;
    deriv(2,0)= 0.0;
    deriv(2,1)= 4.0*s-1.0;
    deriv(3,0)= 4.0*(1.0-2.0*r-s);
    deriv(3,1)=-4.0*r;
    deriv(4,0)= 4.0*s;
    deriv(4,1)= 4.0*r;
    deriv(5,0)=-4.0*s;
    deriv(5,1)= 4.0*(1.0-r-2.0*s);

    break;
  }
  // *********************************************************************
  // 2D serendipity shape functions (quad8)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::serendipity2D:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double rp=1.0+r;
    const double rm=1.0-r;
    const double sp=1.0+s;
    const double sm=1.0-s;
    const double r2=1.0-r*r;
    const double s2=1.0-s*s;

    // values for centernodes are straight forward
    //      0.5*(1-xi*xi)*(1-eta) (0 for xi=+/-1 and eta=+/-1/0
    //                             0 for xi=0    and eta= 1
    //                             1 for xi=0    and eta=-1    )
    // use shape functions on centernodes to zero out the corner node
    // shape functions on the centernodes
    // (0.5 is the value of the linear shape function in the centernode)
    //
    //  0.25*(1-xi)*(1-eta)-0.5*funct[neighbor1]-0.5*funct[neighbor2]

    val[0]=0.25*(rm*sm-(r2*sm+s2*rm));
    val[1]=0.25*(rp*sm-(r2*sm+s2*rp));
    val[2]=0.25*(rp*sp-(s2*rp+r2*sp));
    val[3]=0.25*(rm*sp-(r2*sp+s2*rm));
    val[4]=0.5*r2*sm;
    val[5]=0.5*s2*rp;
    val[6]=0.5*r2*sp;
    val[7]=0.5*s2*rm;

    deriv(0,0)= 0.25*sm*(2*r+s);
    deriv(0,1)= 0.25*rm*(r+2*s);
    deriv(1,0)= 0.25*sm*(2*r-s);
    deriv(1,1)= 0.25*rp*(2*s-r);
    deriv(2,0)= 0.25*sp*(2*r+s);
    deriv(2,1)= 0.25*rp*(r+2*s);
    deriv(3,0)= 0.25*sp*(2*r-s);
    deriv(3,1)= 0.25*rm*(2*s-r);
    deriv(4,0)=-sm*r;
    deriv(4,1)=-0.5*rm*rp;
    deriv(5,0)= 0.5*sm*sp;
    deriv(5,1)=-rp*s;
    deriv(6,0)=-sp*r;
    deriv(6,1)= 0.5*rm*rp;
    deriv(7,0)=-0.5*sm*sp;
    deriv(7,1)=-rm*s;

    break;
  }
  // *********************************************************************
  // 2D standard biquadratic shape functions (quad9)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case MortarElement::biquad2D:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double rp=1.0+r;
    const double rm=1.0-r;
    const double sp=1.0+s;
    const double sm=1.0-s;
    const double r2=1.0-r*r;
    const double s2=1.0-s*s;
    const double rh=0.5*r;
    const double sh=0.5*s;
    const double rs=rh*sh;
    const double rhp=r+0.5;
    const double rhm=r-0.5;
    const double shp=s+0.5;
    const double shm=s-0.5;

    val[0]= rs*rm*sm;
    val[1]=-rs*rp*sm;
    val[2]= rs*rp*sp;
    val[3]=-rs*rm*sp;
    val[4]=-sh*sm*r2;
    val[5]= rh*rp*s2;
    val[6]= sh*sp*r2;
    val[7]=-rh*rm*s2;
    val[8]= r2*s2;

    deriv(0,0)=-rhm*sh*sm;
    deriv(0,1)=-shm*rh*rm;
    deriv(1,0)=-rhp*sh*sm;
    deriv(1,1)= shm*rh*rp;
    deriv(2,0)= rhp*sh*sp;
    deriv(2,1)= shp*rh*rp;
    deriv(3,0)= rhm*sh*sp;
    deriv(3,1)=-shp*rh*rm;
    deriv(4,0)= 2.0*r*sh*sm;
    deriv(4,1)= shm*r2;
    deriv(5,0)= rhp*s2;
    deriv(5,1)=-2.0*s*rh*rp;
    deriv(6,0)=-2.0*r*sh*sp;
    deriv(6,1)= shp*r2;
    deriv(7,0)= rhm*s2;
    deriv(7,1)= 2.0*s*rh*rm;
    deriv(8,0)=-2.0*r*s2;
    deriv(8,1)=-2.0*s*r2;

    break;
  }
  // *********************************************************************
  // 1D dual linear shape functions (line2)
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case MortarElement::lindual1D:
  {
    val[0] = 0.5*(1-3*xi[0]);
    val[1] = 0.5*(1+3*xi[0]);
    deriv(0,0) = -1.5;
    deriv(1,0) =  1.5;
    break;
  }
  // *********************************************************************
  // 2D dual linear shape functions (tri3)
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case MortarElement::lindual2D:
  {
    val[0] = 3-4*xi[0]-4*xi[1];
    val[1] = 4*xi[0]-1;
    val[2] = 4*xi[1]-1;
    deriv(0,0) = -4.0; deriv(0,1) = -4.0;
    deriv(1,0) =  4.0; deriv(1,1) =  0.0;
    deriv(2,0) =  0.0; deriv(2,1) =  4.0;
    break;
  }
  // *********************************************************************
  // 1D dual quadratic shape functions (line3)
  // 2D dual bilinear shape functions (quad4)
  // 2D dual quadratic shape functions (tri6)
  // 2D dual serendipity shape functions (quad8)
  // 2D dual biquadratic shape functions (quad9)
  // (used for interpolation of Lagrange mutliplier field)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case MortarElement::quaddual1D:
  case MortarElement::bilindual2D:
  case MortarElement::quaddual2D:
  case MortarElement::serendipitydual2D:
  case MortarElement::biquaddual2D:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::Integrator integrator(Shape());

    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
      EvaluateShape(gpc, val, deriv, nnodes);
      detg = Jacobian(gpc);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=integrator.Weight(i)*val[j]*val[k]*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*val[j]*detg;
        }
    }

    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);

    // get solution matrix with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    EvaluateShape(xi, val, deriv, nnodes);

    LINALG::SerialDenseVector valtemp(nnodes,true);
    LINALG::SerialDenseMatrix derivtemp(nnodes,2,true);
    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        valtemp[i]+=ae(i,j)*val[j];
        derivtemp(i,0)+=ae(i,j)*deriv(j,0);
        derivtemp(i,1)+=ae(i,j)*deriv(j,1);
      }

    val=valtemp;
    deriv=derivtemp;

    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (const replacing linear, line2)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case MortarElement::lindual1D_edge0:
  {
    val[0] = 0;
    val[1] = 1;
    deriv(0,0) = 0;
    deriv(1,0) = 0;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (const replacing linear, line2)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case MortarElement::lindual1D_edge1:
  {
    val[0] = 1;
    val[1] = 0;
    deriv(0,0) = 0;
    deriv(1,0) = 0;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear replacing quad, line3)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (only form a basis and have to be adapted for distorted elements)
  // *********************************************************************
  case MortarElement::dual1D_base_for_edge0:
  {
    val[0] = xi[0];
    val[1] = 1-xi[0];
    deriv(0,0) =  1;
    deriv(1,0) = -1;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear replacing quad, line3)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (only form a basis and have to be adapted for distorted elements)
  // *********************************************************************
  case MortarElement::dual1D_base_for_edge1:
  {
    val[0] = -xi[0];
    val[1] = 1+xi[0];
    deriv(0,0) = -1;
    deriv(1,0) =  1;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear replacing quad, line3)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case MortarElement::quaddual1D_edge0:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();

    // empty shape function vals + derivs
    LINALG::SerialDenseVector valquad(nnodes);
    LINALG::SerialDenseMatrix derivquad(nnodes,1);
    LINALG::SerialDenseVector vallin(nnodes-1);
    LINALG::SerialDenseMatrix derivlin(nnodes-1,1);
    LINALG::SerialDenseVector valtemp(nnodes,true);
    LINALG::SerialDenseMatrix derivtemp(nnodes,1,true);

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::Integrator integrator(Shape());

    LINALG::SerialDenseMatrix me(nnodes-1,nnodes-1,true);
    LINALG::SerialDenseMatrix de(nnodes-1,nnodes-1,true);

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), 0.0};
      ShapeFunctions(MortarElement::quad1D,gpc,valquad,derivquad);
      ShapeFunctions(MortarElement::dual1D_base_for_edge0,gpc,vallin,derivlin);
      detg = Jacobian(gpc);

      for (int j=1;j<nnodes;++j)
        for (int k=1;k<nnodes;++k)
        {
          me(j-1,k-1)+=integrator.Weight(i)*vallin[j-1]*valquad[k]*detg;
          de(j-1,k-1)+=(j==k)*integrator.Weight(i)*valquad[k]*detg;
        }
    }

    // invert bi-ortho matrix me
    // CAUTION: This is a non-symmetric inverse operation!
    double detme = me(0,0)*me(1,1)-me(0,1)*me(1,0);
    LINALG::SerialDenseMatrix meold(nnodes-1,nnodes-1);
    meold=me;
    me(0,0) =  1/detme*meold(1,1);
    me(0,1) = -1/detme*meold(0,1);
    me(1,0) = -1/detme*meold(1,0);
    me(1,1) =  1/detme*meold(0,0);

    // get solution matrix with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes-1,nnodes-1);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    ShapeFunctions(MortarElement::dual1D_base_for_edge0,xi,vallin,derivlin);
    for (int i=1;i<nnodes;++i)
      for (int j=1;j<nnodes;++j)
      {
        valtemp[i]+=ae(i-1,j-1)*vallin[j-1];
        derivtemp(i,0)+=ae(i-1,j-1)*derivlin(j-1,0);
      }

    val[0] = 0.0;
    val[1] = valtemp[1];
    val[2] = valtemp[2];
    deriv(0,0) =  0.0;
    deriv(1,0) = derivtemp(1,0);
    deriv(2,0) = derivtemp(2,0);

    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear replacing quad, line3)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case MortarElement::quaddual1D_edge1:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();

    // empty shape function vals + derivs
    LINALG::SerialDenseVector valquad(nnodes);
    LINALG::SerialDenseMatrix derivquad(nnodes,1);
    LINALG::SerialDenseVector vallin(nnodes-1);
    LINALG::SerialDenseMatrix derivlin(nnodes-1,1);
    LINALG::SerialDenseVector valtemp(nnodes,true);
    LINALG::SerialDenseMatrix derivtemp(nnodes,1,true);

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::Integrator integrator(Shape());

    LINALG::SerialDenseMatrix me(nnodes-1,nnodes-1,true);
    LINALG::SerialDenseMatrix de(nnodes-1,nnodes-1,true);

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), 0.0};
      ShapeFunctions(MortarElement::quad1D,gpc,valquad,derivquad);
      ShapeFunctions(MortarElement::dual1D_base_for_edge1,gpc,vallin,derivlin);
      detg = Jacobian(gpc);

      for (int j=0;j<nnodes-1;++j)
        for (int k=0;k<nnodes-1;++k)
        {
          me(j,k)+=integrator.Weight(i)*vallin[j]*valquad[2*k]*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*valquad[2*k]*detg;
        }
    }

    // invert bi-ortho matrix me
    // CAUTION: This is a non-symmetric inverse operation!
    double detme = me(0,0)*me(1,1)-me(0,1)*me(1,0);
    LINALG::SerialDenseMatrix meold(nnodes-1,nnodes-1);
    meold=me;
    me(0,0) =  1/detme*meold(1,1);
    me(0,1) = -1/detme*meold(0,1);
    me(1,0) = -1/detme*meold(1,0);
    me(1,1) =  1/detme*meold(0,0);

    // get solution matrix with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes-1,nnodes-1);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    ShapeFunctions(MortarElement::dual1D_base_for_edge1,xi,vallin,derivlin);
    for (int i=0;i<nnodes-1;++i)
      for (int j=0;j<nnodes-1;++j)
      {
        valtemp[2*i]+=ae(i,j)*vallin[j];
        derivtemp(2*i,0)+=ae(i,j)*derivlin(j,0);
      }

    val[0] = valtemp[0];
    val[1] = 0.0;
    val[2] = valtemp[2];
    deriv(0,0) = derivtemp(0,0);
    deriv(1,0) = 0.0;
    deriv(2,0) = derivtemp(2,0);

    break;
  }
  // *********************************************************************
  // Unkown shape function type
  // *********************************************************************
  default:
    dserror("ERROR: Unknown shape function type identifier");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate shape functions                                  popp 01/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                                      LINALG::SerialDenseMatrix& deriv, const int& valdim)
{
  if (!xi)
    dserror("ERROR: EvaluateShape called with xi=NULL");

  switch(Shape())
  {
  // 2D linear case (2noded line element)
  case DRT::Element::line2:
  {
    if (valdim!=2) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::lin1D,xi,val,deriv);
  break;
  }
  // 2D quadratic case (3noded line element)
  case DRT::Element::line3:
  {
    if (valdim!=3) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::quad1D,xi,val,deriv);
    break;
  }
  // 3D linear case (3noded triangular element)
  case DRT::Element::tri3:
  {
    if (valdim!=3) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::lin2D,xi,val,deriv);
    break;
  }
  // 3D bilinear case (4noded quadrilateral element)
  case DRT::Element::quad4:
  {
    if (valdim!=4) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::bilin2D,xi,val,deriv);
    break;
  }
  // 3D quadratic case (6noded triangular element)
  case DRT::Element::tri6:
  {
    if (valdim!=6) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::quad2D,xi,val,deriv);
    break;
  }
  // 3D serendipity case (8noded quadrilateral element)
  case DRT::Element::quad8:
  {
    if (valdim!=8) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::serendipity2D,xi,val,deriv);
    break;
  }
  // 3D biquadratic case (9noded quadrilateral element)
  case DRT::Element::quad9:
  {
    if (valdim!=9) dserror("ERROR: Inconsistency in EvluateShape");
    ShapeFunctions(MortarElement::biquad2D,xi,val,deriv);
    break;
  }
  // unknown case
  default:
    dserror("ERROR: EvaluateShape called for unknown MortarElement type");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate dual shape functions                             popp 12/07|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::EvaluateShapeDual(const double* xi, LINALG::SerialDenseVector& val,
                                          LINALG::SerialDenseMatrix& deriv, const int& valdim)
{
  if (!xi)
    dserror("ERROR: EvaluateShapeDual called with xi=NULL");
  if (!IsSlave())
    dserror("ERROR: EvaluateShapeDual called for master element");

  // get node number and node pointers
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: EvaluateShapeDual: Null pointer!");

  switch(Shape())
  {
  // 2D linear case (2noded line element)
  case DRT::Element::line2:
  {
    if (valdim!=2) dserror("ERROR: Inconsistency in EvluateShape");
    // check for boundary nodes
    MortarNode* mymrtrnode0 = static_cast<MortarNode*> (mynodes[0]);
    MortarNode* mymrtrnode1 = static_cast<MortarNode*> (mynodes[1]);
    if (!mymrtrnode0) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    if (!mymrtrnode1) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    bool isonbound0 = mymrtrnode0->IsOnBound();
    bool isonbound1 = mymrtrnode1->IsOnBound();

    // both nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1)
      ShapeFunctions(MortarElement::lindual1D,xi,val,deriv);

    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1)
      ShapeFunctions(MortarElement::lindual1D_edge0,xi,val,deriv);

    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1)
      ShapeFunctions(MortarElement::lindual1D_edge1,xi,val,deriv);

    // both nodes are on boundary: infeasible case
    else
      dserror("ERROR: EvaluateShapeDual1D: Element with 2 boundary nodes!");

    break;
  }

  // 2D quadratic case (3noded line element)
  case DRT::Element::line3:
  {
    if (valdim!=3) dserror("ERROR: Inconsistency in EvluateShape");
    // check for boundary nodes
    MortarNode* mymrtrnode0 = static_cast<MortarNode*> (mynodes[0]);
    MortarNode* mymrtrnode1 = static_cast<MortarNode*> (mynodes[1]);
    MortarNode* mymrtrnode2 = static_cast<MortarNode*> (mynodes[2]);
    if (!mymrtrnode0) dserror("ERROR: EvaluateShapeDual: Null pointer!");
    if (!mymrtrnode1) dserror("ERROR: EvaluateShapeDual: Null pointer!");
    if (!mymrtrnode2) dserror("ERROR: EvaluateShapeDual: Null pointer!");
    bool isonbound0 = mymrtrnode0->IsOnBound();
    bool isonbound1 = mymrtrnode1->IsOnBound();
    bool isonbound2 = mymrtrnode2->IsOnBound();

    // all 3 nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1 && !isonbound2)
      ShapeFunctions(MortarElement::quaddual1D,xi,val,deriv);

    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1 && !isonbound2)
      ShapeFunctions(MortarElement::quaddual1D_edge0,xi,val,deriv);

    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1 && !isonbound2)
      ShapeFunctions(MortarElement::quaddual1D_edge1,xi,val,deriv);

    // node 2 is on boundary: infeasible case
    else if (isonbound2)
      dserror("ERROR: EvaluateShapeDual: Middle boundary node");

    // nodes 0 and 1 are on boundary: infeasible case
    else
      dserror("ERROR: EvaluateShapeDual: Element with 2 boundary nodes");

    break;
  }

  // 3D cases
  case DRT::Element::tri3:
  case DRT::Element::quad4:
  case DRT::Element::tri6:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    // check for boundary nodes
    bool bound = false;
    for (int i=0;i<NumNode();++i)
    {
      MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);
      if (!mymrtrnode) dserror("ERROR: EvaluateShapeDual: Null pointer!");
      bound += mymrtrnode->IsOnBound();
    }

    // all nodes are interior: use unmodified dual shape functions
    if (!bound)
    {
      if (Shape()==tri3)       ShapeFunctions(MortarElement::lindual2D,xi,val,deriv);
      else if (Shape()==quad4) ShapeFunctions(MortarElement::bilindual2D,xi,val,deriv);
      else if (Shape()==tri6)  ShapeFunctions(MortarElement::quaddual2D,xi,val,deriv);
      else if (Shape()==quad8) ShapeFunctions(MortarElement::serendipitydual2D,xi,val,deriv);
      else /*Shape()==quad9*/  ShapeFunctions(MortarElement::biquaddual2D,xi,val,deriv);
    }

    // some nodes are on slave boundary
    else
      dserror("ERROR: EvaluateShapeDual: boundary mod. not yet impl. for 3D contact!");

    break;
  }

  // unknown case
  default:
    dserror("ERROR: EvaluateShapeDual called for unknown element type");
  }

  return true;
}

#endif  // #ifdef CCADISCRET
