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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "mortar_element.H"
#include "mortar_node.H"
#include "mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

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
  // 1D modified standard shape functions (const replacing linear, line2)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case MortarElement::lin1D_edge0:
  {
    val[0] = 0;
    val[1] = 1;
    deriv(0,0) = 0;
    deriv(1,0) = 0;
    break;
  }
  // *********************************************************************
  // 1D modified standard shape functions (const replacing linear, line2)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case MortarElement::lin1D_edge1:
  {
    val[0] = 1;
    val[1] = 0;
    deriv(0,0) = 0;
    deriv(1,0) = 0;
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
  // 2D standard bilinear shape functions (quad4)
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
  // 1D modified standard shape functions (linear replacing quad, line3)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case MortarElement::quad1D_edge0:
  {
    val[0] = 0;
    val[1] = xi[0];
    val[2] = 1-xi[0];
    deriv(0,0) =  0;
    deriv(1,0) =  1;
    deriv(2,0) = -1;
    break;
  }
  // *********************************************************************
  // 1D modified standard shape functions (linear replacing quad, line3)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case MortarElement::quad1D_edge1:
  {
    val[0] = -xi[0];
    val[1] = 0;
    val[2] = 1+xi[0];
    deriv(0,0) = -1;
    deriv(1,0) =  0;
    deriv(2,0) =  1;
    break;
  }
  // *********************************************************************
  // 1D linear part of standard quadratic shape functions (line3)
  // (used for linear interpolation of std LM field in 2D quadratic mortar)
  // *********************************************************************
  case MortarElement::quad1D_only_lin:
  {
    val[0] = 0.5*(1-xi[0]);
    val[1] = 0.5*(1+xi[0]);
    val[2] = 0.0;
    deriv(0,0) = -0.5;
    deriv(1,0) =  0.5;
    deriv(2,0) =  0.0;
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
  // 2D modified quadratic shape functions (tri6)
  // (used in combination with quadr dual LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::quad2D_modified:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double t1 = 1.0-r-s;
    const double t2 = r;
    const double t3 = s;

    LINALG::SerialDenseVector valtmp(NumNode(),1);
    LINALG::SerialDenseMatrix derivtmp(NumNode(),2);

    valtmp[0] = t1*(2.0*t1-1.0);
    valtmp[1] = t2*(2.0*t2-1.0);
    valtmp[2] = t3*(2.0*t3-1.0);
    valtmp[3] = 4.0*t2*t1;
    valtmp[4] = 4.0*t2*t3;
    valtmp[5] = 4.0*t3*t1;

    derivtmp(0,0)= -3.0+4.0*(r+s);
    derivtmp(0,1)= -3.0+4.0*(r+s);
    derivtmp(1,0)= 4.0*r-1.0;
    derivtmp(1,1)= 0.0;
    derivtmp(2,0)= 0.0;
    derivtmp(2,1)= 4.0*s-1.0;
    derivtmp(3,0)= 4.0*(1.0-2.0*r-s);
    derivtmp(3,1)=-4.0*r;
    derivtmp(4,0)= 4.0*s;
    derivtmp(4,1)= 4.0*r;
    derivtmp(5,0)=-4.0*s;
    derivtmp(5,1)= 4.0*(1.0-r-2.0*s);

    // define constant modification factor 1/5
    // (NOTE: lower factors, e.g. 1/12 would be sufficient here
    // as well, but in order to be globally continuous for mixed
    // meshes with tet10/hex20 elements, we always choose 1/5.)
    double fac = 1.0/5.0;

    // apply constant modification at vertex nodes and PoU
    val[0] = valtmp[0]+(valtmp[3]+valtmp[5])*fac;
    val[1] = valtmp[1]+(valtmp[3]+valtmp[4])*fac;
    val[2] = valtmp[2]+(valtmp[4]+valtmp[5])*fac;
    val[3] = valtmp[3]*(1.0-2.0*fac);
    val[4] = valtmp[4]*(1.0-2.0*fac);
    val[5] = valtmp[5]*(1.0-2.0*fac);

    deriv(0,0)= derivtmp(0,0)+(derivtmp(3,0)+derivtmp(5,0))*fac;
    deriv(0,1)= derivtmp(0,1)+(derivtmp(3,1)+derivtmp(5,1))*fac;
    deriv(1,0)= derivtmp(1,0)+(derivtmp(3,0)+derivtmp(4,0))*fac;
    deriv(1,1)= derivtmp(1,1)+(derivtmp(3,1)+derivtmp(4,1))*fac;
    deriv(2,0)= derivtmp(2,0)+(derivtmp(4,0)+derivtmp(5,0))*fac;
    deriv(2,1)= derivtmp(2,1)+(derivtmp(4,1)+derivtmp(5,1))*fac;
    deriv(3,0)= derivtmp(3,0)*(1.0-2.0*fac);
    deriv(3,1)= derivtmp(3,1)*(1.0-2.0*fac);
    deriv(4,0)= derivtmp(4,0)*(1.0-2.0*fac);
    deriv(4,1)= derivtmp(4,1)*(1.0-2.0*fac);
    deriv(5,0)= derivtmp(5,0)*(1.0-2.0*fac);
    deriv(5,1)= derivtmp(5,1)*(1.0-2.0*fac);

    break;
  }
  // *********************************************************************
  // 2D modified (hierarchical) quadratic shape functions (tri6)
  // (used in combination with linear dual LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::quad2D_hierarchical:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double t1 = 1.0-r-s;
    const double t2 = r;
    const double t3 = s;

    val[0] = t1;
    val[1] = t2;
    val[2] = t3;
    val[3] = 4.0*t2*t1;
    val[4] = 4.0*t2*t3;
    val[5] = 4.0*t3*t1;

    deriv(0,0)= -1.0;
    deriv(0,1)= -1.0;
    deriv(1,0)=  1.0;
    deriv(1,1)=  0.0;
    deriv(2,0)=  0.0;
    deriv(2,1)=  1.0;
    deriv(3,0)= 4.0*(1.0-2.0*r-s);
    deriv(3,1)=-4.0*r;
    deriv(4,0)= 4.0*s;
    deriv(4,1)= 4.0*r;
    deriv(5,0)=-4.0*s;
    deriv(5,1)= 4.0*(1.0-r-2.0*s);

    break;
  }
  // *********************************************************************
  // 2D linear part of standard quadratic shape functions (tri6)
  // (used for linear interpolation of std LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::quad2D_only_lin:
  {
    val[0] = 1-xi[0]-xi[1];
    val[1] = xi[0];
    val[2] = xi[1];
    val[3] = 0.0;
    val[4] = 0.0;
    val[5] = 0.0;
    
    deriv(0,0) = -1.0; deriv(0,1) = -1.0;
    deriv(1,0) =  1.0; deriv(1,1) =  0.0;
    deriv(2,0) =  0.0; deriv(2,1) =  1.0;
    deriv(3,0) =  0.0; deriv(3,1) =  0.0;
    deriv(4,0) =  0.0; deriv(4,1) =  0.0;
    deriv(5,0) =  0.0; deriv(5,1) =  0.0;
        
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
  // 2D modified serendipity shape functions (quad8)
  // (used in combination with quadr dual LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::serendipity2D_modified:
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

    LINALG::SerialDenseVector valtmp(NumNode(),1);
    LINALG::SerialDenseMatrix derivtmp(NumNode(),2);

    valtmp[0]=0.25*(rm*sm-(r2*sm+s2*rm));
    valtmp[1]=0.25*(rp*sm-(r2*sm+s2*rp));
    valtmp[2]=0.25*(rp*sp-(s2*rp+r2*sp));
    valtmp[3]=0.25*(rm*sp-(r2*sp+s2*rm));
    valtmp[4]=0.5*r2*sm;
    valtmp[5]=0.5*s2*rp;
    valtmp[6]=0.5*r2*sp;
    valtmp[7]=0.5*s2*rm;

    derivtmp(0,0)= 0.25*sm*(2*r+s);
    derivtmp(0,1)= 0.25*rm*(r+2*s);
    derivtmp(1,0)= 0.25*sm*(2*r-s);
    derivtmp(1,1)= 0.25*rp*(2*s-r);
    derivtmp(2,0)= 0.25*sp*(2*r+s);
    derivtmp(2,1)= 0.25*rp*(r+2*s);
    derivtmp(3,0)= 0.25*sp*(2*r-s);
    derivtmp(3,1)= 0.25*rm*(2*s-r);
    derivtmp(4,0)=-sm*r;
    derivtmp(4,1)=-0.5*rm*rp;
    derivtmp(5,0)= 0.5*sm*sp;
    derivtmp(5,1)=-rp*s;
    derivtmp(6,0)=-sp*r;
    derivtmp(6,1)= 0.5*rm*rp;
    derivtmp(7,0)=-0.5*sm*sp;
    derivtmp(7,1)=-rm*s;

    // define constant modification factor 1/5
    double fac = 1.0/5.0;

    // apply constant modification at vertex nodes and PoU
    val[0]=valtmp[0]+(valtmp[4]+valtmp[7])*fac;
    val[1]=valtmp[1]+(valtmp[4]+valtmp[5])*fac;
    val[2]=valtmp[2]+(valtmp[5]+valtmp[6])*fac;
    val[3]=valtmp[3]+(valtmp[6]+valtmp[7])*fac;
    val[4]=valtmp[4]*(1.0-2.0*fac);
    val[5]=valtmp[5]*(1.0-2.0*fac);
    val[6]=valtmp[6]*(1.0-2.0*fac);
    val[7]=valtmp[7]*(1.0-2.0*fac);

    deriv(0,0)=derivtmp(0,0)+(derivtmp(4,0)+derivtmp(7,0))*fac;
    deriv(0,1)=derivtmp(0,1)+(derivtmp(4,1)+derivtmp(7,1))*fac;
    deriv(1,0)=derivtmp(1,0)+(derivtmp(4,0)+derivtmp(5,0))*fac;
    deriv(1,1)=derivtmp(1,1)+(derivtmp(4,1)+derivtmp(5,1))*fac;
    deriv(2,0)=derivtmp(2,0)+(derivtmp(5,0)+derivtmp(6,0))*fac;
    deriv(2,1)=derivtmp(2,1)+(derivtmp(5,1)+derivtmp(6,1))*fac;
    deriv(3,0)=derivtmp(3,0)+(derivtmp(6,0)+derivtmp(7,0))*fac;
    deriv(3,1)=derivtmp(3,1)+(derivtmp(6,1)+derivtmp(7,1))*fac;
    deriv(4,0)=derivtmp(4,0)*(1.0-2.0*fac);
    deriv(4,1)=derivtmp(4,1)*(1.0-2.0*fac);
    deriv(5,0)=derivtmp(5,0)*(1.0-2.0*fac);
    deriv(5,1)=derivtmp(5,1)*(1.0-2.0*fac);
    deriv(6,0)=derivtmp(6,0)*(1.0-2.0*fac);
    deriv(6,1)=derivtmp(6,1)*(1.0-2.0*fac);
    deriv(7,0)=derivtmp(7,0)*(1.0-2.0*fac);
    deriv(7,1)=derivtmp(7,1)*(1.0-2.0*fac);

    break;
  }
  // *********************************************************************
  // 2D modified (hierarchical) serendipity shape functions (quad8)
  // (used in combination with linear dual LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::serendipity2D_hierarchical:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double rp=1.0+r;
    const double rm=1.0-r;
    const double sp=1.0+s;
    const double sm=1.0-s;
    const double r2=1.0-r*r;
    const double s2=1.0-s*s;

    val[0]=0.25*rm*sm;
    val[1]=0.25*rp*sm;
    val[2]=0.25*rp*sp;
    val[3]=0.25*rm*sp;
    val[4]=0.5*r2*sm;
    val[5]=0.5*s2*rp;
    val[6]=0.5*r2*sp;
    val[7]=0.5*s2*rm;

    deriv(0,0)=-0.25*sm;
    deriv(0,1)=-0.25*rm;
    deriv(1,0)= 0.25*sm;
    deriv(1,1)=-0.25*rp;
    deriv(2,0)= 0.25*sp;
    deriv(2,1)= 0.25*rp;
    deriv(3,0)=-0.25*sp;
    deriv(3,1)= 0.25*rm;
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
  // 2D bilinear part of serendipity quadratic shape functions (quad8)
  // (used for linear interpolation of std LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::serendipity2D_only_lin:
  {
    val[0] = 0.25*(1-xi[0])*(1-xi[1]);
    val[1] = 0.25*(1+xi[0])*(1-xi[1]);
    val[2] = 0.25*(1+xi[0])*(1+xi[1]);
    val[3] = 0.25*(1-xi[0])*(1+xi[1]);
    val[4] = 0.0;
    val[5] = 0.0;
    val[6] = 0.0;
    val[7] = 0.0;
    
    deriv(0,0) = -0.25*(1-xi[1]); deriv(0,1) = -0.25*(1-xi[0]);
    deriv(1,0) =  0.25*(1-xi[1]); deriv(1,1) = -0.25*(1+xi[0]);
    deriv(2,0) =  0.25*(1+xi[1]); deriv(2,1) =  0.25*(1+xi[0]);
    deriv(3,0) = -0.25*(1+xi[1]); deriv(3,1) =  0.25*(1-xi[0]);
    deriv(4,0) =  0.0;            deriv(4,1) =  0.0;
    deriv(5,0) =  0.0;            deriv(5,1) =  0.0;
    deriv(6,0) =  0.0;            deriv(6,1) =  0.0;
    deriv(7,0) =  0.0;            deriv(7,1) =  0.0;
    
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
  // 2D standard biquadratic shape functions (quad9)
  // (used in combination with quadr dual LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::biquad2D_modified:
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

    LINALG::SerialDenseVector valtmp(NumNode(),1);
    LINALG::SerialDenseMatrix derivtmp(NumNode(),2);
    
    valtmp[0]= rs*rm*sm;
    valtmp[1]=-rs*rp*sm;
    valtmp[2]= rs*rp*sp;
    valtmp[3]=-rs*rm*sp;
    valtmp[4]=-sh*sm*r2;
    valtmp[5]= rh*rp*s2;
    valtmp[6]= sh*sp*r2;
    valtmp[7]=-rh*rm*s2;
    valtmp[8]= r2*s2;

    derivtmp(0,0)=-rhm*sh*sm;
    derivtmp(0,1)=-shm*rh*rm;
    derivtmp(1,0)=-rhp*sh*sm;
    derivtmp(1,1)= shm*rh*rp;
    derivtmp(2,0)= rhp*sh*sp;
    derivtmp(2,1)= shp*rh*rp;
    derivtmp(3,0)= rhm*sh*sp;
    derivtmp(3,1)=-shp*rh*rm;
    derivtmp(4,0)= 2.0*r*sh*sm;
    derivtmp(4,1)= shm*r2;
    derivtmp(5,0)= rhp*s2;
    derivtmp(5,1)=-2.0*s*rh*rp;
    derivtmp(6,0)=-2.0*r*sh*sp;
    derivtmp(6,1)= shp*r2;
    derivtmp(7,0)= rhm*s2;
    derivtmp(7,1)= 2.0*s*rh*rm;
    derivtmp(8,0)=-2.0*r*s2;
    derivtmp(8,1)=-2.0*s*r2;
    
    // define constant modification factor
    // (CURRENTLY NOT USED -> ZERO)
    double fac = 0.0;

    // apply constant modification at vertex nodes and PoU
    val[0]=valtmp[0]+(valtmp[4]+valtmp[7])*fac+0.5*valtmp[8]*fac;
    val[1]=valtmp[1]+(valtmp[4]+valtmp[5])*fac+0.5*valtmp[8]*fac;
    val[2]=valtmp[2]+(valtmp[5]+valtmp[6])*fac+0.5*valtmp[8]*fac;
    val[3]=valtmp[3]+(valtmp[6]+valtmp[7])*fac+0.5*valtmp[8]*fac;
    val[4]=valtmp[4]*(1.0-2.0*fac);
    val[5]=valtmp[5]*(1.0-2.0*fac);
    val[6]=valtmp[6]*(1.0-2.0*fac);
    val[7]=valtmp[7]*(1.0-2.0*fac);
    val[8]=valtmp[8]*(1.0-4.0*0.5*fac);

    deriv(0,0)=derivtmp(0,0)+(derivtmp(4,0)+derivtmp(7,0))*fac+0.5*derivtmp(8,0)*fac;
    deriv(0,1)=derivtmp(0,1)+(derivtmp(4,1)+derivtmp(7,1))*fac+0.5*derivtmp(8,1)*fac;
    deriv(1,0)=derivtmp(1,0)+(derivtmp(4,0)+derivtmp(5,0))*fac+0.5*derivtmp(8,0)*fac;
    deriv(1,1)=derivtmp(1,1)+(derivtmp(4,1)+derivtmp(5,1))*fac+0.5*derivtmp(8,1)*fac;
    deriv(2,0)=derivtmp(2,0)+(derivtmp(5,0)+derivtmp(6,0))*fac+0.5*derivtmp(8,0)*fac;
    deriv(2,1)=derivtmp(2,1)+(derivtmp(5,1)+derivtmp(6,1))*fac+0.5*derivtmp(8,1)*fac;
    deriv(3,0)=derivtmp(3,0)+(derivtmp(6,0)+derivtmp(7,0))*fac+0.5*derivtmp(8,0)*fac;
    deriv(3,1)=derivtmp(3,1)+(derivtmp(6,1)+derivtmp(7,1))*fac+0.5*derivtmp(8,1)*fac;
    deriv(4,0)=derivtmp(4,0)*(1.0-2.0*fac);
    deriv(4,1)=derivtmp(4,1)*(1.0-2.0*fac);
    deriv(5,0)=derivtmp(5,0)*(1.0-2.0*fac);
    deriv(5,1)=derivtmp(5,1)*(1.0-2.0*fac);
    deriv(6,0)=derivtmp(6,0)*(1.0-2.0*fac);
    deriv(6,1)=derivtmp(6,1)*(1.0-2.0*fac);
    deriv(7,0)=derivtmp(7,0)*(1.0-2.0*fac);
    deriv(7,1)=derivtmp(7,1)*(1.0-2.0*fac);
    deriv(8,0)=derivtmp(8,0)*(1.0-4.0*0.5*fac);
    deriv(8,1)=derivtmp(8,1)*(1.0-4.0*0.5*fac);

    break;
  }
  // *********************************************************************
  // 2D standard biquadratic shape functions (quad9)
  // (used in combination with linear dual LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::biquad2D_hierarchical:
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
    //const double rs=rh*sh;
    const double rhp=r+0.5;
    const double rhm=r-0.5;
    const double shp=s+0.5;
    const double shm=s-0.5;

    val[0]=0.25*rm*sm;
    val[1]=0.25*rp*sm;
    val[2]=0.25*rp*sp;
    val[3]=0.25*rm*sp;
    val[4]=-sh*sm*r2;
    val[5]= rh*rp*s2;
    val[6]= sh*sp*r2;
    val[7]=-rh*rm*s2;
    val[8]= r2*s2;

    deriv(0,0)=-0.25*sm;
    deriv(0,1)=-0.25*rm;
    deriv(1,0)= 0.25*sm;
    deriv(1,1)=-0.25*rp;
    deriv(2,0)= 0.25*sp;
    deriv(2,1)= 0.25*rp;
    deriv(3,0)=-0.25*sp;
    deriv(3,1)= 0.25*rm;
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
  // 2D bilinear part of biquadratic quadratic shape functions (quad9)
  // (used for linear interpolation of std LM field in 3D quadratic mortar)
  // *********************************************************************
  case MortarElement::biquad2D_only_lin:
  {
    val[0] = 0.25*(1-xi[0])*(1-xi[1]);
    val[1] = 0.25*(1+xi[0])*(1-xi[1]);
    val[2] = 0.25*(1+xi[0])*(1+xi[1]);
    val[3] = 0.25*(1-xi[0])*(1+xi[1]);
    val[4] = 0.0;
    val[5] = 0.0;
    val[6] = 0.0;
    val[7] = 0.0;
    val[8] = 0.0;
    
    deriv(0,0) = -0.25*(1-xi[1]); deriv(0,1) = -0.25*(1-xi[0]);
    deriv(1,0) =  0.25*(1-xi[1]); deriv(1,1) = -0.25*(1+xi[0]);
    deriv(2,0) =  0.25*(1+xi[1]); deriv(2,1) =  0.25*(1+xi[0]);
    deriv(3,0) = -0.25*(1+xi[1]); deriv(3,1) =  0.25*(1-xi[0]);
    deriv(4,0) =  0.0;            deriv(4,1) =  0.0;
    deriv(5,0) =  0.0;            deriv(5,1) =  0.0;
    deriv(6,0) =  0.0;            deriv(6,1) =  0.0;
    deriv(7,0) =  0.0;            deriv(7,1) =  0.0;
    deriv(8,0) =  0.0;            deriv(8,1) =  0.0;
    
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
  // (used for interpolation of Lagrange mutliplier field)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case MortarElement::quaddual1D:
  case MortarElement::bilindual2D:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());

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

    // check whether this is a 1D or 2D mortar element
    int dim = 2;
    if (shape==MortarElement::quaddual1D) dim = 1;

    // evaluate dual shape functions
    LINALG::SerialDenseVector valtemp(nnodes,true);
    LINALG::SerialDenseMatrix derivtemp(nnodes,dim,true);
    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        valtemp[i]+=ae(i,j)*val[j];
        derivtemp(i,0)+=ae(i,j)*deriv(j,0);
        if (dim==2) derivtemp(i,1)+=ae(i,j)*deriv(j,1);
      }

    val=valtemp;
    deriv=derivtemp;

    break;
  }
  // *********************************************************************
  // 2D dual quadratic shape functions (tri6)
  // 2D dual serendipity shape functions (quad8)
  // 2D dual biquadratic shape functions (quad9)
  // (used for interpolation of Lagrange mutliplier field)
  // (including adaption process for distorted elements)
  // (including modification of displacement shape functions)
  // *********************************************************************
  case MortarElement::quaddual2D:
  case MortarElement::serendipitydual2D:
  case MortarElement::biquaddual2D:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();

    // empty shape function vals + derivs
    LINALG::SerialDenseVector valquad(nnodes);
    LINALG::SerialDenseMatrix derivquad(nnodes,2);

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());
    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
      EvaluateShape(gpc, valquad, derivquad, nnodes, true);
      detg = Jacobian(gpc);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=integrator.Weight(i)*valquad[j]*valquad[k]*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*valquad[j]*detg;
        }
    }

    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);

    // get solution matrix with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    EvaluateShape(xi, valquad, derivquad, nnodes, true);
    val.Zero(); deriv.Zero();

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        val[i]+=ae(i,j)*valquad[j];
        deriv(i,0)+=ae(i,j)*derivquad(j,0);
        deriv(i,1)+=ae(i,j)*derivquad(j,1);
      }

    break;
  }
  // *********************************************************************
  // 2D dual quadratic shape functions (tri6)
  // 2D dual serendipity shape functions (quad8)
  // 2D dual biquadratic shape functions (quad9)
  // (used for LINEAR interpolation of Lagrange mutliplier field)
  // (including adaption process for distorted elements)
  // (including modification of displacement shape functions)
  // *********************************************************************
  case MortarElement::quaddual2D_only_lin:
  case MortarElement::serendipitydual2D_only_lin:
  case MortarElement::biquaddual2D_only_lin:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();

    // empty shape function vals + derivs
    LINALG::SerialDenseVector valquad(nnodes);
    LINALG::SerialDenseMatrix derivquad(nnodes,2);

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());
    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
      if (shape==quaddual2D_only_lin)
        ShapeFunctions(MortarElement::quad2D_only_lin,gpc,valquad,derivquad);
      else if (shape==serendipitydual2D_only_lin)
        ShapeFunctions(MortarElement::serendipity2D_only_lin,gpc,valquad,derivquad);
      else if (shape==biquaddual2D_only_lin)
        ShapeFunctions(MortarElement::biquad2D_only_lin,gpc,valquad,derivquad);
      detg = Jacobian(gpc);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=integrator.Weight(i)*valquad[j]*valquad[k]*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*valquad[j]*detg;
        }
    }

    // how many non-zero nodes
    int nnodeslin = 0;
    if (shape==quaddual2D_only_lin)             nnodeslin = 3;
    else if (shape==serendipitydual2D_only_lin) nnodeslin = 4;
    else if (shape==biquaddual2D_only_lin)      nnodeslin = 4;
    
    // reduce me to non-zero nodes before inverting
    LINALG::SerialDenseMatrix melin(nnodeslin,nnodeslin,true);
    for (int j=0;j<nnodeslin;++j)
      for (int k=0;k<nnodeslin;++k)
        melin(j,k) = me(j,k);
    
    // invert bi-ortho matrix melin
    LINALG::SymmetricInverse(melin,nnodeslin);

    // re-inflate inverse of melin to full size
    LINALG::SerialDenseMatrix invme(nnodes,nnodes,true);
    for (int j=0;j<nnodeslin;++j)
      for (int k=0;k<nnodeslin;++k)
        invme(j,k) = melin(j,k);
    
    // get solution matrix with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);
    ae.Multiply('N','N',1.0,de,invme,0.0);

    // evaluate dual shape functions at loc. coord. xi
    if (shape==quaddual2D_only_lin)
      ShapeFunctions(MortarElement::quad2D_only_lin,xi,valquad,derivquad);
    else if (shape==serendipitydual2D_only_lin)
      ShapeFunctions(MortarElement::serendipity2D_only_lin,xi,valquad,derivquad);
    else if (shape==biquaddual2D_only_lin)
      ShapeFunctions(MortarElement::biquad2D_only_lin,xi,valquad,derivquad);
    val.Zero(); deriv.Zero();

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        val[i]+=ae(i,j)*valquad[j];
        deriv(i,0)+=ae(i,j)*derivquad(j,0);
        deriv(i,1)+=ae(i,j)*derivquad(j,1);
      }

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
    MORTAR::ElementIntegrator integrator(Shape());

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
    MORTAR::ElementIntegrator integrator(Shape());

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
 |  Evaluate displacement shape functions                     popp 01/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                                          LINALG::SerialDenseMatrix& deriv, const int& valdim,
                                          bool dualquad3d)
{
  if (!xi)
    dserror("ERROR: EvaluateShape called with xi=NULL");

  // get node number and node pointers
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
  
  // check for boundary nodes
  bool bound = false;
  for (int i=0;i<NumNode();++i)
  {
    MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);
    if (!mymrtrnode) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }
  
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
    if (dualquad3d && !bound)     ShapeFunctions(MortarElement::quad2D_modified,xi,val,deriv);
    else if (dualquad3d && bound) ShapeFunctions(MortarElement::quad2D_hierarchical,xi,val,deriv);
    else                          ShapeFunctions(MortarElement::quad2D,xi,val,deriv);
    break;
  }
  // 3D serendipity case (8noded quadrilateral element)
  case DRT::Element::quad8:
  {
    if (valdim!=8) dserror("ERROR: Inconsistency in EvluateShape");
    if (dualquad3d && !bound)     ShapeFunctions(MortarElement::serendipity2D_modified,xi,val,deriv);
    else if (dualquad3d && bound) ShapeFunctions(MortarElement::serendipity2D_hierarchical,xi,val,deriv);
    else                          ShapeFunctions(MortarElement::serendipity2D,xi,val,deriv);
    break;
  }
  // 3D biquadratic case (9noded quadrilateral element)
  case DRT::Element::quad9:
  {
    if (valdim!=9) dserror("ERROR: Inconsistency in EvluateShape");
    if (dualquad3d && !bound)     ShapeFunctions(MortarElement::biquad2D_modified,xi,val,deriv);
    else if (dualquad3d && bound) ShapeFunctions(MortarElement::biquad2D_hierarchical,xi,val,deriv);
    else                          ShapeFunctions(MortarElement::biquad2D,xi,val,deriv);
    break;
  }
  // unknown case
  default:
    dserror("ERROR: EvaluateShape called for unknown MortarElement type");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate Lagrange multiplier shape functions              popp 12/07|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::EvaluateShapeLagMult(const INPAR::MORTAR::ShapeFcn& lmtype,
                                                 const double* xi, LINALG::SerialDenseVector& val,
                                                 LINALG::SerialDenseMatrix& deriv, const int& valdim)
{
  if (!xi)
    dserror("ERROR: EvaluateShapeLagMult called with xi=NULL");
  if (!IsSlave())
    dserror("ERROR: EvaluateShapeLagMult called for master element");

  // dual LM shape functions or not
  bool dual = false;
  if (lmtype==INPAR::MORTAR::shape_dual || lmtype==INPAR::MORTAR::shape_petrovgalerkin)
    dual = true;
  
  // get node number and node pointers
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");

  switch(Shape())
  {
  // 2D linear case (2noded line element)
  case DRT::Element::line2:
  {
    if (valdim!=2) dserror("ERROR: Inconsistency in EvluateShape");
    // check for boundary nodes
    MortarNode* mymrtrnode0 = static_cast<MortarNode*> (mynodes[0]);
    MortarNode* mymrtrnode1 = static_cast<MortarNode*> (mynodes[1]);
    if (!mymrtrnode0) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    if (!mymrtrnode1) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    bool isonbound0 = mymrtrnode0->IsOnBound();
    bool isonbound1 = mymrtrnode1->IsOnBound();

    // both nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1)
    {
      if (dual) ShapeFunctions(MortarElement::lindual1D,xi,val,deriv);
      else      ShapeFunctions(MortarElement::lin1D,xi,val,deriv);
    }

    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1)
    {
      if (dual) ShapeFunctions(MortarElement::lindual1D_edge0,xi,val,deriv);
      else      ShapeFunctions(MortarElement::lin1D_edge0,xi,val,deriv);
    }

    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1)
    {
      if (dual) ShapeFunctions(MortarElement::lindual1D_edge1,xi,val,deriv);
      else      ShapeFunctions(MortarElement::lin1D_edge1,xi,val,deriv);
    }

    // both nodes are on boundary: infeasible case
    else
      dserror("ERROR: EvaluateShapeLagMult: Element with 2 boundary nodes!");

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
    if (!mymrtrnode0) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    if (!mymrtrnode1) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    if (!mymrtrnode2) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    bool isonbound0 = mymrtrnode0->IsOnBound();
    bool isonbound1 = mymrtrnode1->IsOnBound();
    bool isonbound2 = mymrtrnode2->IsOnBound();

    // all 3 nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1 && !isonbound2)
    {
      if (dual) ShapeFunctions(MortarElement::quaddual1D,xi,val,deriv);
      else      ShapeFunctions(MortarElement::quad1D,xi,val,deriv);
    }

    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1 && !isonbound2)
    {
      if (dual) ShapeFunctions(MortarElement::quaddual1D_edge0,xi,val,deriv);
      else      ShapeFunctions(MortarElement::quad1D_edge0,xi,val,deriv);
    }

    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1 && !isonbound2)
    {
      if (dual) ShapeFunctions(MortarElement::quaddual1D_edge1,xi,val,deriv);
      else      ShapeFunctions(MortarElement::quad1D_edge1,xi,val,deriv);
    }

    // node 2 is on boundary: infeasible case
    else if (isonbound2)
      dserror("ERROR: EvaluateShapeLagMult: Middle boundary node");

    // nodes 0 and 1 are on boundary: infeasible case
    else
      dserror("ERROR: EvaluateShapeLagMult: Element with 2 boundary nodes");

    break;
  }

  // 3D cases
  case DRT::Element::tri3:
  case DRT::Element::quad4:
  case DRT::Element::tri6:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
      // dual Lagrange multipliers
      if (dual)
      {
        if (Shape()==tri3)       ShapeFunctions(MortarElement::lindual2D,xi,val,deriv);
        else if (Shape()==quad4) ShapeFunctions(MortarElement::bilindual2D,xi,val,deriv);
        else if (Shape()==tri6)  ShapeFunctions(MortarElement::quaddual2D,xi,val,deriv);
        else if (Shape()==quad8) ShapeFunctions(MortarElement::serendipitydual2D,xi,val,deriv);
        else /*Shape()==quad9*/  ShapeFunctions(MortarElement::biquaddual2D,xi,val,deriv);
      }
      
      // standard Lagrange multipliers
      else
      {
        if (Shape()==tri3)       ShapeFunctions(MortarElement::lin2D,xi,val,deriv);
        else if (Shape()==quad4) ShapeFunctions(MortarElement::bilin2D,xi,val,deriv);
        else if (Shape()==tri6)  ShapeFunctions(MortarElement::quad2D,xi,val,deriv);
        else if (Shape()==quad8) ShapeFunctions(MortarElement::serendipity2D,xi,val,deriv);
        else /*Shape()==quad9*/  ShapeFunctions(MortarElement::biquad2D,xi,val,deriv);
      }
      
    break;
  }

  // unknown case
  default:
    dserror("ERROR: EvaluateShapeLagMult called for unknown element type");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate Lagrange multiplier shape functions              popp 12/07|
 |  THIS IS A SPECIAL VERSION FOR 3D QUADRATIC MORTAR WITH LIN LM!      |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::EvaluateShapeLagMultLin(const INPAR::MORTAR::ShapeFcn& lmtype,
                                                    const double* xi, LINALG::SerialDenseVector& val,
                                                    LINALG::SerialDenseMatrix& deriv, const int& valdim)
{
  if (!xi)
    dserror("ERROR: EvaluateShapeLagMultLin called with xi=NULL");
  if (!IsSlave())
    dserror("ERROR: EvaluateShapeLagMultLin called for master element");

  // check for feasible element types (line3,tri6, quad8 or quad9)
  if (Shape()!=DRT::Element::line3 &&
      Shape()!=DRT::Element::tri6 &&
      Shape()!=DRT::Element::quad8 &&
      Shape()!=DRT::Element::quad9)
    dserror("ERROR: Linear LM interpolation only for quadratic finite elements");
  
  // dual shape functions or not
  bool dual = false;
  if (lmtype==INPAR::MORTAR::shape_dual || lmtype==INPAR::MORTAR::shape_petrovgalerkin)
    dual = true;
  
  // get node number and node pointers
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");

  // check for boundary nodes
  bool bound = false;
  for (int i=0;i<NumNode();++i)
  {
    MortarNode* mymrtrnode = static_cast<MortarNode*> (mynodes[i]);
    if (!mymrtrnode) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // all nodes are interior: use unmodified shape functions
  if (!bound)
  {
    dserror("ERROR: You should not be here...");
  }

  switch(Shape())
  {
    // 2D quadratic case (quadratic line)
    case DRT::Element::line3:
    {
      // the middle node is defined as slave boundary (=master)
      // dual Lagrange multipliers
      if (dual)
      {
        dserror("ERROR: Quad->Lin modification of dual LM shape functions not yet implemented");
        ShapeFunctions(MortarElement::quaddual1D_only_lin,xi,val,deriv);
      }

      // standard Lagrange multipliers
      else
      {
        ShapeFunctions(MortarElement::quad1D_only_lin,xi,val,deriv);
      }

      break;
    }

    // 3D quadratic cases (quadratic triangle, biquadratic and serendipity quad)
    case DRT::Element::tri6:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    {
      // the edge nodes are defined as slave boundary (=master)
      // dual Lagrange multipliers
      if (dual)
      {
        //dserror("ERROR: Quad->Lin modification of dual LM shape functions not yet implemented");
        if (Shape()==tri6)       ShapeFunctions(MortarElement::quaddual2D_only_lin,xi,val,deriv);
        else if (Shape()==quad8) ShapeFunctions(MortarElement::serendipitydual2D_only_lin,xi,val,deriv);
        else /*Shape()==quad9*/  ShapeFunctions(MortarElement::biquaddual2D_only_lin,xi,val,deriv);
      }

      // standard Lagrange multipliers
      else
      {
        if (Shape()==tri6)       ShapeFunctions(MortarElement::quad2D_only_lin,xi,val,deriv);
        else if (Shape()==quad8) ShapeFunctions(MortarElement::serendipity2D_only_lin,xi,val,deriv);
        else /*Shape()==quad9*/  ShapeFunctions(MortarElement::biquad2D_only_lin,xi,val,deriv);
      }
  
      break;
    }
      
    // unknown case
    default:
      dserror("ERROR: EvaluateShapeLagMult called for unknown element type");
  }

  return true;
}
/*----------------------------------------------------------------------*
 |  1D/2D shape function linearizations repository            popp 05/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::ShapeFunctionLinearizations(MORTAR::MortarElement::ShapeType shape,
                                   std::vector<std::vector<std::map<int,double> > >& derivdual)
{
  switch(shape)
  {
  // *********************************************************************
  // 1D dual quadratic shape functions (line3)
  // 2D dual bilinear shape functions (quad4)
  // 2D dual biquadratic shape functions (quad9)
  // (used for interpolation of Lagrange multiplier field)
  // (linearization necessary due to adaption for distorted elements !!!)
  // *********************************************************************
  case MORTAR::MortarElement::quaddual1D:
  case MORTAR::MortarElement::bilindual2D:
  case MORTAR::MortarElement::biquaddual2D:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    GetNodalCoords(coord);

    // prepare computation with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,2,true);
    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

    // two-dim arrays of maps for linearization of me/de
    std::vector<std::vector<std::map<int,double> > > derivme(nnodes,std::vector<std::map<int,double> >(nnodes));
    std::vector<std::vector<std::map<int,double> > > derivde(nnodes,std::vector<std::map<int,double> >(nnodes));

    // build me, de, derivme, derivde
    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
      EvaluateShape(gpc, val, deriv, nnodes);
      detg = Jacobian(gpc);

      // directional derivative of Jacobian
      std::map<int,double> testmap;
      typedef std::map<int,double>::const_iterator CI;
      DerivJacobian(gpc,testmap);

      // loop over all entries of me/de
      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          double facme = integrator.Weight(i)*val[j]*val[k];
          double facde = (j==k)*integrator.Weight(i)*val[j];

          me(j,k)+=facme*detg;
          de(j,k)+=facde*detg;

          // loop over all directional derivatives
          for (CI p=testmap.begin();p!=testmap.end();++p)
          {
            derivme[j][k][p->first] += facme*(p->second);
            derivde[j][k][p->first] += facde*(p->second);
          }
        }
    }

    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);

    // get solution matrix ae with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // build linearization of ae and store in derivdual
    // (this is done according to a quite complex formula, which
    // we get from the linearization of the biorthogonality condition:
    // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
    typedef std::map<int,double>::const_iterator CI;

    // loop over all entries of ae (index i,j)
    for (int i=0;i<nnodes;++i)
    {
      for (int j=0;j<nnodes;++j)
      {
        // compute Lin(Ae) according to formula above
        for (int l=0;l<nnodes;++l) // loop over sum l
        {
          // part1: Lin(De)*Inv(Me)
          for (CI p=derivde[i][l].begin();p!=derivde[i][l].end();++p)
            derivdual[i][j][p->first] += me(l,j)*(p->second);

          // part2: Ae*Lin(Me)*Inv(Me)
          for (int k=0;k<nnodes;++k) // loop over sum k
          {
            for (CI p=derivme[k][l].begin();p!=derivme[k][l].end();++p)
              derivdual[i][j][p->first] -= ae(i,k)*me(l,j)*(p->second);
          }
        }
      }
    }

    // std::cout linearization of Ae
    //std::cout << "Analytical A-derivative of Element: " << Id() << endl;
    //for (int i=0;i<nnodes;++i)
    //  for (int j=0;j<nnodes;++j)
    //    for (CI p=derivdual[i][j].begin();p!=derivdual[i][j].end();++p)
    //      std::cout << "A" << i << j << " " << p->first << " " << p->second << endl;

    /*
#ifdef DEBUG
    // *******************************************************************
    // FINITE DIFFERENCE check of Lin(Ae)
    // *******************************************************************

    std::cout << "FD Check for A-derivative of Element: " << Id() << endl;
    Epetra_SerialDenseMatrix aeref(ae);
    double delta = 1e-8;
    int thedim=3;
    if (shape==MORTAR::MortarElement::quaddual1D) thedim=2;

    for (int dim=0;dim<thedim;++dim)
    {
      for (int node=0;node<nnodes;++node)
      {
        // apply FD
        DRT::Node** mynodes = Nodes();
        CoNode* mycnode = static_cast<CoNode*> (mynodes[node]);
        mycnode->xspatial()[dim] += delta;

        LINALG::SerialDenseVector val1(nnodes);
        LINALG::SerialDenseMatrix deriv1(nnodes,2,true);
        LINALG::SerialDenseMatrix me1(nnodes,nnodes,true);
        LINALG::SerialDenseMatrix de1(nnodes,nnodes,true);

        // build me, de
        for (int i=0;i<integrator.nGP();++i)
        {
          double gpc1[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
          EvaluateShape(gpc1, val1, deriv1, nnodes);
          detg = Jacobian(gpc1);

          for (int j=0;j<nnodes;++j)
            for (int k=0;k<nnodes;++k)
            {
              double facme1 = integrator.Weight(i)*val1[j]*val1[k];
              double facde1 = (j==k)*integrator.Weight(i)*val1[j];

              me1(j,k)+=facme1*detg;
              de1(j,k)+=facde1*detg;
            }
        }

        // invert bi-ortho matrix me
        LINALG::SymmetricInverse(me1,nnodes);

        // get solution matrix ae with dual parameters
        LINALG::SerialDenseMatrix ae1(nnodes,nnodes);
        ae1.Multiply('N','N',1.0,de1,me1,0.0);
        int col= mycnode->Dofs()[dim];

        std::cout << "A-Derivative: " << col << endl;

        // FD solution
        for (int i=0;i<nnodes;++i)
          for (int j=0;j<nnodes;++j)
          {
            double val = (ae1(i,j)-aeref(i,j))/delta;
            std::cout << "A" << i << j << " " << val << endl;
          }

        // undo FD
        mycnode->xspatial()[dim] -= delta;
      }
    }
    // *******************************************************************
#endif // #ifdef DEBUG
    */

    break;
  }
  // *********************************************************************
  // 2D dual quadratic shape functions (tri6)
  // 2D dual serendipity shape functions (quad8)
  // (used for interpolation of Lagrange multiplier field)
  // (linearization necessary due to adaption for distorted elements !!!)
  // (including modification of displacement shape functions)
  // *********************************************************************
  case MORTAR::MortarElement::quaddual2D:
  case MORTAR::MortarElement::serendipitydual2D:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    GetNodalCoords(coord);

    // prepare computation with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());
    LINALG::SerialDenseVector val(nnodes);
    LINALG::SerialDenseMatrix deriv(nnodes,2,true);
    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

    // two-dim arrays of maps for linearization of me/de
    std::vector<std::vector<std::map<int,double> > > derivme(nnodes,std::vector<std::map<int,double> >(nnodes));
    std::vector<std::vector<std::map<int,double> > > derivde(nnodes,std::vector<std::map<int,double> >(nnodes));

    // build me, de, derivme, derivde
    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
      EvaluateShape(gpc, val, deriv, nnodes, true);
      detg = Jacobian(gpc);

      // directional derivative of Jacobian
      std::map<int,double> testmap;
      typedef std::map<int,double>::const_iterator CI;
      DerivJacobian(gpc,testmap);

      // loop over all entries of me/de
      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          double facme = integrator.Weight(i)*val[j]*val[k];
          double facde = (j==k)*integrator.Weight(i)*val[j];

          me(j,k)+=facme*detg;
          de(j,k)+=facde*detg;

          // loop over all directional derivatives
          for (CI p=testmap.begin();p!=testmap.end();++p)
          {
            derivme[j][k][p->first] += facme*(p->second);
            derivde[j][k][p->first] += facde*(p->second);
          }
        }
    }

    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);

    // get solution matrix ae with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // build linearization of ae and store in derivdual
    // (this is done according to a quite complex formula, which
    // we get from the linearization of the biorthogonality condition:
    // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
    typedef std::map<int,double>::const_iterator CI;

    // loop over all entries of ae (index i,j)
    for (int i=0;i<nnodes;++i)
    {
      for (int j=0;j<nnodes;++j)
      {
        // compute Lin(Ae) according to formula above
        for (int l=0;l<nnodes;++l) // loop over sum l
        {
          // part1: Lin(De)*Inv(Me)
          for (CI p=derivde[i][l].begin();p!=derivde[i][l].end();++p)
            derivdual[i][j][p->first] += me(l,j)*(p->second);

          // part2: Ae*Lin(Me)*Inv(Me)
          for (int k=0;k<nnodes;++k) // loop over sum k
          {
            for (CI p=derivme[k][l].begin();p!=derivme[k][l].end();++p)
              derivdual[i][j][p->first] -= ae(i,k)*me(l,j)*(p->second);
          }
        }
      }
    }

    // std::cout linearization of Ae
    //std::cout << "Analytical A-derivative of Element: " << Id() << endl;
    //for (int i=0;i<nnodes;++i)
    //  for (int j=0;j<nnodes;++j)
    //    for (CI p=derivdual[i][j].begin();p!=derivdual[i][j].end();++p)
    //      std::cout << "A" << i << j << " " << p->first << " " << p->second << endl;
    /*
#ifdef DEBUG
    // *******************************************************************
    // FINITE DIFFERENCE check of Lin(Ae)
    // *******************************************************************

    std::cout << "FD Check for A-derivative of Element: " << Id() << endl;
    Epetra_SerialDenseMatrix aeref(ae);
    double delta = 1e-8;
    int thedim=3;

    for (int dim=0;dim<thedim;++dim)
    {
      for (int node=0;node<nnodes;++node)
      {
        // apply FD
        DRT::Node** mynodes = Nodes();
        MortarNode* mycnode = static_cast<MortarNode*> (mynodes[node]);
        mycnode->xspatial()[dim] += delta;

        LINALG::SerialDenseVector val1(nnodes);
        LINALG::SerialDenseMatrix deriv1(nnodes,2,true);
        LINALG::SerialDenseMatrix me1(nnodes,nnodes,true);
        LINALG::SerialDenseMatrix de1(nnodes,nnodes,true);

        // build me, de
        for (int i=0;i<integrator.nGP();++i)
        {
          double gpc1[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
          EvaluateShape(gpc1, val1, deriv1, nnodes, true);
          detg = Jacobian(gpc1);

          for (int j=0;j<nnodes;++j)
            for (int k=0;k<nnodes;++k)
            {
              double facme1 = integrator.Weight(i)*val1[j]*val1[k];
              double facde1 = (j==k)*integrator.Weight(i)*val1[j];

              me1(j,k)+=facme1*detg;
              de1(j,k)+=facde1*detg;
            }
        }

        // invert bi-ortho matrix me
        LINALG::SymmetricInverse(me1,nnodes);

        // get solution matrix ae with dual parameters
        LINALG::SerialDenseMatrix ae1(nnodes,nnodes);
        ae1.Multiply('N','N',1.0,de1,me1,0.0);
        int col= mycnode->Dofs()[dim];

        std::cout << "A-Derivative: " << col << endl;

        // FD solution
        for (int i=0;i<nnodes;++i)
          for (int j=0;j<nnodes;++j)
          {
            double val = (ae1(i,j)-aeref(i,j))/delta;
            std::cout << "A" << i << j << " " << val << endl;
          }

        // undo FD
        mycnode->xspatial()[dim] -= delta;
      }
    }
    // *******************************************************************
#endif // #ifdef DEBUG
    */
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (linearization necessary due to adaption for distorted elements !!!)
  // *********************************************************************
  case MORTAR::MortarElement::quaddual1D_edge0:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    GetNodalCoords(coord);

    // empty shape function vals + derivs
    LINALG::SerialDenseVector valquad(nnodes);
    LINALG::SerialDenseMatrix derivquad(nnodes,1);
    LINALG::SerialDenseVector vallin(nnodes-1);
    LINALG::SerialDenseMatrix derivlin(nnodes-1,1);

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());

    LINALG::SerialDenseMatrix me(nnodes-1,nnodes-1,true);
    LINALG::SerialDenseMatrix de(nnodes-1,nnodes-1,true);

    // two-dim arrays of maps for linearization of me/de
    std::vector<std::vector<std::map<int,double> > > derivme(nnodes,std::vector<std::map<int,double> >(nnodes));
    std::vector<std::vector<std::map<int,double> > > derivde(nnodes,std::vector<std::map<int,double> >(nnodes));

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), 0.0};
      ShapeFunctions(MORTAR::MortarElement::quad1D,gpc,valquad,derivquad);
      ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge0,gpc,vallin,derivlin);
      detg = Jacobian(gpc);

      // directional derivative of Jacobian
      std::map<int,double> testmap;
      typedef std::map<int,double>::const_iterator CI;
      DerivJacobian(gpc,testmap);

      // loop over all entries of me/de
      for (int j=1;j<nnodes;++j)
        for (int k=1;k<nnodes;++k)
        {
          double facme = integrator.Weight(i)*vallin[j-1]*valquad[k];
          double facde = (j==k)*integrator.Weight(i)*valquad[k];

          me(j-1,k-1)+=facme*detg;
          de(j-1,k-1)+=facde*detg;

          // loop over all directional derivatives
          for (CI p=testmap.begin();p!=testmap.end();++p)
          {
            derivme[j-1][k-1][p->first] += facme*(p->second);
            derivde[j-1][k-1][p->first] += facde*(p->second);
          }
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

    // build linearization of ae and store in derivdual
    // (this is done according to a quite complex formula, which
    // we get from the linearization of the biorthogonality condition:
    // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
    typedef std::map<int,double>::const_iterator CI;

    // loop over all entries of ae (index i,j)
    for (int i=1;i<nnodes;++i)
    {
      for (int j=1;j<nnodes;++j)
      {
        // compute Lin(Ae) according to formula above
        for (int l=1;l<nnodes;++l) // loop over sum l
        {
          // part1: Lin(De)*Inv(Me)
          for (CI p=derivde[i-1][l-1].begin();p!=derivde[i-1][l-1].end();++p)
            derivdual[i][j][p->first] += me(l-1,j-1)*(p->second);

          // part2: Ae*Lin(Me)*Inv(Me)
          for (int k=1;k<nnodes;++k) // loop over sum k
          {
            for (CI p=derivme[k-1][l-1].begin();p!=derivme[k-1][l-1].end();++p)
              derivdual[i][j][p->first] -= ae(i-1,k-1)*me(l-1,j-1)*(p->second);
          }
        }
      }
    }

    // std::cout linearization of Ae
    //std::cout << "Analytical A-derivative of Element: " << Id() << endl;
    //for (int i=1;i<nnodes;++i)
    //  for (int j=1;j<nnodes;++j)
    //    for (CI p=derivdual[i][j].begin();p!=derivdual[i][j].end();++p)
    //      std::cout << "A" << i << j << " " << p->first << " " << p->second << endl;
    /*
#ifdef DEBUG
    // *******************************************************************
    // FINITE DIFFERENCE check of Lin(Ae)
    // *******************************************************************

    std::cout << "FD Check for A-derivative of Element: " << Id() << endl;
    LINALG::SerialDenseMatrix aeref(ae);
    double delta = 1e-8;

    for (int dim=0;dim<2;++dim)
    {
      for (int node=0;node<nnodes;++node)
      {
        // apply FD
        coord(dim,node)+=delta;

        // empty shape function vals + derivs
        LINALG::SerialDenseVector valquad1(nnodes);
        LINALG::SerialDenseMatrix derivquad1(nnodes,1);
        LINALG::SerialDenseVector vallin1(nnodes-1);
        LINALG::SerialDenseMatrix derivlin1(nnodes-1,1);
        //LINALG::SerialDenseVector valtemp1(nnodes);
        //LINALG::SerialDenseMatrix derivtemp1(nnodes,1);
        LINALG::SerialDenseMatrix me1(nnodes-1,nnodes-1,true);
        LINALG::SerialDenseMatrix de1(nnodes-1,nnodes-1,true);

        // build me, de
        for (int i=0;i<integrator.nGP();++i)
        {
          double gpc1[2] = {integrator.Coordinate(i), 0.0};
          ShapeFunctions(MORTAR::MortarElement::quad1D,gpc1,valquad1,derivquad1);
          ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge0,gpc1,vallin1,derivlin1);
          detg = Jacobian(valquad1,derivquad1,coord);

          for (int j=1;j<nnodes;++j)
            for (int k=1;k<nnodes;++k)
            {
              double facme1 = integrator.Weight(i)*vallin1[j-1]*valquad1[k];
              double facde1 = (j==k)*integrator.Weight(i)*valquad1[k];

              me1(j-1,k-1)+=facme1*detg;
              de1(j-1,k-1)+=facde1*detg;
            }
        }

        // invert bi-ortho matrix me1
        double detme1 = me1(0,0)*me1(1,1)-me1(0,1)*me1(1,0);
        LINALG::SerialDenseMatrix meold(nnodes-1,nnodes-1);
        meold=me1;
        me1(0,0) =  1/detme1*meold(1,1);
        me1(0,1) = -1/detme1*meold(0,1);
        me1(1,0) = -1/detme1*meold(1,0);
        me1(1,1) =  1/detme1*meold(0,0);

        // get solution matrix ae with dual parameters
        LINALG::SerialDenseMatrix ae1(nnodes-1,nnodes-1);
        ae1.Multiply('N','N',1.0,de1,me1,0.0);

        DRT::Node** mynodes = Nodes();
        CoNode* mycnode = static_cast<CoNode*> (mynodes[node]);
        int col= mycnode->Dofs()[dim];

        std::cout << "A-Derivative: " << col << endl;

        // FD solution
        for (int i=1;i<nnodes;++i)
          for (int j=1;j<nnodes;++j)
          {
            double val = (ae1(i-1,j-1)-aeref(i-1,j-1))/delta;
            std::cout << "A" << i << j << " " << val << endl;
          }

        // undo FD
        coord(dim,node)-=delta;
      }
    }
    // *******************************************************************
#endif // #ifdef DEBUG
    */
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (linearization necessary due to adaption for distorted elements !!!)
  // *********************************************************************
  case MORTAR::MortarElement::quaddual1D_edge1:
  {
    // establish fundamental data
    double detg = 0.0;
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    GetNodalCoords(coord);

    // empty shape function vals + derivs
    LINALG::SerialDenseVector valquad(nnodes);
    LINALG::SerialDenseMatrix derivquad(nnodes,1);
    LINALG::SerialDenseVector vallin(nnodes-1);
    LINALG::SerialDenseMatrix derivlin(nnodes-1,1);

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::ElementIntegrator integrator(Shape());

    LINALG::SerialDenseMatrix me(nnodes-1,nnodes-1,true);
    LINALG::SerialDenseMatrix de(nnodes-1,nnodes-1,true);

    // two-dim arrays of maps for linearization of me/de
    std::vector<std::vector<std::map<int,double> > > derivme(nnodes,std::vector<std::map<int,double> >(nnodes));
    std::vector<std::vector<std::map<int,double> > > derivde(nnodes,std::vector<std::map<int,double> >(nnodes));

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), 0.0};
      ShapeFunctions(MORTAR::MortarElement::quad1D,gpc,valquad,derivquad);
      ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge1,gpc,vallin,derivlin);
      detg = Jacobian(gpc);

      // directional derivative of Jacobian
      std::map<int,double> testmap;
      typedef std::map<int,double>::const_iterator CI;
      DerivJacobian(gpc,testmap);

      // loop over all entries of me/de
      for (int j=0;j<nnodes-1;++j)
        for (int k=0;k<nnodes-1;++k)
        {
          double facme = integrator.Weight(i)*vallin[j]*valquad[2*k];
          double facde = (j==k)*integrator.Weight(i)*valquad[2*k];

          me(j,k)+=facme*detg;
          de(j,k)+=facde*detg;

          // loop over all directional derivatives
          for (CI p=testmap.begin();p!=testmap.end();++p)
          {
            derivme[j][k][p->first] += facme*(p->second);
            derivde[j][k][p->first] += facde*(p->second);
          }
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

    // build linearization of ae and store in derivdual
    // (this is done according to a quite complex formula, which
    // we get from the linearization of the biorthogonality condition:
    // Lin (Me * Ae = De) -> Lin(Ae)=Lin(De)*Inv(Me)-Ae*Lin(Me)*Inv(Me) )
    typedef std::map<int,double>::const_iterator CI;

    // loop over all entries of ae (index i,j)
    for (int i=0;i<nnodes-1;++i)
    {
      for (int j=0;j<nnodes-1;++j)
      {
        // compute Lin(Ae) according to formula above
        for (int l=0;l<nnodes-1;++l) // loop over sum l
        {
          // part1: Lin(De)*Inv(Me)
          for (CI p=derivde[i][l].begin();p!=derivde[i][l].end();++p)
            derivdual[i][j][p->first] += me(l,j)*(p->second);

          // part2: Ae*Lin(Me)*Inv(Me)
          for (int k=0;k<nnodes-1;++k) // loop over sum k
          {
            for (CI p=derivme[k][l].begin();p!=derivme[k][l].end();++p)
              derivdual[i][j][p->first] -= ae(i,k)*me(l,j)*(p->second);
          }
        }
      }
    }

    // std::cout linearization of Ae
    //std::cout << "Analytical A-derivative of Element: " << Id() << endl;
    //for (int i=0;i<nnodes-1;++i)
    //  for (int j=0;j<nnodes-1;++j)
    //    for (CI p=derivdual[i][j].begin();p!=derivdual[i][j].end();++p)
    //      std::cout << "A" << i << j << " " << p->first << " " << p->second << endl;
    /*
#ifdef DEBUG
    // *******************************************************************
    // FINITE DIFFERENCE check of Lin(Ae)
    // *******************************************************************

    std::cout << "FD Check for A-derivative of Element: " << Id() << endl;
    LINALG::SerialDenseMatrix aeref(ae);
    double delta = 1e-8;

    for (int dim=0;dim<2;++dim)
    {
      for (int node=0;node<nnodes;++node)
      {
        // apply FD
        coord(dim,node)+=delta;

        // empty shape function vals + derivs
        LINALG::SerialDenseVector valquad1(nnodes);
        LINALG::SerialDenseMatrix derivquad1(nnodes,1);
        LINALG::SerialDenseVector vallin1(nnodes-1);
        LINALG::SerialDenseMatrix derivlin1(nnodes-1,1);
        //LINALG::SerialDenseVector valtemp1(nnodes);
        //LINALG::SerialDenseMatrix derivtemp1(nnodes,1);
        LINALG::SerialDenseMatrix me1(nnodes-1,nnodes-1,true);
        LINALG::SerialDenseMatrix de1(nnodes-1,nnodes-1,true);

        // build me, de
        for (int i=0;i<integrator.nGP();++i)
        {
          double gpc1[2] = {integrator.Coordinate(i), 0.0};
          ShapeFunctions(MORTAR::MortarElement::quad1D,gpc1,valquad1,derivquad1);
          ShapeFunctions(MORTAR::MortarElement::dual1D_base_for_edge1,gpc1,vallin1,derivlin1);
          detg = Jacobian(valquad1,derivquad1,coord);

          for (int j=0;j<nnodes-1;++j)
            for (int k=0;k<nnodes-1;++k)
            {
              double facme1 = integrator.Weight(i)*vallin1[j]*valquad1[2*k];
              double facde1 = (j==k)*integrator.Weight(i)*valquad1[2*k];

              me1(j,k)+=facme1*detg;
              de1(j,k)+=facde1*detg;
            }
        }

        // invert bi-ortho matrix me1
        double detme1 = me1(0,0)*me1(1,1)-me1(0,1)*me1(1,0);
        LINALG::SerialDenseMatrix meold(nnodes-1,nnodes-1);
        meold=me1;
        me1(0,0) =  1/detme1*meold(1,1);
        me1(0,1) = -1/detme1*meold(0,1);
        me1(1,0) = -1/detme1*meold(1,0);
        me1(1,1) =  1/detme1*meold(0,0);

        // get solution matrix ae with dual parameters
        LINALG::SerialDenseMatrix ae1(nnodes-1,nnodes-1);
        ae1.Multiply('N','N',1.0,de1,me1,0.0);

        DRT::Node** mynodes = Nodes();
        CoNode* mycnode = static_cast<CoNode*> (mynodes[node]);
        int col= mycnode->Dofs()[dim];

        std::cout << "A-Derivative: " << col << endl;

        // FD solution
        for (int i=0;i<nnodes-1;++i)
          for (int j=0;j<nnodes-1;++j)
          {
            double val = (ae1(i,j)-aeref(i,j))/delta;
            std::cout << "A" << i << j << " " << val << endl;
          }

        // undo FD
        coord(dim,node)-=delta;
      }
    }
    // *******************************************************************
#endif // #ifdef DEBUG
    */
    break;
  }
  // *********************************************************************
  // Unknown shape function type
  // *********************************************************************
  default:
    dserror("ERROR: Unknown shape function type identifier");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate 2nd derivative of shape functions                popp 05/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::Evaluate2ndDerivShape(const double* xi,
                                              LINALG::SerialDenseMatrix& secderiv,
                                              const int& valdim)
{
  if (!xi)
    dserror("ERROR: Evaluate2ndDerivShape called with xi=NULL");

  //**********************************************************************
  // IMPORTANT NOTE: In 3D the ordering of the 2nd derivatives is:
  // 1) dxi,dxi 2) deta,deta 3) dxi,deta
  //**********************************************************************

  switch(Shape())
  {
  // 2D linear case (2noded line element)
  case DRT::Element::line2:
  {
    secderiv(0,0) = 0;
    secderiv(1,0) = 0;
    break;
  }
  // 2D quadratic case (3noded line element)
  case DRT::Element::line3:
  {
    secderiv(0,0) =  1;
    secderiv(1,0) =  1;
    secderiv(2,0) = -2;
    break;
  }
  // 3D linear case (3noded triangular element)
  case DRT::Element::tri3:
  {
    secderiv(0,0) =  0; secderiv(0,1) =  0; secderiv(0,2) =  0;
    secderiv(1,0) =  0; secderiv(1,1) =  0; secderiv(1,2) =  0;
    secderiv(2,0) =  0; secderiv(2,1) =  0; secderiv(2,2) =  0;
    break;
  }
  // 3D bilinear case (4noded quadrilateral element)
  case DRT::Element::quad4:
  {
    secderiv(0,0) =  0; secderiv(0,1) =  0; secderiv(0,2) =  0.25;
    secderiv(1,0) =  0; secderiv(1,1) =  0; secderiv(1,2) = -0.25;
    secderiv(2,0) =  0; secderiv(2,1) =  0; secderiv(2,2) =  0.25;
    secderiv(3,0) =  0; secderiv(3,1) =  0; secderiv(3,2) = -0.25;
    break;
  }
  // 3D quadratic case (6noded triangular element)
  case DRT::Element::tri6:
  {
    secderiv(0,0) =  4.0; secderiv(0,1) =  4.0; secderiv(0,2) =  4.0;
    secderiv(1,0) =  4.0; secderiv(1,1) =    0; secderiv(1,2) =    0;
    secderiv(2,0) =    0; secderiv(2,1) =  4.0; secderiv(2,2) =    0;
    secderiv(3,0) = -8.0; secderiv(3,1) =    0; secderiv(3,2) = -4.0;
    secderiv(4,0) =    0; secderiv(4,1) =    0; secderiv(4,2) =  4.0;
    secderiv(5,0) =    0; secderiv(5,1) = -8.0; secderiv(5,2) = -4.0;
    break;
  }
  // 3D serendipity case (8noded quadrilateral element)
  case DRT::Element::quad8:
  {
    const double r=xi[0];
    const double s=xi[1];
    const double rp=1.0+r;
    const double rm=1.0-r;
    const double sp=1.0+s;
    const double sm=1.0-s;

    secderiv(0,0) = 0.5*sm; secderiv(0,1) = 0.5*rm; secderiv(0,2) =  -0.25*(2*r+2*s-1.0);
    secderiv(1,0) = 0.5*sm; secderiv(1,1) = 0.5*rp; secderiv(1,2) =  0.25*(-2*r+2*s-1.0);
    secderiv(2,0) = 0.5*sp; secderiv(2,1) = 0.5*rp; secderiv(2,2) =   0.25*(2*r+2*s+1.0);
    secderiv(3,0) = 0.5*sp; secderiv(3,1) = 0.5*rm; secderiv(3,2) = -0.25*(-2*r+2*s+1.0);
    secderiv(4,0) =    -sm; secderiv(4,1) =      0; secderiv(4,2) =                    r;
    secderiv(5,0) =      0; secderiv(5,1) =    -rp; secderiv(5,2) =                   -s;
    secderiv(6,0) =    -sp; secderiv(6,1) =      0; secderiv(6,2) =                   -r;
    secderiv(7,0) =      0; secderiv(7,1) =    -rm; secderiv(7,2) =                    s;
    break;
  }
  // 3D biquadratic case (9noded quadrilateral element)
  case DRT::Element::quad9:
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
    const double rhp=r+0.5;
    const double rhm=r-0.5;
    const double shp=s+0.5;
    const double shm=s-0.5;

    secderiv(0,0) =     -sh*sm; secderiv(0,1) =     -rh*rm; secderiv(0,2) =     shm*rhm;
    secderiv(1,0) =     -sh*sm; secderiv(1,1) =      rh*rp; secderiv(1,2) =     shm*rhp;
    secderiv(2,0) =      sh*sp; secderiv(2,1) =      rh*rp; secderiv(2,2) =     shp*rhp;
    secderiv(3,0) =      sh*sp; secderiv(3,1) =     -rh*rm; secderiv(3,2) =     shp*rhm;
    secderiv(4,0) =  2.0*sh*sm; secderiv(4,1) =         r2; secderiv(4,2) =  -2.0*r*shm;
    secderiv(5,0) =         s2; secderiv(5,1) = -2.0*rh*rp; secderiv(5,2) =  -2.0*s*rhp;
    secderiv(6,0) = -2.0*sh*sp; secderiv(6,1) =         r2; secderiv(6,2) =  -2.0*r*shp;
    secderiv(7,0) =         s2; secderiv(7,1) =  2.0*rh*rm; secderiv(7,2) =  -2.0*s*rhm;
    secderiv(8,0) =    -2.0*s2; secderiv(8,1) =    -2.0*r2; secderiv(8,2) = 2.0*s*2.0*r;
    break;
  }
  // unknown case
  default:
    dserror("ERROR: Evaluate2ndDerivShape called for unknown element type");
    break;
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Compute directional derivative of dual shape functions    popp 05/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::DerivShapeDual(std::vector<std::vector<std::map<int,double> > >& derivdual)
{
  if (!IsSlave()) dserror("ERROR: DerivShapeDual called for master element");

  // get node number and node pointers
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: DerivShapeDual: Null pointer!");

  switch(Shape())
  {
  // 2D linear case (2noded line element)
  // 3D linear case (3noded triangular element)
  case DRT::Element::line2:
  case DRT::Element::tri3:
  {
    dserror("ERROR: DerivShapeDual called for line2 or tri3!");
    break;
  }

  // 2D quadratic case (3noded line element)
  case DRT::Element::line3:
  {
    // check for boundary nodes
    MortarNode* mycnode0 = static_cast<MortarNode*> (mynodes[0]);
    MortarNode* mycnode1 = static_cast<MortarNode*> (mynodes[1]);
    MortarNode* mycnode2 = static_cast<MortarNode*> (mynodes[2]);
    if (!mycnode0) dserror("ERROR: DerivShapeDual: Null pointer!");
    if (!mycnode1) dserror("ERROR: DerivShapeDual: Null pointer!");
    if (!mycnode2) dserror("ERROR: DerivShapeDual: Null pointer!");
    bool isonbound0 = mycnode0->IsOnBound();
    bool isonbound1 = mycnode1->IsOnBound();
    bool isonbound2 = mycnode2->IsOnBound();

    // all 3 nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1 && !isonbound2)
      ShapeFunctionLinearizations(MORTAR::MortarElement::quaddual1D,derivdual);

    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1 && !isonbound2)
      ShapeFunctionLinearizations(MORTAR::MortarElement::quaddual1D_edge0,derivdual);

    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1 && !isonbound2)
      ShapeFunctionLinearizations(MORTAR::MortarElement::quaddual1D_edge1,derivdual);

    // node 2 is on boundary: infeasible case
    else if (isonbound2)
      dserror("ERROR: DerivShapeDual: Middle boundary node");

    // nodes 0 and 1 are on boundary: infeasible case
    else
      dserror("ERROR: DerivShapeDual: Element with 2 boundary nodes");

    break;
  }

  // all other 3D cases
  case DRT::Element::quad4:
  case DRT::Element::tri6:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    // check for boundary nodes
    bool bound = false;
    for (int i=0;i<NumNode();++i)
    {
      MortarNode* mycnode = static_cast<MortarNode*> (mynodes[i]);
      if (!mycnode) dserror("ERROR: EvaluateShapeLagMult: Null pointer!");
      bound += mycnode->IsOnBound();
    }

    // all nodes are interior: use unmodified dual shape functions
    if (!bound)
    {
      if (Shape()==quad4)      ShapeFunctionLinearizations(MORTAR::MortarElement::bilindual2D,derivdual);
      else if (Shape()==tri6)  ShapeFunctionLinearizations(MORTAR::MortarElement::quaddual2D,derivdual);
      else if (Shape()==quad8) ShapeFunctionLinearizations(MORTAR::MortarElement::serendipitydual2D,derivdual);
      else /*Shape()==quad9*/  ShapeFunctionLinearizations(MORTAR::MortarElement::biquaddual2D,derivdual);
    }

    // some nodes are on slave boundary
    else
      dserror("ERROR: DerivShapeDual: boundary mod. not yet impl. for 3D contact!");

    break;
  }

  // unknown case
  default:
    dserror("ERROR: DerivShapeDual called for unknown element type");
  }

  return true;
}

