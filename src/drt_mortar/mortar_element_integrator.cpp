/*!----------------------------------------------------------------------
\file mortar_element_integrator.cpp
\brief A class to perform Gaussian integration on a mortar element

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
#include "../drt_fem_general/drt_utils_integration.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 08/08|
 *----------------------------------------------------------------------*/
MORTAR::ElementIntegrator::ElementIntegrator(DRT::Element::DiscretizationType eletype)
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
  case DRT::Element::nurbs2:
  case DRT::Element::nurbs3:
  {
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
  case DRT::Element::nurbs4:
  case DRT::Element::nurbs8:
  case DRT::Element::nurbs9:
  {
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
    dserror("ERROR: ElementIntegrator: This contact element type is not implemented!");
  } // switch(eletype)

  return;
}

