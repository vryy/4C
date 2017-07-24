/*!----------------------------------------------------------------------
\file mortar_element_integrator.cpp

\level 2

\maintainer Alexander Seitz

\brief A class to perform Gaussian integration on a mortar element

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
  // 5  points: for integrals on 1D lines                 (1,2,3,4,5)
  // 7  points: for integrals on 2D 1st triangles         (1,3,6,7,12,16,37)
  // 16 points: for integrals on 2D 1st triangles         (1,3,6,7,12,16,37)
  // 9  points: for integrals on 2D 1st order quadrilaterals
  // 25 points: for integrals on 2D 2nd order quadrilaterals
  //**********************************************************************
  Teuchos::RCP<DRT::UTILS::IntegrationPoints2D> rule2d;

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
    rule2d = Teuchos::rcp(new DRT::UTILS::IntegrationPoints2D(DRT::UTILS::intrule_tri_7point));
    break;
  case DRT::Element::tri6:
    rule2d = Teuchos::rcp(new DRT::UTILS::IntegrationPoints2D(DRT::UTILS::intrule_tri_16point));
    break;
  case DRT::Element::quad4:
    rule2d = Teuchos::rcp(new DRT::UTILS::IntegrationPoints2D(DRT::UTILS::intrule_quad_9point));
    break;
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  case DRT::Element::nurbs4:
  case DRT::Element::nurbs9:
    rule2d = Teuchos::rcp(new DRT::UTILS::IntegrationPoints2D(DRT::UTILS::intrule_quad_25point));
    break;
  default:
    dserror("ERROR: ElementIntegrator: This contact element type is not implemented!");
  } // switch(eletype)

  // save Gauss points for all 2D rules
  if (rule2d!=Teuchos::null)
  {
    ngp_ = rule2d->nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=rule2d->qxg[i][0];
      coords_(i,1)=rule2d->qxg[i][1];
      weights_[i]=rule2d->qwgt[i];
    }
  }

  return;
}

