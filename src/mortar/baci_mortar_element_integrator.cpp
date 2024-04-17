/*-----------------------------------------------------------------------*/
/*! \file
\level 2


\brief A class to perform Gaussian integration on a mortar element
*/
/*---------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 08/08|
 *----------------------------------------------------------------------*/
MORTAR::ElementIntegrator::ElementIntegrator(CORE::FE::CellType eletype)
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
  Teuchos::RCP<CORE::FE::IntegrationPoints2D> rule2d;

  switch (eletype)
  {
    case CORE::FE::CellType::line2:
    case CORE::FE::CellType::line3:
    case CORE::FE::CellType::nurbs2:
    case CORE::FE::CellType::nurbs3:
    {
      const CORE::FE::IntegrationPoints1D intpoints(CORE::FE::GaussRule1D::line_5point);
      ngp_ = intpoints.nquad;
      coords_.reshape(nGP(), 2);
      weights_.resize(nGP());
      for (int i = 0; i < nGP(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = 0.0;
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case CORE::FE::CellType::tri3:
      rule2d = Teuchos::rcp(new CORE::FE::IntegrationPoints2D(CORE::FE::GaussRule2D::tri_7point));
      break;
    case CORE::FE::CellType::tri6:
      rule2d = Teuchos::rcp(new CORE::FE::IntegrationPoints2D(CORE::FE::GaussRule2D::tri_16point));
      break;
    case CORE::FE::CellType::quad4:
      rule2d = Teuchos::rcp(new CORE::FE::IntegrationPoints2D(CORE::FE::GaussRule2D::quad_9point));
      break;
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
      rule2d = Teuchos::rcp(new CORE::FE::IntegrationPoints2D(CORE::FE::GaussRule2D::quad_25point));
      break;
    default:
      dserror("ElementIntegrator: This contact element type is not implemented!");
  }  // switch(eletype)

  // save Gauss points for all 2D rules
  if (rule2d != Teuchos::null)
  {
    ngp_ = rule2d->nquad;
    coords_.reshape(nGP(), 2);
    weights_.resize(nGP());
    for (int i = 0; i < nGP(); ++i)
    {
      coords_(i, 0) = rule2d->qxg[i][0];
      coords_(i, 1) = rule2d->qxg[i][1];
      weights_[i] = rule2d->qwgt[i];
    }
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
