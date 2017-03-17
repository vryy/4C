/*---------------------------------------------------------------------*/
/*!
\file contact_integrator_utils.cpp

\brief Utility methods for the contact integration.

\level 2

\maintainer Michael Hiermeier

\date Mar 8, 2017

*/
/*---------------------------------------------------------------------*/

#include "contact_integrator_utils.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_contact/contact_integrator.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::INTEGRATOR::FindFeasibleMasterElement3D(
    MORTAR::MortarElement& sele,
    const std::vector<MORTAR::MortarElement*> & meles,
    bool boundary_ele,
    CoIntegrator& wrapper,
    UniqueProjInfoPair& projInfo )
{
  const unsigned msize = meles.size();
  const unsigned numGP = wrapper.nGP();

  projInfo.clear();
  if ( msize > projInfo.capacity() )
    projInfo.resize( msize );

  for ( unsigned gp=0; gp<numGP; ++gp )
  {
    // get Gauss point in slave element coordinates
    const double eta[2] = {wrapper.Coordinate(gp,0), wrapper.Coordinate(gp,1)};
    const double sxi[2] = {eta[0], eta[1]};

    double mxi[2] = {0.0, 0.0};
    double projalpha = 0.0;
    bool isprojected = false;

    // reset projection info
    int uniqueMaEle = -1;
    double uniqueProjAlpha = 0.0;
    double uniqueMxi[2] = { 0.0, 0.0 };

    // --> find master elements with an unique projection
    for ( unsigned nummaster=0; nummaster<msize; ++nummaster )
    {
      DRT::Element::DiscretizationType mastertype = meles[nummaster]->Shape();
      // project Gauss point onto master element
      MORTAR::MortarProjector::Impl(sele,*meles[nummaster])->ProjectGaussPoint3D(
          sele,sxi,*meles[nummaster],mxi,projalpha);

      bool is_on_mele=true;

      // check GP projection
      const double tol = 0.0;
      switch ( mastertype )
      {
        case DRT::Element::quad4:
        case DRT::Element::quad8:
        case DRT::Element::quad9:
        {
          if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
          {
            is_on_mele=false;
          }
          break;
        }
        default:
        {
          if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
          {
            is_on_mele=false;
          }
          break;
        }
      }

      // gp is valid and the current master element is the first feasible one
      if (is_on_mele==true and !isprojected)
      {
        isprojected = true;
        uniqueMaEle = nummaster;
        uniqueProjAlpha = projalpha;
        uniqueMxi[0] = mxi[0];
        uniqueMxi[1] = mxi[1];
      }
      // found a second master element with a feasible projection
      else if (is_on_mele and isprojected)
      {
        if (projalpha < uniqueProjAlpha)
        {
          uniqueMaEle = nummaster;
          uniqueProjAlpha = projalpha;
          uniqueMxi[0] = mxi[0];
          uniqueMxi[1] = mxi[1];
        }
      }
    }//mele loop

    if ( uniqueMaEle != -1 )
    {
      MORTAR::MortarElement* mele = meles[ uniqueMaEle ];
      if ( projInfo.find( mele ) == projInfo.end() )
      {
        projInfo[ mele ] = UniqueProjInfo( numGP, msize );
      }

      projInfo[ mele ].Insert( gp, uniqueProjAlpha, uniqueMxi );
    }
    else if (uniqueMaEle == -1 && (not boundary_ele))
        std::cout << "*** warning *** Non-boundary element has non-projectable Gauss point \n" ;
  }

  return ( projInfo.size() > 0 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::INTEGRATOR::FindFeasibleMasterElement2D(
    MORTAR::MortarElement& sele,
    const std::vector<MORTAR::MortarElement*> & meles,
    CoIntegrator& wrapper,
    UniqueProjInfoPair& projInfo )
{
  const unsigned msize = meles.size();
  const unsigned numGP = wrapper.nGP();


  projInfo.clear();
  if ( msize > projInfo.capacity() )
    projInfo.resize( msize );

  for ( unsigned gp=0; gp<numGP; ++gp )
  {
    bool kink_projection=false;

    // coordinates and weight
    const double eta = wrapper.Coordinate(gp,0);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    const double sxi[2] = {eta, 0.0};

    // reset projection info
    int uniqueMaEle = -1;
    double uniqueMxi[2] = { 0.0, 0.0 };

    // loop over all Master Elements
    for ( unsigned nummaster=0; nummaster<msize; ++nummaster )
    {
      // project Gauss point onto master element
      double mxi[2] = {0.0, 0.0};
      MORTAR::MortarProjector::Impl(sele,*meles[nummaster])->ProjectGaussPoint2D(sele,sxi,
          *meles[nummaster],mxi);

      if ((mxi[0]>=-1.0) and (mxi[0]<=1.0) and (not kink_projection))
      {
        kink_projection = true;
        uniqueMaEle = nummaster;
        uniqueMxi[0] = mxi[0];
      }

      if ( kink_projection )
        break;
    }

    if ( uniqueMaEle != -1 )
    {
      MORTAR::MortarElement* mele = meles[ uniqueMaEle ];
      if ( projInfo.find( mele ) == projInfo.end() )
      {
        projInfo[ mele ] = UniqueProjInfo( numGP, msize );
      }

      projInfo[ mele ].Insert( gp, uniqueMxi );
    }
  }

  return ( projInfo.size() > 0 );
}
