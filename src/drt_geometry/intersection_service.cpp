/*!----------------------------------------------------------------------
\file intersection_service.cpp

\brief collection of service methods for intersection computations

      
<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "intersection_service.H"
#include "../drt_lib/drt_element.H"



/*----------------------------------------------------------------------*
 |  ML:     computes the cross product                       u.may 08/07|
 |          of 2 vectors c = a x b                                      |
 *----------------------------------------------------------------------*/  
LINALG::Matrix<3,1> GEO::computeCrossProduct(
    const LINALG::Matrix<3,1>& a,
    const LINALG::Matrix<3,1>& b)
{
    LINALG::Matrix<3,1> c;
   
    c(0) = a(1)*b(2) - a(2)*b(1);
    c(1) = a(2)*b(0) - a(0)*b(2);
    c(2) = a(0)*b(1) - a(1)*b(0);
    
    return c;
}



/*----------------------------------------------------------------------*
 |  ICS:    computes an extended axis-aligned bounding box   u.may 06/07|
 |          XAABB for a given element                                   |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,2> GEO::computeFastXAABB( 
    const DRT::Element*                 element,
    const LINALG::SerialDenseMatrix&    xyze,
    const EleGeoType                    eleGeoType)
{
    const int nsd = 3;
    LINALG::Matrix<3,2> XAABB;
    
    // first node
    for(int dim=0; dim<nsd; ++dim)
    {
        XAABB(dim, 0) = xyze(dim, 0) - TOL7;
        XAABB(dim, 1) = xyze(dim, 0) + TOL7;
    }
    // remaining nodes
    for(int i=1; i<element->NumNode(); ++i)
    {
        for(int dim=0; dim<nsd; dim++)
        {
            XAABB(dim, 0) = std::min( XAABB(dim, 0), xyze(dim,i) - TOL7);
            XAABB(dim, 1) = std::max( XAABB(dim, 1), xyze(dim,i) + TOL7);
        }
    }
    
    if(eleGeoType == HIGHERORDER)
    {
      double maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
      for(int dim=1; dim<nsd; ++dim)
        maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );
    
      // subtracts half of the maximal distance to minX, minY, minZ
      // adds half of the maximal distance to maxX, maxY, maxZ 
      const double halfMaxDistance = 0.5*maxDistance;
      for(int dim=0; dim<nsd; ++dim)
      {
          XAABB(dim, 0) -= halfMaxDistance;
          XAABB(dim, 1) += halfMaxDistance;
      }   
    }
    /*
    printf("\n");
    printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
              XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
    printf("\n");
    */
    
    return XAABB;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if two XAABB's intersect                  u.may 06/07|
 *----------------------------------------------------------------------*/
bool GEO::intersectionOfXAABB(  
    const LINALG::Matrix<3,2>&     cutterXAABB, 
    const LINALG::Matrix<3,2>&     xfemXAABB)
{
    
  /*====================================================================*/
  /* bounding box topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (x,y,z) of nodes
   * node 0: (minX, minY, minZ)
   * node 1: (maxX, minY, minZ)
   * node 2: (maxX, maxY, minZ)
   * node 3: (minX, maxY, minZ)
   * node 4: (minX, minY, maxZ)
   * node 5: (maxX, minY, maxZ)
   * node 6: (maxX, maxY, maxZ)
   * node 7: (minX, maxY, maxZ)
   * 
   *                      z
   *                      |           
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/
    
    bool intersection =  false;
    static std::vector< LINALG::Matrix<3,1> > nodes(8, LINALG::Matrix<3,1>());
    
    nodes[0](0) = cutterXAABB(0,0); nodes[0](1) = cutterXAABB(1,0); nodes[0](2) = cutterXAABB(2,0); // node 0   
    nodes[1](0) = cutterXAABB(0,1); nodes[1](1) = cutterXAABB(1,0); nodes[1](2) = cutterXAABB(2,0); // node 1
    nodes[2](0) = cutterXAABB(0,1); nodes[2](1) = cutterXAABB(1,1); nodes[2](2) = cutterXAABB(2,0); // node 2
    nodes[3](0) = cutterXAABB(0,0); nodes[3](1) = cutterXAABB(1,1); nodes[3](2) = cutterXAABB(2,0); // node 3
    nodes[4](0) = cutterXAABB(0,0); nodes[4](1) = cutterXAABB(1,0); nodes[4](2) = cutterXAABB(2,1); // node 4
    nodes[5](0) = cutterXAABB(0,1); nodes[5](1) = cutterXAABB(1,0); nodes[5](2) = cutterXAABB(2,1); // node 5
    nodes[6](0) = cutterXAABB(0,1); nodes[6](1) = cutterXAABB(1,1); nodes[6](2) = cutterXAABB(2,1); // node 6
    nodes[7](0) = cutterXAABB(0,0); nodes[7](1) = cutterXAABB(1,1); nodes[7](2) = cutterXAABB(2,1); // node 7
    
    for (int i = 0; i < 8; i++)
        if(isPositionWithinXAABB(nodes[i], xfemXAABB))
        {
            intersection = true;
            break;
        }
    
    if(!intersection)
    {
        for (int i = 0; i < 12; i++)
        {
            const int index1 = DRT::UTILS::eleNodeNumbering_hex27_lines[i][0];
            const int index2 = DRT::UTILS::eleNodeNumbering_hex27_lines[i][1];
            if(isLineWithinXAABB(nodes[index1], nodes[index2], xfemXAABB))
            {
                intersection = true;
                break;
            }
        }
    }
    
    if(!intersection)
    {
        nodes[0](0) = xfemXAABB(0,0);   nodes[0](1) = xfemXAABB(1,0);   nodes[0](2) = xfemXAABB(2,0);   // node 0   
        nodes[1](0) = xfemXAABB(0,1);   nodes[1](1) = xfemXAABB(1,0);   nodes[1](2) = xfemXAABB(2,0);   // node 1
        nodes[2](0) = xfemXAABB(0,1);   nodes[2](1) = xfemXAABB(1,1);   nodes[2](2) = xfemXAABB(2,0);   // node 2
        nodes[3](0) = xfemXAABB(0,0);   nodes[3](1) = xfemXAABB(1,1);   nodes[3](2) = xfemXAABB(2,0);   // node 3
        nodes[4](0) = xfemXAABB(0,0);   nodes[4](1) = xfemXAABB(1,0);   nodes[4](2) = xfemXAABB(2,1);   // node 4
        nodes[5](0) = xfemXAABB(0,1);   nodes[5](1) = xfemXAABB(1,0);   nodes[5](2) = xfemXAABB(2,1);   // node 5
        nodes[6](0) = xfemXAABB(0,1);   nodes[6](1) = xfemXAABB(1,1);   nodes[6](2) = xfemXAABB(2,1);   // node 6
        nodes[7](0) = xfemXAABB(0,0);   nodes[7](1) = xfemXAABB(1,1);   nodes[7](2) = xfemXAABB(2,1);   // node 7
    
        for (int i = 0; i < 8; i++)
            if(isPositionWithinXAABB(nodes[i], cutterXAABB))
            {
                intersection = true;
                break;
            }
    }   
    
    if(!intersection)
    {
        for (int i = 0; i < 12; i++)
        {
            const int index1 = DRT::UTILS::eleNodeNumbering_hex27_lines[i][0];
            const int index2 = DRT::UTILS::eleNodeNumbering_hex27_lines[i][1];
            if(isLineWithinXAABB(nodes[index1], nodes[index2], cutterXAABB))
            {
                intersection = true;
                break;
            }
        }
    }
    return intersection;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if an element is CARTESIAN, LINEAR and    u.may 07/08|
 |          HIGHERORDER                                                 |
 *----------------------------------------------------------------------*/
void GEO::checkGeoType(
           DRT::Element*                      element,
           const LINALG::SerialDenseMatrix&   xyze_element,
           EleGeoType&                        eleGeoType)
{
  bool cartesian = true;
  int CartesianCount = 0;
  const int dimCoord = 3;
  const DRT::Element::DiscretizationType distype = element->Shape();
  const int eleDim = DRT::UTILS::getDimension(distype);
  
  if(DRT::UTILS::getOrder(distype) ==1)
    eleGeoType = LINEAR;
  else if(DRT::UTILS::getOrder(distype)==2)
    eleGeoType = HIGHERORDER;
  else
    dserror("order of element shapefuntion is not correct");
  
  // check if cartesian
  if(eleDim == 3)
  {
    const vector< vector<int> > eleNodeNumbering = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
    vector< RCP<DRT::Element> >surfaces = element->Surfaces();
    for(int i = 0; i < element->NumSurface(); i++)
    {      
      CartesianCount = 0;
      const DRT::Element* surfaceP = surfaces[i].get();
  
      for(int k = 0; k < dimCoord; k++)
      { 
        int nodeId = eleNodeNumbering[i][0];
        const double nodalcoord =  xyze_element(k,nodeId);
        for(int j = 1; j < surfaceP->NumNode(); j++)
        {
          nodeId = eleNodeNumbering[i][j];
          if(fabs(nodalcoord - xyze_element(k,nodeId)) > TOL7)
          {
            CartesianCount++;
            break;
          } 
        }
      }
      if(CartesianCount > 2)  
      {
        cartesian = false;
        break;
      }
    } // for xfem surfaces
  } // if eleDim == 3
  else if(eleDim == 2)
  {
    CartesianCount = 0;
    for(int k = 0; k < dimCoord; k++)
    { 
      const double nodalcoord =  xyze_element(k,0);
      for(int j = 1; j < element->NumNode(); j++)
      {
        if(fabs(nodalcoord - xyze_element(k,j)) > TOL7)
        {
          CartesianCount++;
          break;
        } 
      }
    }
    if(CartesianCount > 2)  
      cartesian = false;
  }
  else
    dserror("dimension of element is not correct");

  
  
  if(cartesian)
    eleGeoType = CARTESIAN;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if a position is within an XAABB          u.may 06/07|
 *----------------------------------------------------------------------*/
bool GEO::isPositionWithinXAABB(    
    const LINALG::Matrix<3,1>&                    pos,
    const LINALG::Matrix<3,2>&                    XAABB)
{
    bool isWithin = true;
    for (int i=0; i<3; i++)
    {
        const double diffMin = XAABB(i,0) - TOL7;
        const double diffMax = XAABB(i,1) + TOL7;
        
       // printf("nodal value =  %f, min =  %f, max =  %f\n", node[dim], diffMin, diffMax);   
        if((pos(i) < diffMin)||(pos(i) > diffMax)) //check again !!!!!   
        {
            isWithin = false;
            break;
        }
    }
   
    return isWithin;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if a pos is within an XAABB               u.may 06/07|
 *----------------------------------------------------------------------*/
bool GEO::isLineWithinXAABB(    
    const LINALG::Matrix<3,1>&                  pos1,
    const LINALG::Matrix<3,1>&                  pos2,
    const LINALG::Matrix<3,2>&                  XAABB)
{
    const int nsd = 3;
    bool isWithin = true;
    int isd = -1;
    
    for(isd=0; isd<nsd; isd++)
        if(fabs(pos1(isd)-pos2(isd)) > TOL7)
            break;
    
    for(int ksd = 0; ksd < nsd; ksd++)
    {
        if(ksd != isd)
        {
            const double min = XAABB(ksd,0) - TOL7;
            const double max = XAABB(ksd,1) + TOL7;
   
            if((pos1(ksd) < min)||(pos1(ksd) > max))
                isWithin = false;
            
        }
        if(!isWithin)
            break;
    }
        
    if(isWithin && isd > -1)
    {
        isWithin = false;
        const double min = XAABB(isd,0) - TOL7;
        const double max = XAABB(isd,1) + TOL7;
                            
        if( ((pos1(isd) < min) && (pos2(isd) > max)) ||  
            ((pos2(isd) < min) && (pos1(isd) > max)) )
            isWithin = true;
    }
    return isWithin;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within a given element   u.may 06/07|   
 *----------------------------------------------------------------------*/
bool GEO::checkPositionWithinElement(  
    const DRT::Element*                 element,
    const LINALG::SerialDenseMatrix&    xyze,
    const LINALG::Matrix<3,1>&          x)
{
    dsassert(DRT::UTILS::getDimension(element->Shape()) == 3, "only valid for 3 dimensional elements");
    LINALG::Matrix<3,1> xsi(true);
    bool nodeWithinElement = currentToVolumeElementCoordinates(element->Shape(), xyze, x, xsi);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0),xsi(1),xsi(2), residual, TOL14);
    
    nodeWithinElement = checkPositionWithinElementParameterSpace(xsi,element->Shape());

    return nodeWithinElement;
}


/*----------------------------------------------------------------------*
 |  RQI:    searches the nearest point on a surface          u.may 02/08|
 |          element for a given point in physical coordinates           |
 *----------------------------------------------------------------------*/
bool GEO::searchForNearestPointOnSurface(
    const DRT::Element*                     surfaceElement,
    const LINALG::SerialDenseMatrix&        xyze_surfaceElement,
    const LINALG::Matrix<3,1>&              physCoord,
    LINALG::Matrix<2,1>&                    eleCoord,
    LINALG::Matrix<3,1>&                    normal,
    double&                                 distance)
{
  
  CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, physCoord, eleCoord);
  
  const bool pointWithinElement = checkPositionWithinElementParameterSpace(eleCoord, surfaceElement->Shape());
  
  // normal vector at position xsi
  static LINALG::Matrix<3,1> eleNormalAtXsi;
  computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, eleCoord, eleNormalAtXsi);
  
  LINALG::Matrix<3,1> x_surface_phys;
  elementToCurrentCoordinates(surfaceElement->Shape(), xyze_surfaceElement, eleCoord, x_surface_phys);
  // normal pointing away from the surface towards physCoord
  normal.Update(1.0, physCoord, -1.0, x_surface_phys);
  // absolute distance between point and surface
  distance = normal.Norm2();
 
  
  if(fabs(distance) > GEO::TOL7)
  {
    // compute distance with sign
    const double scalarproduct = eleNormalAtXsi(0)*normal(0) + eleNormalAtXsi(1)*normal(1) + eleNormalAtXsi(2)*normal(2);
    const double teiler = eleNormalAtXsi.Norm2() * normal.Norm2();
    const double cosphi = scalarproduct / teiler;
    const double vorzeichen = cosphi/abs(cosphi);
    distance *= vorzeichen;
  }
  else
    distance = 0.0;
  
  return pointWithinElement;
}



#endif  // #ifdef CCADISCRET
