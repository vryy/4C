/*!----------------------------------------------------------------------
\file intersection_service_templates.H

\brief collection of templated service methods for intersection computations

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef INTERSECTION_SERVICE_TEMPLATES_H_
#define INTERSECTION_SERVICE_TEMPLATES_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"

namespace GEO
{


  //! based on this element property, one can speed up geometry algorithms 
  enum EleGeoType {CARTESIAN, LINEAR, HIGHERORDER};



  //! checks if a given surface element lies in a plane
  template <class M>
    static inline bool SurfaceElementIsPlane(   
          const M&                              xyze) 
  {
    static LINALG::Matrix<3,1> elenorm;
    // xyze 3 x 4
  
    if(xyze.N() == 3)
      return true;
    
    elenorm(0)= -(xyze(2,0)*xyze(1,1)) + xyze(1,0)*xyze(2,1) + xyze(2,0)*xyze(1,2) -
                  xyze(2,1)*xyze(1,2)  - xyze(1,0)*xyze(2,2) + xyze(1,1)*xyze(2,2);
    elenorm(1)=   xyze(2,0)*xyze(0,1)  - xyze(0,0)*xyze(2,1) - xyze(2,0)*xyze(0,2) +
                  xyze(2,1)*xyze(0,2)  + xyze(0,0)*xyze(2,2) - xyze(0,1)*xyze(2,2);
    elenorm(2)= -(xyze(1,0)*xyze(0,1)) + xyze(0,0)*xyze(1,1) + xyze(1,0)*xyze(0,2) -
                  xyze(1,1)*xyze(0,2)  - xyze(0,0)*xyze(1,2) + xyze(0,1)*xyze(1,2);
    
    for(int i = 3; i < xyze.N(); i++)
    {
      double distance = (-(elenorm(0)*( xyze(0,0) - xyze(0,i))) + 
                           elenorm(1)*(-xyze(1,0) + xyze(1,i)) - 
                           elenorm(2)*( xyze(2,0) - xyze(2,i)))/(-elenorm.Norm2());
  
      if(fabs(distance) >= GEO::TOL7)
        return false;
    }
    return true;
  }



  //! checks if a given point lies in a plane surface element 
  template <class M1, class M2>
    static inline bool pointInPlaneSurfaceElement(   
          const M1&                              xyze,
          const M2&                              xyze_points) 
  {
    static LINALG::Matrix<3,1> elenorm;
    // xyze 3 x ...
    
    elenorm(0)= -(xyze(2,0)*xyze(1,1)) + xyze(1,0)*xyze(2,1) + xyze(2,0)*xyze(1,2) -
                  xyze(2,1)*xyze(1,2)  - xyze(1,0)*xyze(2,2) + xyze(1,1)*xyze(2,2);
    elenorm(1)=   xyze(2,0)*xyze(0,1)  - xyze(0,0)*xyze(2,1) - xyze(2,0)*xyze(0,2) +
                  xyze(2,1)*xyze(0,2)  + xyze(0,0)*xyze(2,2) - xyze(0,1)*xyze(2,2);
    elenorm(2)= -(xyze(1,0)*xyze(0,1)) + xyze(0,0)*xyze(1,1) + xyze(1,0)*xyze(0,2) -
                  xyze(1,1)*xyze(0,2)  - xyze(0,0)*xyze(1,2) + xyze(0,1)*xyze(1,2);
    
    for(int i = 0; i < 2; i++)
    {
      double distance = (-(elenorm(0)*( xyze(0,0) - xyze_points(0,i))) + 
                           elenorm(1)*(-xyze(1,0) + xyze_points(1,i)) - 
                           elenorm(2)*( xyze(2,0) - xyze_points(2,i)))/(-elenorm.Norm2());
  
      if(fabs(distance) >= GEO::TOL7)
        return false;
    }
    return true;
  }



  //! template of extended aligned bounding boxes
  template <  int ndim,
              DRT::Element::DiscretizationType distype, 
              class M>
    static inline LINALG::Matrix<3,2> computeFastXAABBT( 
      const M&                            xyze,
      const GEO::EleGeoType               eleGeoType)
  {
      LINALG::Matrix<3,2> XAABB;
      
      // first node
      for(int dim=0; dim<3; ++dim)
      {
          XAABB(dim, 0) = xyze(dim, 0) - TOL7;
          XAABB(dim, 1) = xyze(dim, 0) + TOL7;
      }
      // remaining nodes
      const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      for(int i=1; i < numNodes; ++i)
        for(int dim=0; dim < ndim; dim++)
        {
          XAABB(dim, 0) = std::min( XAABB(dim, 0), xyze(dim,i) - TOL7);
          XAABB(dim, 1) = std::max( XAABB(dim, 1), xyze(dim,i) + TOL7);
        }  
      
      /*
      
      if(eleGeoType == HIGHERORDER)
      {
        double maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
        for(int dim = 1; dim < ndim; ++dim)
          maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );
      
        // subtracts a certain percent of maximal distance to minX, minY, minZ
        // adds a certain percent of the maximal distance to maxX, maxY, maxZ 
        const double extension = 0.2*maxDistance;
        for(int dim=0; dim < ndim; ++dim)
        {
            XAABB(dim, 0) -= extension;
            XAABB(dim, 1) += extension;
        }   
      }
      
      */
      /*
      printf("\n");
      printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
                XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
      printf("\n");
      */
      
      return XAABB;
  }


    
  /*!
  \brief Computes a rough overestimating extended  
         axis-aligned bounding box for an element (XAABB)
  \param distype        (in)  distype of element
  \param xyze           (in)  nodal position array (3,numnode)
  \param eleGeoType     (in)  element geometric type CARTESIAN LINEAR or HIGHERORDER
  \return extended axis-aligned bounding box  (XAABB) for an element
   */
  template <class M>
    LINALG::Matrix<3,2> computeFastXAABB( 
      DRT::Element::DiscretizationType    distype, 
      const M&                            xyze,
      const GEO::EleGeoType               eleGeoType)
  {
    switch (distype)
    {
      case DRT::Element::hex8:  return computeFastXAABBT<3, DRT::Element::hex8> (xyze, eleGeoType);
      case DRT::Element::quad4: return computeFastXAABBT<3, DRT::Element::quad4>(xyze, eleGeoType);
      case DRT::Element::hex20: return computeFastXAABBT<3, DRT::Element::hex20>(xyze, eleGeoType);
      case DRT::Element::hex27: return computeFastXAABBT<3, DRT::Element::hex27>(xyze, eleGeoType);
      case DRT::Element::tet4:  return computeFastXAABBT<3, DRT::Element::tet4> (xyze, eleGeoType);
      case DRT::Element::tet10: return computeFastXAABBT<3, DRT::Element::tet10>(xyze, eleGeoType);
      case DRT::Element::line2: return computeFastXAABBT<3, DRT::Element::line2>(xyze, eleGeoType);
      case DRT::Element::line3: return computeFastXAABBT<3, DRT::Element::line3>(xyze, eleGeoType);
      case DRT::Element::tri3:  return computeFastXAABBT<3, DRT::Element::tri3> (xyze, eleGeoType);
      case DRT::Element::tri6:  return computeFastXAABBT<3, DRT::Element::tri6> (xyze, eleGeoType);
      case DRT::Element::quad8: return computeFastXAABBT<3, DRT::Element::quad8>(xyze, eleGeoType);
      case DRT::Element::quad9: return computeFastXAABBT<3, DRT::Element::quad9>(xyze, eleGeoType);
      default:
          std::cout << DistypeToString(distype) << endl;
          dserror("add your distype to this switch!");
    }   
    return LINALG::Matrix<3,2>(true);
  }
    
    
  //! checks if two XAABB's intersect
  template<int dim>
    static inline bool intersectionOfXAABB(  
      const LINALG::Matrix<3,2>&     cutterXAABB, 
      const LINALG::Matrix<3,2>&     xfemXAABB)
  {
    // check XAABB for intersection
    for(int i = 0; i < dim; i++)
    {
      if(! (  ((cutterXAABB(i,0)  > (xfemXAABB(i,0)   -GEO::TOL7))  && ( cutterXAABB(i,0) < (xfemXAABB(i,1)   +GEO::TOL7)) )  || 
              ((cutterXAABB(i,1)  > (xfemXAABB(i,0)   -GEO::TOL7))  && ( cutterXAABB(i,1) < (xfemXAABB(i,1)   +GEO::TOL7)) )  || 
              ((xfemXAABB(i,0)    > (cutterXAABB(i,0) -GEO::TOL7))  && ( xfemXAABB(i,0)   < (cutterXAABB(i,1) +GEO::TOL7)) )  || 
              ((xfemXAABB(i,1)    > (cutterXAABB(i,0) -GEO::TOL7))  && ( xfemXAABB(i,1)   < (cutterXAABB(i,1) +GEO::TOL7)) )  ))
        return false;
    }  
    return true;
  }


}  // namespace GEO


#endif /*INTERSECTION_SERVICE_TEMPLATES_H_*/
#endif /*CCADISCRET*/



