/*!
\file intersection_service.H

\brief collection of service methods for intersection computations

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

#ifdef CCADISCRET
#ifndef INTERSECTION_SERVICE_H_
#define INTERSECTION_SERVICE_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_geometry/vector_definitions.H"
#include "blitz_tiny_operation.H"


namespace GEO
{

    //! based on this element property, one can speed up geometry algorithms 
    enum EleGeoType {CARTESIAN, LINEAR, HIGHERORDER};

    //! tolerance used for residuals in Newton-methods
    const double                 TOL14 = 1e-14;
    
    //! tolerance used for residuals in Newton-methods
    const double                 TOL13 = 1e-13;
    
    //! tolerance used for residuals in Newton-methods
    const double                 TOL10 = 1e-10;

    //! named tolerance for easy search/grep
    const double                 TOL7 = 1e-7;
    //! named tolerance for easy search/grep
    const double                 TOL4 = 1e-4;
    //! named tolerance for easy search/grep
    const double                 TOL3 = 1e-3;
    //! named tolerance for easy search/grep
    const double                 TOL2 = 1e-2;
    //! named tolerance for easy search/grep
    const double                 TOLPLUS8 = 1e8;
    //! named tolerance for easy search/grep
    const double                 TOLPLUS15 = 1e15;   
    //! large number to start computations of nearest distance in tree
    const double                 LARGENUMBER = 1e18;
    
    /*!
    \brief calculate a*a
    \param a
    \return a^2
     */
    template <typename T>
    inline T sqr(const T a)
    {
        return a * a;
    }
    
    
    /*!
    \brief computes the cross product of two BlitzVecs a x b
    \param a (in) : arbitrary BlitzVec     
    \param b (in) : arbitrary BlitzVec
    \return cross product of two BlitzVecs a x b
    */
    BlitzVec computeCrossProduct(
        const BlitzVec& a,
        const BlitzVec& b
        );
    
    
    /*!
    \brief calculate quadratic norm
    \param v (in) : BlitzVec
    \return normalized BlitzVec 
    */   
    inline double Norm2(   
        const BlitzVec&     v)
    {
        blitz::firstIndex i;
        return std::sqrt(blitz::sum(v(i)*v(i)));
    }
     
    
    /*!
    \brief calculate quadratic norm
    \param v (in) : BlitzVec
    \return normalized BlitzVec 
    */   
    inline double Norm2(   
        const blitz::TinyVector<double,1>&     v)
    {
        return abs(v(0));
    }

    
    /*!
    \brief calculate quadratic norm
    \param v (in) : BlitzVec
    \return normalized BlitzVec 
    */   
    inline double Norm2(   
        const blitz::TinyVector<double,3>&     v)
    {
        return sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));
    }
    
    
    /*!
    \brief normalizes a vector, such that the resulting length is 1
    \param v (in/out) : BlitzVec
    */   
    inline void normalizeVectorInPLace(
        BlitzVec&     v 
        )
    {
        v /= Norm2(v);
        return;
    }
    
    
    /*!
    \brief sum up the absolute value of all entries of a vector
    \param v (in) : BlitzVec
    \return sum
    */   
    inline double SumOfFabsEntries(
        const BlitzVec3&     v  
        )
    {
        return abs(v(0)) + abs(v(1)) + abs(v(2));
    }
    
    
    /*!
    \brief  computes the Theorem of Pythagoras
            (a^2 + b^2)^(1/2)
            (modified from NUMERICAL RECIPES)
                   
    \param a    (in)        : a
    \param b    (in)        : b
    \return result of the Theorem of Pythagoras
    */    
    inline double pythagoras(
        const double  a, 
        const double  b
        )
    {
        //note: sqrt(0.0) is properly defined and returns 0.0 as expected
        //      no check for zero needed
        return sqrt(sqr(a)+sqr(b));
    }
    
    
    /*!
    \brief compares two points (overloaded method)
    \param point1       (in)    : first point   (double*)
    \param point2       (in)    : second point  (double*)
    \return true if both points equal each other, false otherwise
    */
    template<int length>
    bool comparePoints(
        const double*     point1,   
        const double*     point2   
        ) 
    {   
        for(int i = 0; i < length; i++)
            if(fabs(point1[i] - point2[i]) > TOL7)
            {
                return false;
            }
        return true;
    }
    
    
    /*!
     \brief compares two points (overloaded method)
 
     \param point1       (in)    : first point   (vector<double>)
     \param point2       (in)    : second point  (double*)
     \return true if both points equal each other, false otherwise
     */
    template <int length>
    bool comparePoints(   
        const vector<double>&     point1,
        const double*             point2) 
    {   
        dsassert(point1.size() == length, "not of proper length!");
        for(unsigned int i = 0; i < point1.size() ; ++i)
            if(fabs(point1[i] - point2[i]) > TOL7)
            {
                return false;
            }
        return true;
    }
     
     
     /*!
     \brief compares two points (overloaded method)
     
     \param point1       (in)    : first point   (vector<double>)
     \param point2       (in)    : second point  (vector<double>)
     \return true if both points equal each other, false otherwise
     */
     template <int length>
     bool comparePoints(   
         const vector<double>& point1,
         const vector<double>& point2) 
     {   
         dsassert(point1.size() == length, "not of proper length!");
         dsassert(point2.size() == length, "not of proper length!");
         for(unsigned int i = 0; i < length ; ++i)
             if(fabs(point1[i] - point2[i]) > TOL7)
             {
                 return false;
             }
         return true;
     }
     
     
     /*!
     \brief compares two points (overloaded method)
     
     \param point1       (in)    : first point   (vector<double>)
     \param point2       (in)    : second point  (vector<double>)
     \param tol          (in)    : tolerance defining whether positions are equal
     \return true if both points equal each other, false otherwise
     */
     template <int length>
     bool comparePoints(   
         const vector<double>& point1,
         const vector<double>& point2,
         const double          tol) 
     {   
         dsassert(point1.size() == length, "not of proper length!");
         dsassert(point2.size() == length, "not of proper length!");
         for(unsigned int i = 0; i < length ; ++i)
         {
             if(fabs(point1[i] - point2[i]) > tol)
             {
                 return false;
             }
         }
         return true;
     }
     
     
     /*!
     \brief compares two points (overloaded method)
 
     \param point1       (in)    : first point   (BlitzVec)
     \param point2       (in)    : second point  (BlitzVec)
     \return true if both points equal each other, false otherwise
     */
     template <int length>
     bool comparePoints(    
         const BlitzVec&     point1,
         const BlitzVec&     point2) 
     {   
         dsassert(point1.size() == length, "not of proper length!");
         dsassert(point2.size() == length, "not of proper length!");
         for(int i = 0; i < length ; i++)
             if(fabs(point1(i) - point2(i)) > TOL7)
             {
                 return false;
             }
         return true;
     }
     
     
     /*!
     \brief compares two points (overloaded method)
 
     \param point1       (in)    : first point   (BlitzVec)
     \param point2       (in)    : second point  (BlitzVec)
     \return true if both points equal each other, false otherwise
     */
     template <int length>
     bool comparePoints(    
         const blitz::TinyVector<double,length>&     point1,
         const blitz::TinyVector<double,length>&     point2) 
     {   
         for(int i = 0; i < length ; i++)
             if(fabs(point1(i) - point2(i)) > TOL7)
             {
                 return false;
             }
         return true;
     }
     
     
    /*!
    \brief Checks if a node is within an 
           extended axis-aligned bounding box (XAABB)

    \param pos        (in)    : node to be checked
    \param XAABB      (in)    : extended axis-aligned bouning box
    \return true if node is within the XAABB or false otherwise
    */
    bool isPositionWithinXAABB(
            const BlitzVec3&                   pos,    
            const BlitzMat3x2&                 XAABB   
            );
    
    
    /*!
    \brief Checks if a straightline with 2 nodes is within an 
           extended axis-aligned bounding box (XAABB)
    \param pos1       (in)    : node to be checked
    \param pos2       (in)    : node to be checked
    \param XAABB      (in)    : extended axis-aligned bouning box
    \return true if node is within the XAABB or false otherwise
    */
    bool isLineWithinXAABB(
            const BlitzVec3&                   pos1, 
            const BlitzVec3&                   pos2, 
            const BlitzMat3x2&                 XAABB 
            );
    
    
    /*!
    \brief Computes a rough overestimating extended  
           axis-aligned bounding box for an element (XAABB)
    \param element        (in)  element
    \param xyze           (in)  nodal position array (3,numnode)
    \param eleGeoType     (in)  element geometric type CARTESIAN LINEAR or HIGHERORDER
    \return extended axis-aligned bounding box  (XAABB) for an element
     */
    BlitzMat3x2 computeFastXAABB( 
        const DRT::Element*   element,    ///< element for which one wants to have a boundaing box
        const BlitzMat&       xyze,       ///< nodal position array (3,numnode)
        const EleGeoType      eleGeoType  ///< element geometric type CARTESIAN LINEAR or HIGHERORDER
        );


    /*!
    \brief checks if two extended axis aligned bounding boxes intersect

    \param cutterXAABB (in)         : XAABB of the cutting element
    \param xfemXAABB   (in)         : XAABB of the xfem element
    \return true if the XAABB's intersect or false otherwise
     */
    bool intersectionOfXAABB(
        const BlitzMat3x2&    cutterXAABB,
        const BlitzMat3x2&    xfemXAABB
        );


    /*!
    \brief checks if an element is Cartesian, linear or higherorder

    \param element        (in)         : element
    \param xyze_element   (in)         : coordinates of the element
    \param eleGeoType     (out)        : element geometric type CARTESIAN LINEAR or HIGHERORDER
     */
    void checkGeoType(
        DRT::Element*                element,
        const BlitzMat               xyze_element,
        EleGeoType&                  eleGeoType
        );

    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
    
    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates
    */  
    bool checkPositionWithinElement(  
        const DRT::Element*                 element,
        const BlitzVec3&                    x
        );   
    
    
    /*!
    \brief checks if a position in element coordinates lies within a certain Element parameter space
    \param eleCoord              (out)  node in element coordinates (r,s,t)
    \param distype               (in)   shape of the element
    \return  true if point lies in the element, false otherwise
    */
    template <class V>
    bool checkPositionWithinElementParameterSpace(
            const V&                               eleCoord, 
            const DRT::Element::DiscretizationType distype)
    {
      bool nodeWithinElement = true;
      switch(distype)
      {
        case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27: 
        case DRT::Element::quad4: case DRT::Element::quad8:  case DRT::Element::quad9:
        case DRT::Element::line2: case DRT::Element::line3:    
        {
          for(int i = 0; i < DRT::UTILS::getDimension(distype); i++)
            if(eleCoord(i) > (1.0+GEO::TOL7) || eleCoord(i) < (-1)*(1.0+GEO::TOL7) ) 
            {
              nodeWithinElement = false;
              break;
            }
          break;
        }
        default:
          dserror("distype not yet implemented");
        }   
        
      return nodeWithinElement;
    }
    
    
    /*!
    \brief transforms a point in element coordinates to a point
           in current coordinates with respect to a given element
    \param element              (in)        : element 
    \param xyze                 (in)        : nodal positions of element
    \param xsi                  (in)        : element coordinates
    \param x                    (out)       : position in physical coordinates (x, y, z)
    */
    template <class V>
    void elementToCurrentCoordinates(   
        const DRT::Element*  element,  
        const BlitzMat&      xyze,    
        const V&             xsi,      
        BlitzVec3&           x      )
    {
        const int numNodes = element->NumNode();
        BlitzVec funct(numNodes);
        
        switch(DRT::UTILS::getDimension(element->Shape()))
        {
            case 1:
            {
                DRT::UTILS::shape_function_1D(funct, xsi(0), element->Shape());
                break;
            }
            case 2:
            {
                DRT::UTILS::shape_function_2D(funct, xsi(0), xsi(1), element->Shape());
                break;
            }
            case 3:
            {
                DRT::UTILS::shape_function_3D(funct, xsi(0), xsi(1), xsi(2), element->Shape());
                break;
            }
            default:
                dserror("dimension of the element is not correct");
        }          
            
        x = 0.0;
        for(int i=0; i<numNodes; i++)
        {
            for(int j=0; j<3; j++)
              x(j) += xyze(j,i) * funct(i);
        }
        
        return;
    }
    
    
    /*!
    \brief transforms a point in element coordinates to a point
           in current coordinates with respect to a given element
    \param element              (in)        : element 
    \param xyze                 (in)        : nodal positions of element
    \param eleCoord             (out)       : element coordinates
    */
    void elementToCurrentCoordinatesInPlace(   
        const DRT::Element*  element,   
        const BlitzMat&      xyze,       
        BlitzVec&            eleCoord   
        );

    
    /*!
    \brief transforms a point in current coordinates to a point
           in element coordinates with respect to a given element       
           The nonlinear system of equation is solved with help of the Newton-method.

    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates (x, y, z)
    \param xsi                  (out)       : node in element coordinates (r, s, t)
    \return bool, indicates, whether positions was found within element coordinates
    */  
    bool currentToVolumeElementCoordinates(
        const DRT::Element*                 element,
        const BlitzVec3&                    x,
        BlitzVec3&                          xsi
        );
    
    
    /*!
    \brief transforms a point in current coordinates to a point
           in element coordinates with respect to a given element
           rounds point within a tolerance around 1 and -1 to 1.0 and -1.0 (required by tetgen)
    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates (x, y, z)
    \param tol                  (out)       : tolerance
    \return                                   vector of element coordinates
    */
    inline BlitzVec3 currentToVolumeElementCoordinatesExact(
            const DRT::Element*         element,  
            const BlitzVec3&            x,        
            const double                tol       
            )
    {
      static BlitzVec3 xsi;
      currentToVolumeElementCoordinates(element, x, xsi);
     
      // rounding 1 and -1 to be exact for the CDT
      // TODO extend for triangles
      for(int j = 0; j < 3; j++)
      {
          if( fabs((fabs(xsi(j))-1.0)) < tol &&  xsi(j) < 0)    xsi(j) = -1.0;
          if( fabs((fabs(xsi(j))-1.0)) < tol &&  xsi(j) > 0)    xsi(j) =  1.0;      
      }
      return xsi;
    } 
    
   
    /*!
    \brief calculate surface element normal vector in gausspoint (r,s)         
    
    \param deriv              (in)        : shape functions derivative in (r,s) -> (2,numnode)
    \param xyze               (in)        : element nodal coordinates in (x,y,z) space -> (3,numnode)
    \param unrm               (out)       : unit normal vector at (r,s)
    */  
    template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
    inline void computeNormalToSurfaceElementViaShape(
            const M1&       deriv,  
            const M2&       xyze,     
            BlitzVec3&      unrm)
    {
        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
          
        // compute dXYZ / drs
        static BlitzMat3x2 dxyzdrs;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                dxyzdrs(i,j) = 0.0;
                for (int k = 0; k < numnode; ++k)
                {
                    dxyzdrs(i,j) += xyze(i,k)*deriv(j,k);
                }
            }
        }

        /* compute covariant metric tensor G for surface element
        **                        | g11   g12 |
        **                    G = |           |
        **                        | g12   g22 |
        ** where (o denotes the inner product, xyz a vector)
        **
        **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
        ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
        **        dr     dr            dr     ds            ds     ds
        */
        static BlitzMat2x2 metrictensor;
        //metrictensor = blitz::sum(dxyzdrs(k,i)*dxyzdrs(k,j),k);
        BLITZTINY::MtM_product<2,2,3>(dxyzdrs,dxyzdrs,metrictensor);
        
        unrm(0) = dxyzdrs(1,0) * dxyzdrs(2,1) - dxyzdrs(2,0) * dxyzdrs(1,1);
        unrm(1) = dxyzdrs(2,0) * dxyzdrs(0,1) - dxyzdrs(0,0) * dxyzdrs(2,1);
        unrm(2) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(1,0) * dxyzdrs(0,1);
        
        // scale
        const double inv_length = 1.0 / sqrt(unrm(0)*unrm(0) + unrm(1)*unrm(1) + unrm(2)*unrm(2));

        unrm *= inv_length;
        
        return;
    }
    
    
    /*!
    \brief calculate element normal in gausspoint (r,s)
    \param xyze_surfaceElement   (in)   nodal position array (3,numnode)
    \param eleCoord              (in)   node in element coordinates (r, s)
    \param unitnormalvec         (out)  unit normal vector
    */
    template<DRT::Element::DiscretizationType DISTYPE>
    inline void computeNormalToSurfaceElement(
        const BlitzMat&          xyze_surfaceElement,               
        const BlitzVec2&         eleCoord,               
        BlitzVec3&               unitnormalvec)
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      static blitz::TinyMatrix<double,2,numnode> deriv_boundary;
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary,eleCoord(0),eleCoord(1),DISTYPE);
      
      computeNormalToSurfaceElementViaShape<DISTYPE>(deriv_boundary, xyze_surfaceElement, unitnormalvec);
      return;
    }
    
    
    /*!
    \brief calculate element normal in gausspoint (r,s)
    \param surfaceElement        (in) :  surface element
    \param xyze_surfaceElement   (in) :  nodal position array (3,numnode)
    \param eleCoord              (in) :  node in element coordinates (r, s)
    \param unitnormalvec         (out):  unit normal vector
    */
    inline void computeNormalToSurfaceElement(
        const DRT::Element*      surfaceElement,    
        const BlitzMat&          xyze_surfaceElement,                
        const BlitzVec2&         eleCoord,               
        BlitzVec3&               unitnormalvec)
    {
      switch (surfaceElement->Shape())
      {
          case DRT::Element::quad4:
              computeNormalToSurfaceElement<DRT::Element::quad4>(xyze_surfaceElement, eleCoord, unitnormalvec);
              break;
          case DRT::Element::quad8:
              computeNormalToSurfaceElement<DRT::Element::quad8>(xyze_surfaceElement, eleCoord, unitnormalvec);
              break;
          case DRT::Element::quad9:
              computeNormalToSurfaceElement<DRT::Element::quad9>(xyze_surfaceElement, eleCoord, unitnormalvec);
              break;
          case DRT::Element::tri3:
              computeNormalToSurfaceElement<DRT::Element::tri3>(xyze_surfaceElement, eleCoord, unitnormalvec);
              break;
          default:
              cout << DRT::DistypeToString(surfaceElement->Shape()) << " not found!" << endl;
              dserror("please add your distype!");
      }
      return;
    }
    
    
    /*!
    \brief  finds the nearest point for a given point on a surface element
            the element coordinates of this point are stored in eleCoord, the normal
            point from the given point to the nearest point on the surface,
            distance returns the absolut length of the normal. 
            If the computed point lies outside the element , the distance which is retuned equals -1
            and the normal is set to zero.
    \param surfaceElement        (in) :   surface element
    \param xyze_surfaceElement   (in) :   nodal position array (3,numnode)
    \param physCoord             (in) :   node in physical coordinates (x, y, z)
    \param eleCoord              (out):   node in element coordinates (r, s)
    \param normal                (out):   normal point from a given point to the surface point
    \param distance              (out):   length of the normal
    */         
    bool searchForNearestPointOnSurface(
            const DRT::Element*                     surfaceElement,  
            const BlitzMat&                         xyze_surfaceElement,     
            const BlitzVec3&                        physCoord,       
            BlitzVec2&                              eleCoord,        
            BlitzVec3&                              normal,          
            double&                                 distance        
            );
    
    
    /*!
     \brief checks if a position in current coordinates lies within a certain surfaceElement           
            The nonlinear system of equation is solved with help of the Newton-method.
     \param surfaceElement        (in) :  surface element
     \param xyze_surfaceElement   (in) :  nodal position array (3,numnode)
     \param physCoord             (in) :  node in physical coordinates (x, y, z)
     \param eleCoord              (out):  node in element coordinates (r, s)
     */  
    void CurrentToSurfaceElementCoordinates(
        const DRT::Element*      surfaceElement,        
        const BlitzMat&          xyze_surfaceElement,   
        const BlitzVec3&         physCoord,             
        BlitzVec2&               eleCoord                  
        );


    /*!
    \brief checks if a position in current coordinates lies within a certain lineElement           
           The nonlinear system of equation is solved with help of the Newton-method.
    \param lineElement            (in)  : line element
    \param xyze_surfaceElement    (in)  : nodal position array (3,numnode)
    \param physCoord              (in)  : node in physical coordinates (x, y, z)
    \param eleCoord               (out) : node in element coordinates (r)
    */
    void CurrentToLineElementCoordinates(
            const DRT::Element*         lineElement,
            const BlitzMat&             xyze_lineElement,
            const BlitzVec3&            physCoord,
            BlitzVec1&                  eleCoord
            );
    

}


#endif  //INTERSECTION_MATH_H_
#endif  // #ifdef CCADISCRET


