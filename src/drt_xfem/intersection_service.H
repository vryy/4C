/*!
\file intersection_service.H

\brief collection of math tools for the interface computation of two
       curved meshes

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

#include <blitz/array.h>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Teuchos_RefCountPtr.hpp"
#include "xfem.H"
#include "blitz_tiny_operation.H"


namespace XFEM
{

enum EleGeoType {CARTESIAN, LINEAR, HIGHERORDER};

    //! tolerance 1e-14 used for residuals in Newton-methods
    const double                 TOL14 = 1e-14;
    
    //! tolerance 1e-13 used for residuals in Newton-methods
    const double                 TOL13 = 1e-13;
    
    //! tolerance 1e-13 used for residuals in Newton-methods
    const double                 TOL10 = 1e-10;

    //! tolerance 1e-7 
    const double                 TOL7 = 1e-7;
    
    //! tolerance 1e-3 
    const double                 TOL4 = 1e-4;
        
    //! tolerance 1e-3 
    const double                 TOL3 = 1e-3;
    
    //! tolerance 1e-2 
    const double                 TOL2 = 1e-2;
    
    //! tolerance 1e8
    const double                 TOLPLUS8 = 1e8;
    
    //! tolerance 1e15 
    const double                 TOLPLUS15 = 1e15;
    
    /*!
     * \brief calculate a*a
     * \param a
     * \return a^2
     */
    template <typename T>
    inline T sqr(T a)
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
    \return normalized BlitzVec 
    */   
    inline void normalizeVectorInPLace(
        BlitzVec&     v ///< vector to normalize
        )
    {
        v /= Norm2(v);
        return;
    }
    
    /*!
    \brief sum up the absulte value of all entries of a vector
    \return sum
    */   
    inline double SumOfFabsEntries(
        BlitzVec3&     v  ///< vector 
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
    \brief  fill array with current nodal positions

    \return array with element nodal positions
    */    
    inline BlitzMat getCurrentNodalPositions(
        const DRT::Element*            ele,                      ///< element with nodal pointers
        const map<int,BlitzVec3>&      currentcutterpositions    ///< current positions of all cutter nodes
        )
    {
        const int numnode = ele->NumNode();
        BlitzMat xyze(3,numnode);
        const DRT::Node*const* nodes = ele->Nodes();
        for (int inode = 0; inode < numnode; ++inode)
        {
          const BlitzVec3 x = currentcutterpositions.find(nodes[inode]->Id())->second;
          xyze(0,inode) = x(0);
          xyze(1,inode) = x(1);
          xyze(2,inode) = x(2);
        }
        return xyze;
    }
    
    
    /*!
    \brief transforms a point in element coordinates to a point
           in current coordinates with respect to a given element

    \param element          (in)        : element
    \param eleCoord         (in)    	: vector of element coordinates
    \return                               vector of current coordinates
    */
    template <class V>
    void elementToCurrentCoordinates(   
        const DRT::Element*  element, 
        const BlitzMat&      xyze,
        const V&             xsi,
        BlitzVec3&           x
        )
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

    \param element          (in)        : element
    \param eleCoord         (in/out)    : vector of element/current coordinates
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
    \param xsi                  (out)       : node in element coordinates
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

    \param element          (in)        : element
    \param xsi              (in/out)    : vector of current coordinates/
    \tparam dim                         : dimension of the element parameter space
    \return                               vector of element coordinates
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
        for(int j = 0; j < 3; j++)
        {
            if( fabs((fabs(xsi(j))-1.0)) < tol &&  xsi(j) < 0)    xsi(j) = -1.0;
            if( fabs((fabs(xsi(j))-1.0)) < tol &&  xsi(j) > 0)    xsi(j) =  1.0;      
        }
        return xsi;
    } 
    
    
    /*!
    \brief compares two points (overloaded method)

    \param point1       (in)    : first point   (double*)
    \param point2       (in)    : second point  (double*)
    \param lenght       (in)    : length of arrays
    \return true if both points equal each other, false otherwise
    */
    template<int length>
    bool comparePoints(
        const double*     point1,
        const double*     point2) 
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
     
     
    template <int length>
    bool comparePoints(   
        const vector<double>& point1,
        const vector<double>& point2,
        const double          tol) 
    {   
        dsassert(point1.size() == length, "not of proper length!");
        dsassert(point2.size() == length, "not of proper length!");
        for(unsigned int i = 0; i < length ; ++i)
            if(fabs(point1[i] - point2[i]) > tol)
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

    \return true if node is within the XAABB or false otherwise
    */
    bool isPositionWithinXAABB(
            const BlitzVec3&                   pos,    ///< node to be checked
            const BlitzMat3x2&                 XAABB   ///< extended axis-aligned bouning box
            );
    
    
    /*!
    \brief Checks if a straightline with 2 nodes is within an 
           extended axis-aligned bounding box (XAABB)

    \return true if node is within the XAABB or false otherwise
    */
    bool isLineWithinXAABB(
            const BlitzVec3&                   pos1, ///< node to be checked
            const BlitzVec3&                   pos2, ///< node to be checked
            const BlitzMat3x2&                 XAABB ///< extended axis-aligned bouning box
            );
    
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.
    
    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates
    */  
    bool checkPositionWithinElement(  
        const DRT::Element*                 element,
        const BlitzVec3&                    x
        );   
    
    
    //! calculate element normal vector in gausspoint (r,s)
    template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
    inline void computeNormalToBoundaryElement(
            const M1&       deriv,                  ///< shape functions derivative in (r,s)
            const M2&       xs,                     ///< boundary element nodal coordinates in (x,y,z) space
            BlitzVec3&      unrm                    ///< unit normal vector at (r,s)
            )
    {
        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
          
        // compute dXYZ / drs
        static BlitzMat3x2 dxyzdrs;
        //dxyzdrs = blitz::sum(xs(i,k)*deriv(j,k),k);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                dxyzdrs(i,j) = 0.0;
                for (int k = 0; k < numnode; ++k)
                {
                    dxyzdrs(i,j) += xs(i,k)*deriv(j,k);
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
    
    //! calculate element normal in gausspoint (r,s)
    template<DRT::Element::DiscretizationType DISTYPE>
    inline void computeNormalToBoundaryElement(
        const BlitzMat&          xs,                ///< boundary element nodal coordinates in (x,y,z) space
        const BlitzVec2&         xsi,               ///< coordinates in (r,s) space, where the normal shall be calculated,
        BlitzVec3&               unitnormalvec
            )
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      static blitz::TinyMatrix<double,2,numnode> deriv_boundary;
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary,xsi(0),xsi(1),DISTYPE);
      
      computeNormalToBoundaryElement<DISTYPE>(deriv_boundary, xs, unitnormalvec);
      return;
    }
    
    
    //! calculate element normal in gausspoint (r,s)
    inline void computeNormalToBoundaryElement(
        const DRT::Element*      surfaceElement,    ///< surface element
        const BlitzMat&          xs,                ///< boundary element nodal coordinates in (x,y,z) space
        const BlitzVec2&         xsi,               ///< coordinates in (r,s) space, where the normal shall be calculated,
        BlitzVec3&               unitnormalvec
            )
    {
      switch (surfaceElement->Shape())
      {
          case DRT::Element::quad4:
              computeNormalToBoundaryElement<DRT::Element::quad4>(xs, xsi, unitnormalvec);
              break;
          case DRT::Element::quad8:
              computeNormalToBoundaryElement<DRT::Element::quad8>(xs, xsi, unitnormalvec);
              break;
          case DRT::Element::quad9:
              computeNormalToBoundaryElement<DRT::Element::quad9>(xs, xsi, unitnormalvec);
              break;
          case DRT::Element::tri3:
              computeNormalToBoundaryElement<DRT::Element::tri3>(xs, xsi, unitnormalvec);
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

    \param surfaceElement       (in)       : surface element
    \param physCoord            (out)      : physical coordinates      
    \param eleCoord             (in)       : element coordinates
    \param normal               (in)       : normal point from a given point to the surface point
    \param distance             (in)       : length of the normal
    */         
    bool searchForNearestPointOnSurface(
            const DRT::Element*                     surfaceElement,
            const BlitzMat&                         xyze_cutter,
            const BlitzVec3&                        physCoord,
            BlitzVec2&                              eleCoord,
            BlitzVec3&                              normal,
            double&                                 distance
            );
    
    
    /*!
       \brief checks if a position in current coordinates lies within a certain surfaceElement           
          The nonlinear system of equation is solved with help of the Newton-method.
    
       \param surfaceElement   (in)        : surfaceElement 
       \param physCoord           	(in)        : node in physical coordinates (x, y, z)
       \return									  node in element coordinates (r, s)
    */  
    void CurrentToSurfaceElementCoordinates(
        const DRT::Element*      surfaceElement,
        const BlitzMat&          xyze_surfaceElement,
        const BlitzVec3&         physCoord,
        BlitzVec2&               xsi
        );
    
    
    /*!
    \brief checks if a position in element coordinates lies within a certain Element parameter space

    \return  true if point lies on the surface element, false otherwise
    */
    template <class V>
    bool checkPositionWithinElementParameterSpace(
            const V&                               eleCoord, ///< position in element coordinates
            const DRT::Element::DiscretizationType distype   ///< shape of the surface element
            )
    {
        if (distype != DRT::Element::line2 and
                distype != DRT::Element::line3 and
                distype != DRT::Element::quad4 and 
                distype != DRT::Element::quad8 and 
                distype != DRT::Element::quad9 and
                distype != DRT::Element::hex8 and
                distype != DRT::Element::hex20 and 
                distype != DRT::Element::hex27)
            dserror("function only defined for rectangular element types at the moment");
        
        bool nodeWithinElement = true;
        
        // loop over r and s (local coordinates)
        for(int i=0; i<DRT::UTILS::getDimension(distype); i++)
            if( (fabs(eleCoord(i))-1.0) > TOL7)     
            {    
                nodeWithinElement = false;
                break;
            }
        
        return nodeWithinElement;
    }
    
    
    /*!
    \brief Computes a rough overestimating extended  
           axis-aligned bounding box for an element (XAABB)

    \param  element          (in)      : element
    \return extended axis-aligned bounding box  (XAABB) for an element
    */
    BlitzMat3x2 computeFastXAABB( 
            const DRT::Element*   element,
            const BlitzMat&       xyze,
            const EleGeoType      eleGeoType      
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
    \brief checks an element is Cartesian, linear or higherorder

    \param element        (in)         : XAABB of the cutting element
    \param xyze_element   (in)         : XAABB of the xfem element
    \param eleGeoType     (in)         : element geometric type CARTESIAN LINEAR or HIGHERORDER
    \return true if the XAABB's intersect or false otherwise
    */
    void checkGeoType(
             DRT::Element*                element,
             const BlitzMat               xyze_element,
             EleGeoType&                  eleGeoType
             );
        
}


#endif  //INTERSECTION_MATH_H_
#endif  // #ifdef CCADISCRET


