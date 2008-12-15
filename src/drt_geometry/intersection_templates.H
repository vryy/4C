/*!----------------------------------------------------------------------
\file intersection_templates.H

\brief collection of templated methods in intersection computations

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
#ifndef INTERSECTION_TEMPLATES_H_
#define INTERSECTION_TEMPLATES_H_


#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_geometry/intersection_math.H"


namespace GEO
{
	
    /*!
    \brief tests if one entry is NaN according to IEEE NaN is not comparable to itself
     */
    template <int length>
	  static inline bool isNaN(
	      const LINALG::Matrix<length,1>&		xsi)
	  {
	    for(int i = 0; i < length; i++)
	      if( !( xsi(i)==xsi(i)) ) // is NaN
	        return true;
				
	    return false; 
	  }

    /*!
    \brief checks if intersection point lies within a surface and line element
    */
    template<	DRT::Element::DiscretizationType surftype,
    			class V1, 
    			class V2, 
    			class V3>
    static inline bool intersectionPointInElements(
        const V1&        							xsi,
        const V2&        							upLimit,
        const V3&        							loLimit) 
    {	
	    switch (surftype)
	    {
	      case DRT::Element::quad4: case DRT::Element::quad8: case DRT::Element::quad9:
	      {  
	        if( (xsi(0) > (upLimit(0)+GEO::TOL7)) || (xsi(1) > (upLimit(1)+GEO::TOL7)) || (xsi(2) > (upLimit(2)+GEO::TOL7))  ||
	            (xsi(0) < (loLimit(0)-GEO::TOL7)) || (xsi(1) < (loLimit(1)-GEO::TOL7)) || (xsi(2) < (loLimit(2)-GEO::TOL7)))
	          return false;
	        else
	          return true;
	      }
	      case DRT::Element::tri3: case DRT::Element::tri6:
	      {
	        if((xsi(0) > (upLimit(0)+GEO::TOL7)) || (xsi(1) > (upLimit(1)+GEO::TOL7)) || (xsi(2) > (upLimit(2)+GEO::TOL7))  ||
	           (xsi(0) < (loLimit(0)-GEO::TOL7)) || (xsi(1) < (loLimit(1)-GEO::TOL7)) || (xsi(2) < (loLimit(2)-GEO::TOL7))  ||
	           (xsi(1) > (1.0 - xsi(0) + GEO::TOL7) ) )   // checks if in triangle 
	          return false;
	        else
	          return true;
	      }
	      default:
	        dserror("surface type not implemented");
	        return false;
	      };	
    }
    
       
    
    /*!
    \brief updates the systemmatrix at the corresponding element coordinates
           for the computation of curve surface intersections
    */
    template< DRT::Element::DiscretizationType surftype,
              DRT::Element::DiscretizationType linetype,
              class M1,
              class M2,
              class V>
      static inline void updateAForCSI(
          LINALG::Matrix<3,3>&              A,                   ///< system matrix
          const V&                          xsi,                 ///< vector of element coordinates (r,s,t)
          const M1&                         xyze_surfaceElement, ///< nodal positions of surface element
          const M2&                         xyze_lineElement     ///< nodal positions of line element
    )
    {
      const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
      const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;
    
      A.Clear();
      static LINALG::Matrix<2,numNodesSurface> surfaceDeriv1;
      DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surftype);
      
      for(int inode=0; inode<numNodesSurface; inode++)
        for(int isd=0; isd<3; isd++)
        {
          A(isd,0) += xyze_surfaceElement(isd,inode) * surfaceDeriv1(0,inode);
          A(isd,1) += xyze_surfaceElement(isd,inode) * surfaceDeriv1(1,inode);
        }
    
    
      static LINALG::Matrix<2, numNodesLine> lineDeriv1;
      DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), linetype);
      
      for(int inode=0; inode<numNodesLine; inode++)
        for(int isd=0; isd<3; isd++)
          A(isd,2) -= xyze_lineElement(isd,inode) * lineDeriv1(0,inode);
    }
    
    
    
    /*!
    \brief updates the rhs at the corresponding element coordinates
           for the computation of curve surface intersections
    */
    template<DRT::Element::DiscretizationType surftype,
             DRT::Element::DiscretizationType linetype,
             class M1,
             class M2,
             class V>
    static inline void updateRHSForCSI(
        LINALG::Matrix<3,1>&              b,                   ///< right-hand-side
        const V&                          xsi,                 ///< vector of element coordinates (r,s,t)
        const M1&                         xyze_surfaceElement, ///< nodal positions of surface element
        const M2&                         xyze_lineElement     ///< nodal positions of line element
        )
    {
      const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
      const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;
    
      b.Clear();
    
      static LINALG::Matrix<numNodesSurface,1> surfaceFunct;
      DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surftype);
      for(int i=0; i<numNodesSurface; i++)
      {
        for(int dim=0; dim<3; dim++)
          b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);
      }
    
      static LINALG::Matrix<numNodesLine,1> lineFunct;
      DRT::UTILS::shape_function_1D(lineFunct, xsi(2), linetype);
      for(int i=0; i<numNodesLine; i++)
      {
        for(int dim=0; dim<3; dim++)
          b(dim) += xyze_lineElement(dim,i) * lineFunct(i);
      }
    }
    
    
    /*!
    \brief solves a singular system of equations
    \return true if resulting system is singular , false otherwise
    */
    template<DRT::Element::DiscretizationType surftype,
             DRT::Element::DiscretizationType linetype,
             class M1,
             class M2,
             class V>
    bool computeSingularCSI(
        V&                      xsi,                 ///< vector of element coordinates (r,s,t)
        const M1&               xyze_surfaceElement, ///< nodal positions of surface element
        const M2&               xyze_lineElement     ///< nodal positions of line element
        )
    {
      int iter = 0;
      const int maxiter = 5;
      double residual = 1.0;
      static LINALG::Matrix<3,3> A;
      static LINALG::Matrix<3,1> b;
      static LINALG::Matrix<3,1> dx;
    
      updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);
    
      while(residual > GEO::TOL13)
      {
        updateAForCSI<surftype,linetype>( A, xsi, xyze_surfaceElement, xyze_lineElement);
    
        dx.Clear();
        if(GEO::solveLinearSystemWithSVD<3>(A, b, dx, GEO::TOL14))
        {
          xsi += dx;
          return false;
        }
    
        xsi += dx;
        updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);
        residual = b.Norm2();
        iter++;
    
        if(  iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || isNaN<3>(xsi)  )
        	return true;
      }
      return true;
    }
    
    
    /*!
    \brief computes an interseticon point between a curve and a surface - templated part
    
        The nonlinear system of equation is solved with help of the Newton-method.
    
    \param xyze_surfaceElement      (in)    : nodal coordinates of surface element
    \param xyze_lineElement         (in)    : nodal coordinates of line element
    \param xsi                      (in/out): starting value/vector of element coordinates
    \param upLimit                  (in)    : upper search interval boundary
    \param loLimit                  (in)    : lower search interval boundary
    \param doSVD                    (in)    : compute SVD if not Cartesian surface element and linear line elments present
    return true if an intersection point was found, otherwise false
    */
    template<DRT::Element::DiscretizationType surftype,
             DRT::Element::DiscretizationType linetype,
             class M1,
             class M2,
             class V1, class V2, class V3>
    bool computeCurveSurfaceIntersectionT(
        const M1&                     	xyze_surfaceElement,
        const M2&                   	xyze_lineElement,
        V1&                           	xsi,
        const V2&        				upLimit,
        const V3&        				loLimit,
        const bool                     	doSVD
    )
    {
    
      int iter = 0;
      const int maxiter = 20;
      double residual = 1.0;
      static LINALG::Matrix<3,3> A;
      static LINALG::Matrix<3,1> b;
      static LINALG::Matrix<3,1> dx;
    
      updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);
    
      while(residual > GEO::TOL13)
      {
        updateAForCSI<surftype,linetype>( A, xsi, xyze_surfaceElement, xyze_lineElement);
        const bool singular = !GEO::gaussElimination<true,3>(A, b, dx, GEO::TOL14);
    
        if(singular && !doSVD)
          return false;
        else if(singular && doSVD)
        {
          // xsi contains the solution
          if(computeSingularCSI<surftype,linetype>(xsi, xyze_surfaceElement, xyze_lineElement))
            return false;
        	
          dx.Clear();
        }
    
        xsi += dx;
        updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);
        residual = b.Norm2();
        iter++;
    
        // printf("xsi = %20.16f   %20.16f   %20.16f  iter = %d res = %20.16f\n", xsi(0), xsi(1), xsi(2), iter, residual  );
        // has to to be 8 , otherwise not a number is reached               // detect not a number according to IEEE NaN is not comparable to itself
        if(	iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || isNaN<3>(xsi)    )
          return false;
      }    
      //printf("xsi = %f   %f   %f  iter = %d res = %20.16f\n", xsi(0), xsi(1), xsi(2), iter, residual  );
      return intersectionPointInElements<surftype>(xsi, upLimit, loLimit);
    }



    /*!
    \brief computes an interseticon point between a curve and a surface
    
        The nonlinear system of equation is solved with help of the Newton-method.
    
    \param surfaceElement           (in)    : surface element
    \param xyze_surfaceElement      (in)    : nodal coordinates of surface element
    \param lineElement              (in)    : line element
    \param xyze_lineElement         (in)    : nodal coordinates of line element
    \param xsi                      (in/out): starting value/vector of element coordinates
    \param upLimit                  (in)    : upper search interval boundary
    \param loLimit                  (in)    : lower search interval boundary
    \param doSVD                    (in)    : compute SVD if not Cartesian surface element and linear line elments present
    return true if an intersection point was found, otherwise false
    */
    template<class M1, class M2, class V1, class V2, class V3>
    bool computeCurveSurfaceIntersection(
        const DRT::Element*               surfaceElement,
        const M1&                         xyze_surfaceElement,
        const DRT::Element*               lineElement,
        const M2&                         xyze_lineElement,
        const V2&        				  upLimit,
        const V3&        				  loLimit,
        V1&                               xsi,
        const bool                        doSVD
    ) 
    {
      if (lineElement->Shape() == DRT::Element::line2)
      {
        switch (surfaceElement->Shape())
        {
        case DRT::Element::quad4:
          return computeCurveSurfaceIntersectionT<DRT::Element::quad4,DRT::Element::line2>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
        case DRT::Element::quad8:
          return computeCurveSurfaceIntersectionT<DRT::Element::quad8,DRT::Element::line2>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
        case DRT::Element::quad9:
          return computeCurveSurfaceIntersectionT<DRT::Element::quad9,DRT::Element::line2>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
        case DRT::Element::tri3:
          return computeCurveSurfaceIntersectionT<DRT::Element::tri3 ,DRT::Element::line2>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
        case DRT::Element::tri6:
          return computeCurveSurfaceIntersectionT<DRT::Element::tri6 ,DRT::Element::line2>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
        default:
          dserror("template not instatiated yet");
        return false;
      };
    }
    else if (lineElement->Shape() == DRT::Element::line3)
    {
      switch (surfaceElement->Shape())
      {
      case DRT::Element::quad4:
        return computeCurveSurfaceIntersectionT<DRT::Element::quad4,DRT::Element::line3>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
      case DRT::Element::quad8:
        return computeCurveSurfaceIntersectionT<DRT::Element::quad8,DRT::Element::line3>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
      case DRT::Element::quad9:
        return computeCurveSurfaceIntersectionT<DRT::Element::quad9,DRT::Element::line3>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
      case DRT::Element::tri3:
        return computeCurveSurfaceIntersectionT<DRT::Element::tri3 ,DRT::Element::line3>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
      case DRT::Element::tri6:
        return computeCurveSurfaceIntersectionT<DRT::Element::tri6 ,DRT::Element::line3>(xyze_surfaceElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD);
      default:
        dserror("template not instatiated yet");
          return false;
        };
      }
      return true;
    }
    
    
    
    
    /*!
    \brief checks if intersection point lies within a surface element
    */
    template<	DRT::Element::DiscretizationType surftype,
    			class V >
    static inline bool intersectionPointInElementsRCI(
        const V&        							xsi) 
    {	
      switch (surftype)
	    {
	      case DRT::Element::quad4: case DRT::Element::quad8: case DRT::Element::quad9:
	      {  
	        if((xsi(0) > ( 1.0 + GEO::TOL7)) || (xsi(1) > ( 1.0 + GEO::TOL7)) ||
	           (xsi(0) < (-1.0 - GEO::TOL7)) || (xsi(1) < (-1.0 - GEO::TOL7))   )
	          return false;
	        else
	          return true;
	      }
	      case DRT::Element::tri3: case DRT::Element::tri6:
	      {
	        // checks if in triangle 
	        if((xsi(0) > ( 1.0 + GEO::TOL7)) || (xsi(1) > (1.0 - xsi(0) + GEO::TOL7) )  ||  
	           (xsi(0) < ((-1.0)*GEO::TOL7)) || (xsi(1) < ((-1.0)*GEO::TOL7)) )
	          return false;
	        else
	          return true;
	      }
	      default:
	        dserror("surface type not implemented");
	
	    };	
	    return false;
    }
    
    
    
    
    /*!
    \brief  updates the system matrix for the intersection computation
            of an average normal and the corresponding cutter element
            to lift a point of the linearized interface onto the curved interface
  
    \param A                    (out)       :   system matrix
    \param xsi                  (in)        :   intersection point
    \param normal               (in)        :   nodes of the intersection normal
    \param surfaceElement       (in)        :   surface element
    \param xyze_surfaceElement  (in)        :   nodal coordinates of surface element
    \param onBoundary           (in)        :   true, if intersection normal lies on the
                                                element boundary
    */
    template< DRT::Element::DiscretizationType surftype,
              class M,
              class V>
    static inline void updateAForRCINormal(
        LINALG::Matrix<3,3>&                        A,
        const V&                  					xsi,
        const vector<LINALG::Matrix<3,1> >&         normal,
        const M&            						xyze_surfaceElement,
        const bool                                  onBoundary
        ) 
    {
      const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;

      A.Clear();
      static LINALG::Matrix<2,numNodesSurface> surfaceDeriv1;
      DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surftype);

      if(!onBoundary)
      {
        for(int i=0; i<numNodesSurface; i++)
          for(int dim=0; dim<3; dim++)
          {
            A(dim,0) += xyze_surfaceElement(dim,i) * surfaceDeriv1(0, i);
            A(dim,1) += xyze_surfaceElement(dim,i) * surfaceDeriv1(1, i);
          }

        for(int dim=0; dim<3; dim++)
          A(dim,2) -= 0.5*( -normal[0](dim) + normal[1](dim));
      }
      else
      {
        const int numNodesLine = 3;
        static LINALG::Matrix<1,3> lineDeriv1;
        // LINALG::SerialDenseMatrix lineDeriv1(1,numNodesLine);
        DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), DRT::Element::line3);

        for(int i=0; i<numNodesSurface; i++)
          for(int dim=0; dim<3; dim++)
          {
            A(dim,0) += xyze_surfaceElement(dim,i) * surfaceDeriv1(0,i);
            A(dim,1) += xyze_surfaceElement(dim,i) * surfaceDeriv1(1,i);
          }

        for(int i = 0; i < numNodesLine; i++)
          for(int dim=0; dim<3; dim++)
          {
            int index = i;
            if(i > 1)   index = 4;
            A(dim,2) -= normal[index](dim) * lineDeriv1(0,i);
          }
      }
    }



    /*!
    \brief  updates the right-hand-side for the intersection computation
            of an average normal and the corresponding cutter element
            to lift a point of the linearized interface onto the curved interface
  
    \param b                    (out)       :   right-hand-side
    \param xsi                  (in)        :   intersection point
    \param normal               (in)        :   nodes of the intersection normal
    \param xyze_surfaceElement  (in)        :   nodal coordinates of surface element
    \param onBoundary           (in)        :   true, if intersection normal lies on the
                                                element boundary
    */
  	template< DRT::Element::DiscretizationType surftype,
  	          class M,
  	          class V>
  	static inline void updateRHSForRCINormal(
  	    LINALG::Matrix<3,1>&                	b,
  	    const V&            					xsi,
  	    const vector<LINALG::Matrix<3,1> >&   	normal,
  	    const M&      							xyze_surfaceElement,
  	    const bool                            	onBoundary) 
  	{
  	  const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
  	  static LINALG::Matrix<numNodesSurface,1> surfaceFunct;
  	  DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surftype);

  	  b.Clear();
  	  if(!onBoundary)
  	  {
  	    for(int i=0; i<numNodesSurface; i++)
  	      for(int dim=0; dim<3; dim++)
  	        b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);


  	    for(int dim=0; dim<3; dim++)
  	      b(dim) += 0.5*(normal[0](dim)*(1.0 - xsi(2)) + normal[1](dim)*(1.0 + xsi(2)));
  	  }
  	  else
  	  {
  	    const int numNodesLine = 3;
  	    static LINALG::Matrix<3,1> lineFunct;
  	    DRT::UTILS::shape_function_1D(lineFunct, xsi(2), DRT::Element::line3);

  	    for(int i=0; i<numNodesSurface; i++)
  	      for(int dim=0; dim<3; dim++)
  	        b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);



  	    for(int i=0; i<numNodesLine; i++)
  	      for(int dim=0; dim<3; dim++)
  	      {
  	        int index = i;
  	        if(i > 1)   index = 4;

  	        b(dim) += normal[index](dim) * lineFunct(i);
  	      }
  	  }
  	}
	
	
	
    /*!
    \brief solves a singular system of equations
    
  	\return true if resulting system is singular , false otherwise
    */
    template<DRT::Element::DiscretizationType surftype,
             class M,
             class V>
    bool computeSingularRCI(
        V&                      				xsi,                 	///< vector of element coordinates (r,s,t)
        const M&               					xyze_cutterElement, 	///< nodal positions of surface element
        const vector<LINALG::Matrix<3,1> >&   	normal,     			///< nodal positions of line element
        const bool 								onBoundary)
    {
      int iter = 0;
      const int maxiter = 5;
      double residual = 1.0;
      static LINALG::Matrix<3,3> A;
      static LINALG::Matrix<3,1> b;
      static LINALG::Matrix<3,1> dx;

      updateRHSForRCINormal<surftype>( b, xsi, normal, xyze_cutterElement, onBoundary);

      while(residual > GEO::TOL13)
      {
        updateAForRCINormal<surftype>( A, xsi, normal, xyze_cutterElement, onBoundary);

        dx.Clear();
        if(GEO::solveLinearSystemWithSVD<3>(A, b, dx, GEO::TOL14))
        {
          xsi += dx;
          return false;
        }

        xsi += dx;
        updateRHSForRCINormal<surftype>( b, xsi, normal, xyze_cutterElement, onBoundary);
        residual = b.Norm2();
        iter++;

        if(  iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || isNaN<3>(xsi)  )
          return true;
      }
      return true;
    }




	/*----------------------------------------------------------------------*
	 |  RCI:    computes the intersection between a              u.may 08/07|
	 |          line and a surface                                          |
	 *----------------------------------------------------------------------*/
    template< DRT::Element::DiscretizationType surftype,
			        class V,
			        class M>
    bool computeRecoveryNormalT(
            const M&          						xyze_cutterElement,
            const vector<LINALG::Matrix<3,1> >&   	normal,
            V&                      				xsi,
            const bool                             	onBoundary) 
    {
	    int                         iter = 0;
	    const int                   maxiter = 20;
	    double                      residual = 1.0;
	    
	    // TODO remove
	    bool doSVD = true;

	    static LINALG::Matrix<3,3> 	A;
	    static LINALG::Matrix<3,1> 	b;
	    static LINALG::Matrix<3,1> 	dx;

	    updateRHSForRCINormal<surftype>( b, xsi, normal, xyze_cutterElement, onBoundary);
	    // put normal in linalg matrix and updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);

	    while(residual > GEO::TOL13)
	    {
	      updateAForRCINormal<surftype>( A, xsi, normal, xyze_cutterElement, onBoundary);
	      const bool singular = !GEO::gaussElimination<true,3>(A, b, dx, GEO::TOL14);

	      if(singular && !doSVD)
	        return false;
	      else if(singular && doSVD)
	      {
	        cout << "SINGULAR " << endl;
	        if(computeSingularRCI<surftype>(xsi, xyze_cutterElement, normal, onBoundary))
	          return false; 

	        dx.Clear();
	      }

	      xsi += dx;
	      //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx(0), dx(1), dx(2));
	      if(	iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || isNaN<3>(xsi)) 
	        return false;

	      updateRHSForRCINormal<surftype>( b, xsi, normal, xyze_cutterElement, onBoundary);
	      residual = b.Norm2();
	      cout << "residual = " << residual << endl;
	      iter++;   
	    }
	    return intersectionPointInElementsRCI<surftype>(xsi);
	}




    /*!
    \brief  computes the recovery of the curved interface
            for one xfem element after a Contrained Delaunay
            Tetrahedralization.
  
    \param xsi                      (out)   :   intersection point
    \param normal                   (in)    :   nodes of the intersection normal
    \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element
    \param onBoundary               (in)    :   indicates whether the normal lies on the xfem element
                                                boundary or not
    \return true if normal intersects with cutter element or false otherwise
    */
    template<	class M, 
    			    class V>
    bool computeRecoveryNormal(
        const DRT::Element*                   cutterElement,
        const M&                              xyze_cutterElement,
       	const vector<LINALG::Matrix<3,1> >&   normal,
        V&                                    xsi,
        const bool                            onBoundary) 
    {
        switch (cutterElement->Shape())
        {
        case DRT::Element::quad4:
          return computeRecoveryNormalT<DRT::Element::quad4>(xyze_cutterElement, normal, xsi, onBoundary);
        case DRT::Element::quad8:
          return computeRecoveryNormalT<DRT::Element::quad8>(xyze_cutterElement, normal, xsi, onBoundary);
        case DRT::Element::quad9:
          return computeRecoveryNormalT<DRT::Element::quad9>(xyze_cutterElement, normal, xsi, onBoundary);
        case DRT::Element::tri3:
          return computeRecoveryNormalT<DRT::Element::tri3>(xyze_cutterElement, normal, xsi, onBoundary);
        case DRT::Element::tri6:
          return computeRecoveryNormalT<DRT::Element::tri6>(xyze_cutterElement, normal, xsi, onBoundary);
        default:
          dserror("surftype not yet implemented");
        return false;
      };
    
      return true;
    }




    /*!
    \brief  updates the system matrix for the intersection computation
            of an average normal and the corresponding cutter element
            to lift a point of the linearized interface onto the curved interface
  
    \param A                        (out)       :   system matrix
    \param xsi                      (in)        :   intersection point
    \param plane                    (in)        :   nodes of the intersection normal
    \param xyze_lineElement         (in)        :   nodal coordinates of line element
    */
    template< DRT::Element::DiscretizationType linetype,
              class V,
              class M>
    static inline void updateAForRCIPlane(
            LINALG::Matrix<3,3>&                    A,
            const V&                                xsi,
            const vector<LINALG::Matrix<3,1> >&     plane,
            const M&                                xyze_lineElement)
    {
      const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;
      const int numNodesSurface = 4;

      static LINALG::Matrix<2,numNodesSurface> surfaceDeriv;
      DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv, xsi(0),  xsi(1), DRT::Element::quad4);
      static LINALG::Matrix<1,numNodesLine> lineDeriv;
      DRT::UTILS::shape_function_1D_deriv1(lineDeriv, xsi(2), linetype);

      dsassert((int)plane.size() >= numNodesSurface, "plane array has to have size numNodesSurface ( = 4)!");

      A.Clear();
      for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
          A(dim,0) += plane[i](dim) * surfaceDeriv(0,i);
          A(dim,1) += plane[i](dim) * surfaceDeriv(1,i);
        }

      for(int i=0; i<numNodesLine; i++)
        for(int dim=0; dim<3; dim++)
          A(dim,2) -= xyze_lineElement(dim,i) * lineDeriv(0,i);
    }




    /*!
    \brief  updates the right-hand-side for the intersection computation
            of the cutter element edge with a plane through the midpoint
            of the face edge to lift a point of the linearized interface
            onto the curved interface
  
    \param b                        (out)       :   right-hand-side
    \param xsi                      (in)        :   intersection point
    \param plane                    (in)        :   nodes of the intersection normal
    \param xyze_lineElement         (in)        :   nodal coordinates of line element
    */
    template< DRT::Element::DiscretizationType linetype,
              class V,
              class M>
    static inline void updateRHSForRCIPlane(
        LINALG::Matrix<3,1>&                  	b,
        const V&            					xsi,
        const vector<LINALG::Matrix<3,1> >&   	plane,
        const M&            					xyze_lineElement) 
    {
      const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;
      const int numNodesSurface = 4;

      static LINALG::Matrix<numNodesSurface,1> surfaceFunct;
      DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), DRT::Element::quad4 );
      static LINALG::Matrix<numNodesLine,1> lineFunct;
      DRT::UTILS::shape_function_1D(lineFunct, xsi(2), linetype);

      dsassert((int)plane.size() >= numNodesSurface, "plane array has to have size numNodesSurface ( = 4)!");

      b.Clear();
      for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
          b(dim) -= plane[i](dim) * surfaceFunct(i);

      for(int i=0; i<numNodesLine; i++)
        for(int dim=0; dim<3; dim++)
          b(dim) +=  xyze_lineElement(dim,i)  * lineFunct(i);
    }




    /*!
    \brief  computes the recovery of the curved interface
            for one xfem element after a Contrained Delaunay
            Tetrahedralization. 
  
    \param lineIndex                (in)    :   line index
    \param currentcutterpositions   (in)    :   current cutter node positions
    \param xsi                      (out)   :   intersection point
    \param plane                    (in)    :   nodes of the intersection normal
    \param cutterElement            (in)    :   cutter element
    \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element
  
    \return true if plane intersects with the line element or false otherwise
    */
    template< DRT::Element::DiscretizationType linetype,
              class V,
              class M>
    bool computeRecoveryPlaneT(
    		    const M&                                xyze_lineElement,
            const vector<LINALG::Matrix<3,1> >&     plane,
            V&                                      xsi) 
    {
      int                         iter = 0;
      const int                   maxiter = 20;
      double                      residual = 1.0;

      static LINALG::Matrix<3,3> A;
      static LINALG::Matrix<3,1> b;
      static LINALG::Matrix<3,1> dx;

      updateRHSForRCIPlane<linetype>( b, xsi, plane, xyze_lineElement);

      while(residual > GEO::TOL13)
      {
        updateAForRCIPlane<linetype>( A, xsi, plane, xyze_lineElement);

        // check if singular computatio necessary
        if(!gaussElimination<true,3>(A, b, dx, GEO::TOL13))
          return false;

        if(iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || isNaN<3>(xsi)) 
          return false;

        xsi += dx;
        updateRHSForRCIPlane<linetype>( b, xsi, plane, xyze_lineElement);
        residual = b.Norm2();
        iter++;
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
      }

      if( (xsi(2) > ( 1.0 + GEO::TOL7))  || (xsi(2) < (-1.0 - GEO::TOL7)) )     // planes coordinate may be bigger than 1
        return false; 

      //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
      return true;
    }



   /*!
    \brief  computes the recovery of the curved interface
            for one xfem element after a Contrained Delaunay
            Tetrahedralization. 
  
    \param lineIndex                (in)    :   line index
    \param currentcutterpositions   (in)    :   current cutter node positions
    \param xsi                      (out)   :   intersection point
    \param plane                    (in)    :   nodes of the intersection normal
    \param cutterElement            (in)    :   cutter element
    \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element
  
    \return true if plane intersects with the line element or false otherwise
    */
    template<	class M, 
    			    class V>
    bool computeRecoveryPlane(
        const DRT::Element*                 lineElement,
        const M&                            xyze_lineElement,
       	const vector<LINALG::Matrix<3,1> >& plane,
        V&                                  xsi
    ) 
    {
      switch (lineElement->Shape())
      {
      case DRT::Element::line2:
        return computeRecoveryPlaneT<DRT::Element::line2>(xyze_lineElement, plane, xsi);
      case DRT::Element::line3:
        return computeRecoveryPlaneT<DRT::Element::line3>(xyze_lineElement, plane, xsi);
      default:
        dserror("linetype not yet implemented");
        return false;
      };
      return true;
    }


}  /*namespace GEO*/



#endif /*INTERSECTION_TEMPLATES_H_*/
#endif /*CCADISCRET*/





