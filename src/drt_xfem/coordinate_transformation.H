/*----------------------------------------------------------------------*/
/*!
\file coordinate_transformation.H

\brief routines doing coordinate transformation between various coordinate systems
       during integration by subtriangulation

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef COORDINATE_TRANSFORMATION_H
#define COORDINATE_TRANSFORMATION_H

#include <math.h>
#include <blitz/array.h>
#include "../drt_xfem/xfem.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_utils.H"


namespace XFEM
{
  
  //! map position from eta to xi space
  template<XFEM::AssemblyType ASSTYPE>
  XFEM::PosXiDomain mapEtaToXi3D(
          const XFEM::DomainIntCell& cell,
          const XFEM::PosEtaDomain&  pos_eta)
  {   
      // get cell node coordinates
      const BlitzMat xyze_cell(cell.NodalPosXiDomainBlitz());

      // create shape function vectors 
      const BlitzVec funct(DRT::UTILS::shape_function_3D(pos_eta(0),pos_eta(1),pos_eta(2),cell.Shape()));

      // init blitz indices
      blitz::firstIndex i;    // Placeholder for the first index
      blitz::secondIndex j;   // Placeholder for the second index
      
      // translate position into from cell coordinates to element coordinates
      const BlitzVec pos_xsi_domain(blitz::sum(funct(j)*xyze_cell(i,j),j));
      return pos_xsi_domain;
  }
  //! map position from eta to xi space
  template<> 
  XFEM::PosXiDomain mapEtaToXi3D<XFEM::standard_assembly>(
          const XFEM::DomainIntCell& cell,
          const XFEM::PosEtaDomain&  pos_eta)
  {   
      return pos_eta;
  }

  
  
  //! get determinante between eta and xi space
  template<XFEM::AssemblyType ASSTYPE>
  double detEtaToXi3D(
          const XFEM::DomainIntCell& cell,
          const XFEM::PosEtaDomain&  pos_eta)
  {   
      // get node coordinates
      const BlitzMat xyze_cell(cell.NodalPosXiDomainBlitz());

      // init blitz indices
      blitz::firstIndex i;    // Placeholder for the first index
      blitz::secondIndex j;   // Placeholder for the second index
      blitz::thirdIndex k;    // Placeholder for the third index
  
      // get shape function derivative
      const BlitzMat deriv(DRT::UTILS::shape_function_3D_deriv1(pos_eta(0),pos_eta(1),pos_eta(2),cell.Shape()));
      // get jacobian matrix ...
      const BlitzMat xjm(blitz::sum(deriv(i,k)*xyze_cell(j,k),k));
      // ...and its determinante
      const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                         xjm(0,1)*xjm(1,2)*xjm(2,0)+
                         xjm(0,2)*xjm(1,0)*xjm(2,1)-
                         xjm(0,2)*xjm(1,1)*xjm(2,0)-
                         xjm(0,0)*xjm(1,2)*xjm(2,1)-
                         xjm(0,1)*xjm(1,0)*xjm(2,2);

      return det;
  }
  //! get determinante between eta and xi space
  template<>
  double detEtaToXi3D<XFEM::standard_assembly>(
          const XFEM::DomainIntCell& cell,
          const XFEM::PosEtaDomain&  pos_eta)
  {   
      return 1.0;
  }
  
  
  /////////////////////////////////////////////////////
  // 2d
  //
  //! map position from eta^\boundary to xi^\boundary space
  XFEM::PosXiBoundary mapEtaBToXiB(
          const XFEM::BoundaryIntCell& cell,
          const XFEM::PosEtaBoundary&  pos_eta)
  {   
      // get cell node coordinates
      const BlitzMat xyze_cell(cell.NodalPosXiBoundaryBlitz());

      // create shape function vectors 
      const BlitzVec funct(DRT::UTILS::shape_function_2D(pos_eta(0),pos_eta(1),cell.Shape()));

      // init blitz indices
      blitz::firstIndex i;    // Placeholder for the first index
      blitz::secondIndex j;   // Placeholder for the second index
      
      // translate position into from cell coordinates to element coordinates
      const BlitzVec pos_xsi_boundary(blitz::sum(funct(j)*xyze_cell(i,j),j));
      return pos_xsi_boundary;
  }
  
  //! map position from eta^\boundary to xi^\domain space
  XFEM::PosXiDomain mapEtaBToXiD(
          const XFEM::BoundaryIntCell& cell,
          const XFEM::PosEtaBoundary&  pos_eta)
  {   
      // get cell node coordinates in xi_domain
      const BlitzMat xyze_cell(cell.NodalPosXiDomainBlitz());

      // create shape function vectors 
      const BlitzVec funct(DRT::UTILS::shape_function_2D(pos_eta(0),pos_eta(1),cell.Shape()));

      // init blitz indices
      blitz::firstIndex i;    // Placeholder for the first index
      blitz::secondIndex j;   // Placeholder for the second index
      
      // translate position into from cell coordinates to element coordinates
      const BlitzVec pos_xsi_boundary(blitz::sum(funct(j)*xyze_cell(i,j),j));
      return pos_xsi_boundary;
  }

  
  
//! get determinante between eta^\boundary and xi^\boundary space
double detEtaBToXiB(
        const XFEM::BoundaryIntCell& cell,
        const XFEM::PosEtaDomain&  pos_eta)
{   
    // get node coordinates
    const BlitzMat xyze_cell(cell.NodalPosXiBoundaryBlitz());

    // init blitz indices
    blitz::firstIndex i; // Placeholder for the first index
    blitz::secondIndex j; // Placeholder for the second index
    blitz::thirdIndex k; // Placeholder for the third index

    // get shape function derivative
    const BlitzMat deriv(DRT::UTILS::shape_function_2D_deriv1(pos_eta(0), pos_eta(1), cell.Shape()));
    // get jacobian matrix ...
    const BlitzMat xjm(blitz::sum(deriv(i, k)*xyze_cell(j, k), k));
    // ...and its determinante
    const double det = xjm(0, 0)*xjm(1, 1) - xjm(0, 1)*xjm(1, 0);

    return det;
}

}



#endif

#endif
