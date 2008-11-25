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
#include "xfem_enums.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

//! geometry related things
namespace GEO
{
  
  //! map position from eta to xi space
  template<XFEM::AssemblyType ASSTYPE>
  inline void mapEtaToXi3D(
      const GEO::DomainIntCell&   cell,
      const LINALG::Matrix<3,1>&  pos_eta_domain,
      LINALG::Matrix<3,1>&        pos_xsi_domain
      )
  { 
    const int maxnodecell = 27;
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    if (numnodecell > maxnodecell)
      dserror("increase maxnodecell");
    
    // get cell node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(cell.NodalPosXiDomain());
    
    // create shape function vectors
    static LINALG::Matrix<maxnodecell,1> funct;
    DRT::UTILS::shape_function_3D(funct,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    pos_xsi_domain.Clear();
    for (int j = 0; j < numnodecell; ++j)
      for (int isd = 0; isd < 3; ++isd)
        pos_xsi_domain(isd) += funct(j)*xyze_cell(isd,j);

    return;
  }
  //! map position from eta to xi space
  template<> 
  inline void mapEtaToXi3D<XFEM::standard_assembly>(
      const GEO::DomainIntCell&,
      const LINALG::Matrix<3,1>&  pos_eta_domain,
      LINALG::Matrix<3,1>&        pos_xsi_domain
      )
  {   
    pos_xsi_domain = pos_eta_domain;
    return;
  }  
  
  
  //! get determinante between eta and xi space
  template<XFEM::AssemblyType ASSTYPE>
  inline double detEtaToXi3D(
      const GEO::DomainIntCell&   cell,
      const LINALG::Matrix<3,1>&  pos_eta_domain
      )
  { 
    const int maxnodecell = 27;
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    if (numnodecell > maxnodecell)
      dserror("increase maxnodecell");
    
    // get node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(cell.NodalPosXiDomain());
    
    // get shape function derivative
    static LINALG::Matrix<3,maxnodecell> deriv;
    DRT::UTILS::shape_function_3D_deriv1(deriv,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),cell.Shape());
    // get jacobian matrix ...
    static LINALG::Matrix<3,3> xjm;
    xjm.Clear();
    for (int inode = 0; inode < numnodecell; ++inode)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          xjm(i,j) += deriv(i,inode)*xyze_cell(j,inode);

    // ...and return its determinante
    return xjm.Determinant();
  }
  //! get determinante between eta and xi space
  template<>
  inline double detEtaToXi3D<XFEM::standard_assembly>(
      const GEO::DomainIntCell&,
      const LINALG::Matrix<3,1>&
      )
  {   
    return 1.0;
  }
  
  /////////////////////////////////////////////////////
  // 2d
  //
  //! map position from eta^boundary to xi^boundary space
  inline void mapEtaBToXiB(
      const GEO::BoundaryIntCell& cell,
      const LINALG::Matrix<2,1>&  pos_eta_boundary,
      LINALG::Matrix<2,1>&        pos_xsi_boundary
  )
  {   
    const int maxnodecell = 9;
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    if (numnodecell > maxnodecell)
      dserror("increase maxnodecell");
    
    // get cell node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(cell.NodalPosXiBoundary());
    
    // create shape function vectors
    static LINALG::Matrix<maxnodecell,1> funct;
    DRT::UTILS::shape_function_2D(funct,pos_eta_boundary(0),pos_eta_boundary(1),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    pos_xsi_boundary.Clear();
    for (int j = 0; j < numnodecell; ++j)
      for (int i = 0; i < 2; ++i)
        pos_xsi_boundary(i) += xyze_cell(i,j)*funct(j);

    return;
  }
  
  //! map position from eta^boundary to xi^domain space
  inline void mapEtaBToXiD(
      const GEO::BoundaryIntCell&  cell,
      const LINALG::Matrix<2,1>&   pos_eta_boundary,
      LINALG::Matrix<3,1>&         pos_xsi_domain
  )
  {   
    const int maxnodecell = 9;
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    if (numnodecell > maxnodecell)
      dserror("increase maxnodecell");
        
    // get cell node coordinates in xi_domain
    const LINALG::SerialDenseMatrix& xyze_cell(cell.NodalPosXiDomain());
    
    // create shape function vectors
    static LINALG::Matrix<maxnodecell,1> funct;
    DRT::UTILS::shape_function_2D(funct, pos_eta_boundary(0),pos_eta_boundary(1),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    pos_xsi_domain.Clear();
    for (int j = 0; j < numnodecell; ++j)
      for (int i = 0; i < 3; ++i)
        pos_xsi_domain(i) += xyze_cell(i,j)*funct(j);
    
    return;
  }
  
  
  
  //! get determinante between \f$\eta^{\Gamma (e)}\f$ and \f$\xi^{\Gamma (e)}\f$ space
  inline double detEtaBToXiB(
      const GEO::BoundaryIntCell& cell,
      const LINALG::Matrix<2,1>&  pos_eta_boundary)
  {   
    const int maxnodecell = 9;
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    if (numnodecell > maxnodecell)
      dserror("increase maxnode");
    
    // get node coordinates
    const LINALG::SerialDenseMatrix& xiBe_cell(cell.NodalPosXiBoundary());
    
    // get shape function derivative
    static LINALG::Matrix<2,maxnodecell> deriv;
    DRT::UTILS::shape_function_2D_deriv1(deriv, pos_eta_boundary(0), pos_eta_boundary(1), cell.Shape());
    
    // get jacobian matrix ...
    static LINALG::Matrix<2,2> xjm;
    xjm.Clear();
    for (int k = 0; k < numnodecell; ++k)
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          xjm(i,j) += deriv(i, k)*xiBe_cell(j, k);
    
    // ...and return its determinante
    return xjm.Determinant();
  }
}


#endif

#endif
