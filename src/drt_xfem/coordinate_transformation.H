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
#include "../drt_geometry/blitz_tiny_operation.H"
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
    // get cell node coordinates
    const LINALG::SerialDenseMatrix* xyze_cell(cell.NodalPosXiDomain());
    
    // create shape function vectors
    //LINALG::SerialDenseVector funct(DRT::UTILS::getNumberOfElementNodes(cell.Shape()));
    static LINALG::SerialDenseVector funct(27);
    DRT::UTILS::shape_function_3D(funct,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    const int numnodes = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    for (int isd = 0; isd < 3; ++isd)
    {
      pos_xsi_domain(isd) = 0.0;
      for (int j = 0; j < numnodes; ++j)
      {
        pos_xsi_domain(isd) += funct(j)*(*xyze_cell)(isd,j);
      }
    }
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
    pos_xsi_domain(0) = pos_eta_domain(0);
    pos_xsi_domain(1) = pos_eta_domain(1);
    pos_xsi_domain(2) = pos_eta_domain(2);
    return;
  }  
  
  
  //! get determinante between eta and xi space
  template<XFEM::AssemblyType ASSTYPE>
  inline double detEtaToXi3D(
      const GEO::DomainIntCell&   cell,
      const LINALG::Matrix<3,1>&  pos_eta_domain
      )
  {   
    // get node coordinates
    const LINALG::SerialDenseMatrix* xyze_cell(cell.NodalPosXiDomain());
    
    const int numnode = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get shape function derivative
    //LINALG::SerialDenseMatrix deriv(DRT::UTILS::getNumberOfElementNodes(cell.Shape()));
    static LINALG::SerialDenseMatrix deriv(3,27);
    DRT::UTILS::shape_function_3D_deriv1(deriv,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),cell.Shape());
    // get jacobian matrix ...
    LINALG::Matrix<3,3> xjm;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        xjm(i,j) = 0.0;
        for (int inode = 0; inode < numnode; ++inode)
        {
          xjm(i,j) += deriv(i,inode)*(*xyze_cell)(j,inode);
        }
      }
    }

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
      const LINALG::Matrix<2,1>&  pos_eta,
      LINALG::Matrix<2,1>&         pos_xsi_boundary
  )
  {   
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get cell node coordinates
    const LINALG::SerialDenseMatrix* xyze_cell(cell.NodalPosXiBoundary());
    
    // create shape function vectors
    LINALG::SerialDenseVector funct(numnodecell);
    DRT::UTILS::shape_function_2D(funct,pos_eta(0),pos_eta(1),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    //const LINALG::SerialDenseVector pos_xsi_boundary(blitz::sum(funct(j)*xyze_cell(i,j),j));
    if (numnodecell != 3) dserror("numnode of integrationcell is not 3");
    //GEO::PosXiBoundary pos_xsi_boundary;
    BLITZTINY::MV_product<2,3>(*xyze_cell, funct, pos_xsi_boundary);
    return;
  }
  
  //! map position from eta^boundary to xi^domain space
  inline void mapEtaBToXiD(
      const GEO::BoundaryIntCell&  cell,
      const LINALG::Matrix<2,1>&   pos_eta_boundary,
      LINALG::Matrix<3,1>&         pos_xsi_domain
  )
  {   
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get cell node coordinates in xi_domain
    const LINALG::SerialDenseMatrix* xyze_cell(cell.NodalPosXiDomain());
    
    // create shape function vectors
    LINALG::SerialDenseVector funct(numnodecell);
    DRT::UTILS::shape_function_2D(funct, pos_eta_boundary(0),pos_eta_boundary(1),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    //const LINALG::SerialDenseVector pos_xsi_boundary = xyze_cell(i,j)*funct(j);
    if (numnodecell != 3) dserror("numnode of integrationcell is not 3");
    BLITZTINY::MV_product<3,3>(*xyze_cell, funct, pos_xsi_domain);
//    elevec1.GEMV('N', 3, numnodecell, 1.0, xyze_cell.A(), xyze_cell.LDA(), fa.A(), 1.0, elevec1.A());
    
    return;
  }
  
  
  
  //! get determinante between \f$\eta^{\Gamma (e)}\f$ and \f$\xi^{\Gamma (e)}\f$ space
  inline double detEtaBToXiB(
      const GEO::BoundaryIntCell& cell,
      const LINALG::Matrix<2,1>&  pos_eta_boundary)
  {   
    // get node coordinates
    const LINALG::SerialDenseMatrix* xiBe_cell(cell.NodalPosXiBoundary());
    
    const int numnode = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get shape function derivative
    LINALG::SerialDenseMatrix deriv(2,numnode);
    DRT::UTILS::shape_function_2D_deriv1(deriv, pos_eta_boundary(0), pos_eta_boundary(1), cell.Shape());
    
    // get jacobian matrix ...
    LINALG::Matrix<2,2> xjm;
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        xjm(i,j) = 0.0;
        for (int k = 0; k < numnode; ++k)
        {
          xjm(i,j) += deriv(i, k)*(*xiBe_cell)(j, k);
        }
      }
    }
    
    // ...and return its determinante
    return xjm.Determinant();
  }
}


#endif

#endif
