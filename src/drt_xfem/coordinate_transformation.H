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
      const GEO::DomainIntCell& cell,
      const PosEtaDomain&  pos_eta_domain,
      PosXiDomain&  pos_xsi_domain
      )
  {   
    // get cell node coordinates
    const BlitzMat* xyze_cell(cell.NodalPosXiDomainBlitz());
    
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
      const PosEtaDomain&  pos_eta_domain,
      PosXiDomain&  pos_xsi_domain
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
      const GEO::DomainIntCell& cell,
      const GEO::PosEtaDomain&  pos_eta
      )
  {   
    // get node coordinates
    const BlitzMat* xyze_cell(cell.NodalPosXiDomainBlitz());
    
    const int numnode = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get shape function derivative
    //LINALG::SerialDenseMatrix deriv(DRT::UTILS::getNumberOfElementNodes(cell.Shape()));
    static LINALG::SerialDenseMatrix deriv(3,27);
    DRT::UTILS::shape_function_3D_deriv1(deriv,pos_eta(0),pos_eta(1),pos_eta(2),cell.Shape());
    // get jacobian matrix ...
    LINALG::FixedSizeSerialDenseMatrix<3,3> xjm;
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
      const GEO::PosEtaDomain&
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
      const GEO::PosEtaBoundary&  pos_eta,
      GEO::PosXiBoundary&         pos_xsi_boundary
  )
  {   
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get cell node coordinates
    const BlitzMat* xyze_cell(cell.NodalPosXiBoundaryBlitz());
    
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
      const GEO::PosEtaBoundary&   pos_eta,
      GEO::PosXiDomain&            pos_xsi_domain
  )
  {   
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get cell node coordinates in xi_domain
    const BlitzMat* xyze_cell(cell.NodalPosXiDomainBlitz());
    
    // create shape function vectors
    LINALG::SerialDenseVector funct(numnodecell);
    DRT::UTILS::shape_function_2D(funct, pos_eta(0),pos_eta(1),cell.Shape());
    
    // translate position into from cell coordinates to element coordinates
    //const LINALG::SerialDenseVector pos_xsi_boundary(blitz::sum(funct(j)*xyze_cell(i,j),j));
    if (numnodecell != 3) dserror("numnode of integrationcell is not 3");
    BLITZTINY::MV_product<3,3>(*xyze_cell, funct, pos_xsi_domain);
    
    return;
  }
  
  
  
  //! get determinante between \f$\eta^{\Gamma (e)}\f$ and \f$\xi^{\Gamma (e)}\f$ space
  inline double detEtaBToXiB(
      const GEO::BoundaryIntCell& cell,
      const GEO::PosEtaBoundary&  pos_eta)
  {   
    // get node coordinates
    const BlitzMat* xiBe_cell(cell.NodalPosXiBoundaryBlitz());
    
    const int numnode = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // get shape function derivative
    LINALG::SerialDenseMatrix deriv(2,numnode);
    DRT::UTILS::shape_function_2D_deriv1(deriv, pos_eta(0), pos_eta(1), cell.Shape());
    
    // get jacobian matrix ...
    LINALG::FixedSizeSerialDenseMatrix<2,2> xjm;
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
  
  //! Calculate inverse of a 3x3 matrix
  template<class M1, class M2>
  inline void Inverse3x3(
      const M1&     mat,      ///< 3x3 matrix
      const double  det,      ///< determinante of mat
      M2&           inverse   ///< inverse 3x3 matrix 
      )
  {
    inverse(0,0) = (  mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2))/det;
    inverse(1,0) = (- mat(1,0)*mat(2,2) + mat(2,0)*mat(1,2))/det;
    inverse(2,0) = (  mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1))/det;
    inverse(0,1) = (- mat(0,1)*mat(2,2) + mat(2,1)*mat(0,2))/det;
    inverse(1,1) = (  mat(0,0)*mat(2,2) - mat(2,0)*mat(0,2))/det;
    inverse(2,1) = (- mat(0,0)*mat(2,1) + mat(2,0)*mat(0,1))/det;
    inverse(0,2) = (  mat(0,1)*mat(1,2) - mat(1,1)*mat(0,2))/det;
    inverse(1,2) = (- mat(0,0)*mat(1,2) + mat(1,0)*mat(0,2))/det;
    inverse(2,2) = (  mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1))/det;
    return;
  }
  
}


#endif

#endif
