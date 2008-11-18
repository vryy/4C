/*!
\file spacetime_boundary.cpp

\brief integration cell classes for domain and boundary integration

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "spacetime_boundary.H"
#include <string>
#include <sstream>
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


/*!
 * \brief create array with physical coordinates based an local coordinates of a parent element
 */
template<class Cell>
static void ComputePhysicalCoordinates(
        const DRT::Element&  ele,  ///< parent element
        const Cell&          cell, ///< integration cell whose coordinates we'd like to transform
        BlitzMat&            physicalCoordinates
        )
{
    const BlitzMat eleCoord(GEO::InitialPositionArrayBlitz(&ele));
    //DRT::UTILS::fillInitialPositionArray(&ele, eleCoord);
    const LINALG::SerialDenseMatrix* nodalPosXiDomain = cell.NodalPosXiDomain();
    
    const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // return value
    //BlitzMat physicalCoordinates(3, nen_cell);
    physicalCoordinates = 0.0;
    // for each cell node, compute physical position
    const int nen_ele = ele.NumNode();
    LINALG::SerialDenseVector funct(nen_ele);
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        DRT::UTILS::shape_function_3D(funct,
                (*nodalPosXiDomain)(0, inen),
                (*nodalPosXiDomain)(1, inen),
                (*nodalPosXiDomain)(2, inen),
                ele.Shape());

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < nen_ele; ++j)
            {
              physicalCoordinates(i,inen) += eleCoord(i, j) * funct(j);
            }
        }
    };
    return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell(
    const int           bele_id,
    const BlitzMat&     posnp,
    const BlitzMat&     posn
    ) :
      bele_id_(bele_id),
      posnp_(posnp),
      posn_(posn),
      xyzt_(getLinearPositionArray(posnp,posn))
{
    return;
}
    
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell() :
      bele_id_(-1),
      posnp_(),
      posn_(),
      xyzt_()
{
    return;   
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
//XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell(
//    const SpaceTimeBoundaryCell& old
//    ) : 
//      bele_id_(old.bele_id_),
//      posnp_(old.posnp_),
//      posn_(old.posn_),
//      xyzt_(old.xyzt_)
//{
//    return;   
//}
        
BlitzMat XFEM::SpaceTimeBoundaryCell::getLinearPositionArray(
    const BlitzMat&      posnp,                 ///< nodal positions at n+1
    const BlitzMat&      posn                   ///< nodal positions at n
    ) const
{
  BlitzMat xyzt(3,8);
  for (int inode = 0; inode != 4; ++inode) // fill n   position
  {
    for (int isd = 0; isd != 3; ++isd)
    {
      xyzt(isd,inode  ) = posn(isd,inode);
    }
  }
  for (int inode = 0; inode != 4; ++inode) // fill n+1 position
  {
    for (int isd = 0; isd != 3; ++isd)
    {
      xyzt(isd,inode+4) = posnp(isd,inode);
    }
  }
  return xyzt;
}
//
// Print method
//
std::string XFEM::SpaceTimeBoundaryCell::toString() const
{
    std::stringstream s;
    s << "SpaceTimeBoundaryCell: " << getBeleId() << endl;
//    MCONST_FOREACH(vector< vector<double> >, coordinate, nodalpos_xi_domain_)
//    {
//        s << "[";
//        MPFOREACH(vector<double>, val, coordinate)
//        {
//            s << *val << " ";
//        };
//        s << "]" << endl;
//    };
    return s.str();
}




     


#endif  // #ifdef CCADISCRET
