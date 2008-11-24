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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell(
    const int           bele_id,
    const LINALG::SerialDenseMatrix&     posnp,
    const LINALG::SerialDenseMatrix&     posn
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
        
LINALG::SerialDenseMatrix XFEM::SpaceTimeBoundaryCell::getLinearPositionArray(
    const LINALG::SerialDenseMatrix&      posnp,
    const LINALG::SerialDenseMatrix&      posn
    ) const
{
  LINALG::SerialDenseMatrix xyzt(3,8);
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
