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
    const int                               bele_id,
    const DRT::Element::DiscretizationType  surf_distype,
    const LINALG::SerialDenseMatrix&        posnp,
    const LINALG::SerialDenseMatrix&        posn
    ) :
      bele_id_(bele_id),
      surf_distype_(surf_distype),
      num_timestep_(2),
      xyzt_(getLinearPositionArray(posnp,posn))
{
  if (surf_distype != DRT::Element::quad4)
    dserror("SpaceTimeBoundaryCell implemented only for surfaces of hex8 solid elements.");
  return;
}
  
    
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell() :
      bele_id_(-1),
      num_timestep_(-1),
      xyzt_()
{
  return;
}


/*----------------------------------------------------------------------*
 * copy constructor
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell(
    const SpaceTimeBoundaryCell& old
    ) :
      bele_id_(old.bele_id_),
      surf_distype_(old.surf_distype_),
      num_timestep_(old.num_timestep_),
      xyzt_(old.xyzt_)
{
    return;
}
   
    
/*----------------------------------------------------------------------*
 |  assignment operatur                                      u.may 04/09|
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell& XFEM::SpaceTimeBoundaryCell::operator = (const XFEM::SpaceTimeBoundaryCell& old)
{  
  bele_id_ = old.bele_id_;
  surf_distype_ = old.surf_distype_;
  num_timestep_ = old.num_timestep_;
  xyzt_ = old.xyzt_;
  return *this;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix XFEM::SpaceTimeBoundaryCell::getLinearPositionArray(
    const LINALG::SerialDenseMatrix&      posnp,
    const LINALG::SerialDenseMatrix&      posn
    ) const
{
  const int numnode_boundary = DRT::UTILS::getNumberOfElementNodes(surf_distype_);
  const int nsd = 3;

  LINALG::SerialDenseMatrix xyzt(nsd,numnode_boundary*num_timestep_);
  for (int inode = 0; inode != numnode_boundary; ++inode) // fill n   position
  {
    for (int isd = 0; isd != nsd; ++isd)
    {
      xyzt(isd,inode  ) = posn(isd,inode);
    }
  }
  for (int inode = 0; inode != numnode_boundary; ++inode) // fill n+1 position
  {
    for (int isd = 0; isd != nsd; ++isd)
    {
      xyzt(isd,inode+numnode_boundary) = posnp(isd,inode);
    }
  }
  return xyzt;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string XFEM::SpaceTimeBoundaryCell::toString() const
{
    std::stringstream s;
    s << "SpaceTimeBoundaryCell: " << getBeleId() << endl;
    return s.str();
}


#endif  // #ifdef CCADISCRET
