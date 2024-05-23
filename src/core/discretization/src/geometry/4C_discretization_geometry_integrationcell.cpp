/*----------------------------------------------------------------------*/
/*! \file

\brief integration cell classes for domain and boundary integration

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE 4C SOON

\level 3

*----------------------------------------------------------------------*/


#include "4C_discretization_geometry_integrationcell.hpp"

#include "4C_discretization_geometry_element_coordtrafo.hpp"
#include "4C_discretization_geometry_element_volume.hpp"
#include "4C_io_gmsh.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::LINALG::Matrix<3, 1> CORE::GEO::IntCell::compute_physical_center_position(
    const CORE::FE::CellType& distype, const CORE::LINALG::SerialDenseMatrix& xyze) const
{
  // center in local coordinates
  const CORE::LINALG::Matrix<3, 1> localcenterpos(CORE::FE::getLocalCenterPosition<3>(distype));
  // center in physical coordinates
  static CORE::LINALG::Matrix<3, 1> pyhsicalcenterpos;
  CORE::GEO::elementToCurrentCoordinates(distype, xyze, localcenterpos, pyhsicalcenterpos);
  return pyhsicalcenterpos;
}


////////////// Integration cell ////////////////////////////////////////

CORE::GEO::IntCell::IntCell(const CORE::FE::CellType& distype)
    : distype_(distype), indomainplus_(false)
{
}

CORE::GEO::IntCell::IntCell(const IntCell& old)
    : distype_(old.distype_), indomainplus_(old.indomainplus_)
{
}

CORE::GEO::IntCell& CORE::GEO::IntCell::operator=(const IntCell& intcell)
{
  distype_ = intcell.distype_;
  return *this;
}


////////////// Boundary integration cell ////////////////////////////////

CORE::GEO::BoundaryIntCell::BoundaryIntCell(const CORE::FE::CellType& distype,
    const int surface_ele_gid, const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix& eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates)
    : IntCell(distype),
      surface_ele_gid_(surface_ele_gid),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xi_boundary_(eleBoundaryCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(compute_physical_center_position(distype, physDomainCoordinates))
{
  indomainplus_ = true;
}

CORE::GEO::BoundaryIntCell::BoundaryIntCell(const CORE::FE::CellType& distype,
    const int surface_ele_gid, const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix& eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool indomainplus)
    : IntCell(distype),
      surface_ele_gid_(surface_ele_gid),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xi_boundary_(eleBoundaryCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(compute_physical_center_position(distype, physDomainCoordinates))
{
  indomainplus_ = indomainplus;
}

CORE::GEO::BoundaryIntCell::BoundaryIntCell(const BoundaryIntCell& old)
    : IntCell(old),
      surface_ele_gid_(old.surface_ele_gid_),
      nodalpos_xi_domain_(old.nodalpos_xi_domain_),
      nodalpos_xi_boundary_(old.nodalpos_xi_boundary_),
      nodalpos_xyz_domain_(old.nodalpos_xyz_domain_)
{
  indomainplus_ = old.indomainplus_;
}

CORE::GEO::BoundaryIntCell& CORE::GEO::BoundaryIntCell::operator=(
    const CORE::GEO::BoundaryIntCell& boundaryintcell)
{
  this->CORE::GEO::IntCell::operator=(boundaryintcell);
  surface_ele_gid_ = boundaryintcell.surface_ele_gid_;
  nodalpos_xi_domain_ = boundaryintcell.nodalpos_xi_domain_;
  nodalpos_xi_boundary_ = boundaryintcell.nodalpos_xi_boundary_;
  nodalpos_xyz_domain_ = boundaryintcell.nodalpos_xyz_domain_;
  indomainplus_ = boundaryintcell.indomainplus_;
  return *this;
}

CORE::GEO::BoundaryIntCell::BoundaryIntCell(CORE::FE::CellType distype, const int& surface_ele_gid)
    : IntCell(distype), surface_ele_gid_(surface_ele_gid), phys_center_(true)
{
  /* intentionally left blank */
}

FOUR_C_NAMESPACE_CLOSE
