/*----------------------------------------------------------------------*/
/*! \file

\brief integration cell classes for domain and boundary integration

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE 4C SOON

\level 3

*----------------------------------------------------------------------*/


#include "4C_fem_geometry_integrationcell.hpp"

#include "4C_fem_geometry_element_coordtrafo.hpp"
#include "4C_fem_geometry_element_volume.hpp"
#include "4C_io_gmsh.hpp"

FOUR_C_NAMESPACE_OPEN


Core::LinAlg::Matrix<3, 1> Core::Geo::IntCell::compute_physical_center_position(
    const Core::FE::CellType& distype, const Core::LinAlg::SerialDenseMatrix& xyze) const
{
  // center in local coordinates
  const Core::LinAlg::Matrix<3, 1> localcenterpos(Core::FE::getLocalCenterPosition<3>(distype));
  // center in physical coordinates
  static Core::LinAlg::Matrix<3, 1> pyhsicalcenterpos;
  Core::Geo::elementToCurrentCoordinates(distype, xyze, localcenterpos, pyhsicalcenterpos);
  return pyhsicalcenterpos;
}


////////////// Integration cell ////////////////////////////////////////

Core::Geo::IntCell::IntCell(const Core::FE::CellType& distype)
    : distype_(distype), indomainplus_(false)
{
}

Core::Geo::IntCell::IntCell(const IntCell& old)
    : distype_(old.distype_), indomainplus_(old.indomainplus_)
{
}

Core::Geo::IntCell& Core::Geo::IntCell::operator=(const IntCell& intcell)
{
  distype_ = intcell.distype_;
  return *this;
}


////////////// Boundary integration cell ////////////////////////////////

Core::Geo::BoundaryIntCell::BoundaryIntCell(const Core::FE::CellType& distype,
    const int surface_ele_gid, const Core::LinAlg::SerialDenseMatrix& xfemEleDomainCoordinates,
    const Core::LinAlg::SerialDenseMatrix& eleBoundaryCoordinates,
    const Core::LinAlg::SerialDenseMatrix& physDomainCoordinates)
    : IntCell(distype),
      surface_ele_gid_(surface_ele_gid),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xi_boundary_(eleBoundaryCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(compute_physical_center_position(distype, physDomainCoordinates))
{
  indomainplus_ = true;
}

Core::Geo::BoundaryIntCell::BoundaryIntCell(const Core::FE::CellType& distype,
    const int surface_ele_gid, const Core::LinAlg::SerialDenseMatrix& xfemEleDomainCoordinates,
    const Core::LinAlg::SerialDenseMatrix& eleBoundaryCoordinates,
    const Core::LinAlg::SerialDenseMatrix& physDomainCoordinates, const bool indomainplus)
    : IntCell(distype),
      surface_ele_gid_(surface_ele_gid),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xi_boundary_(eleBoundaryCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(compute_physical_center_position(distype, physDomainCoordinates))
{
  indomainplus_ = indomainplus;
}

Core::Geo::BoundaryIntCell::BoundaryIntCell(const BoundaryIntCell& old)
    : IntCell(old),
      surface_ele_gid_(old.surface_ele_gid_),
      nodalpos_xi_domain_(old.nodalpos_xi_domain_),
      nodalpos_xi_boundary_(old.nodalpos_xi_boundary_),
      nodalpos_xyz_domain_(old.nodalpos_xyz_domain_)
{
  indomainplus_ = old.indomainplus_;
}

Core::Geo::BoundaryIntCell& Core::Geo::BoundaryIntCell::operator=(
    const Core::Geo::BoundaryIntCell& boundaryintcell)
{
  this->Core::Geo::IntCell::operator=(boundaryintcell);
  surface_ele_gid_ = boundaryintcell.surface_ele_gid_;
  nodalpos_xi_domain_ = boundaryintcell.nodalpos_xi_domain_;
  nodalpos_xi_boundary_ = boundaryintcell.nodalpos_xi_boundary_;
  nodalpos_xyz_domain_ = boundaryintcell.nodalpos_xyz_domain_;
  indomainplus_ = boundaryintcell.indomainplus_;
  return *this;
}

Core::Geo::BoundaryIntCell::BoundaryIntCell(Core::FE::CellType distype, const int& surface_ele_gid)
    : IntCell(distype), surface_ele_gid_(surface_ele_gid), phys_center_(true)
{
  /* intentionally left blank */
}

FOUR_C_NAMESPACE_CLOSE
