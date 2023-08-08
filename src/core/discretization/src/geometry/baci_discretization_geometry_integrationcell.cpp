/*----------------------------------------------------------------------*/
/*! \file

\brief integration cell classes for domain and boundary integration

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

\level 3

*----------------------------------------------------------------------*/


#include "baci_discretization_geometry_integrationcell.H"

#include "baci_discretization_geometry_element_coordtrafo.H"
#include "baci_discretization_geometry_element_volume.H"
#include "baci_io_gmsh.H"
#include "baci_lib_globalproblem.H"


/*----------------------------------------------------------------------*
 * get center in physical coordinates
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 1> CORE::GEO::IntCell::ComputePhysicalCenterPosition(
    const ::DRT::Element::DiscretizationType& distype,
    const CORE::LINALG::SerialDenseMatrix& xyze) const
{
  // center in local coordinates
  const CORE::LINALG::Matrix<3, 1> localcenterpos(
      CORE::DRT::UTILS::getLocalCenterPosition<3>(distype));
  // center in physical coordinates
  static CORE::LINALG::Matrix<3, 1> pyhsicalcenterpos;
  CORE::GEO::elementToCurrentCoordinates(distype, xyze, localcenterpos, pyhsicalcenterpos);
  return pyhsicalcenterpos;
}


////////////// Integration cell ////////////////////////////////////////

/*----------------------------------------------------------------------*
 * constructor int cell                                                 *
 *----------------------------------------------------------------------*/
CORE::GEO::IntCell::IntCell(const ::DRT::Element::DiscretizationType& distype)
    : distype_(distype), indomainplus_(false)
{
}



/*----------------------------------------------------------------------*
 * copy constructor                                                     *
 *----------------------------------------------------------------------*/
CORE::GEO::IntCell::IntCell(const IntCell& old)
    : distype_(old.distype_), indomainplus_(old.indomainplus_)
{
}



/*----------------------------------------------------------------------*
 * assignment operator                                                  *
 *----------------------------------------------------------------------*/
CORE::GEO::IntCell& CORE::GEO::IntCell::operator=(const IntCell& intcell)
{
  distype_ = intcell.distype_;
  return *this;
}



/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
CORE::GEO::IntCell::~IntCell() {}


int CORE::GEO::IntCell::NumNode() const
{
  return CORE::DRT::UTILS::getNumberOfElementNodes(Shape());
}


/*----------------------------------------------------------------------*
 * to std::string                                                            *
 *----------------------------------------------------------------------*/
std::string CORE::GEO::IntCell::toString() const { return ""; }



////////////// Domain integration cell //////////////////////////////////

/*----------------------------------------------------------------------*
 * Constructor Domain integration cell                                  *
 *----------------------------------------------------------------------*/
CORE::GEO::DomainIntCell::DomainIntCell(const ::DRT::Element::DiscretizationType& distype,
    const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates)
    : IntCell(distype),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = false;
}



/*----------------------------------------------------------------------*
 * Constructor Domain integration cell                  rasthofer 07/09 *
 *----------------------------------------------------------------------*/
CORE::GEO::DomainIntCell::DomainIntCell(const ::DRT::Element::DiscretizationType& distype,
    const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool indomainplus)
    : IntCell(distype),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = indomainplus;
}



/*----------------------------------------------------------------------*
 * Constructor                                                          *
 * Domain integration cell == xfem element if not intersected           *
 *----------------------------------------------------------------------*/
CORE::GEO::DomainIntCell::DomainIntCell(const ::DRT::Element::DiscretizationType& distype,
    const CORE::LINALG::SerialDenseMatrix& xyze_ele)
    : IntCell(distype),
      nodalpos_xi_domain_(CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype)),
      nodalpos_xyz_domain_(xyze_ele),
      phys_center_(ComputePhysicalCenterPosition(distype, xyze_ele))
{
  indomainplus_ = false;
}



/*----------------------------------------------------------------------*
 * Copy Constructor Domain integration cell                             *
 *----------------------------------------------------------------------*/
CORE::GEO::DomainIntCell::DomainIntCell(const DomainIntCell& old)
    : IntCell(old),
      nodalpos_xi_domain_(old.nodalpos_xi_domain_),
      nodalpos_xyz_domain_(old.nodalpos_xyz_domain_),
      phys_center_(old.phys_center_)
{
  indomainplus_ = old.indomainplus_;
}


/*----------------------------------------------------------------------*
 * assignment operator                                                  *
 *----------------------------------------------------------------------*/
CORE::GEO::DomainIntCell& CORE::GEO::DomainIntCell::operator=(
    const CORE::GEO::DomainIntCell& domainintcell)
{
  this->CORE::GEO::IntCell::operator=(domainintcell);
  nodalpos_xi_domain_ = domainintcell.nodalpos_xi_domain_;
  nodalpos_xyz_domain_ = domainintcell.nodalpos_xyz_domain_;
  phys_center_ = domainintcell.phys_center_;
  indomainplus_ = domainintcell.indomainplus_;
  return *this;
}



/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
CORE::GEO::DomainIntCell::~DomainIntCell() {}

/*----------------------------------------------------------------------*
 * output                                                           *
 *----------------------------------------------------------------------*/
static std::string PosToString(double x, double y, double z)
{
  std::ostringstream s;
  s << "(" << std::setw(14) << std::scientific << x << "," << std::setw(14) << std::scientific << y
    << "," << std::setw(14) << std::scientific << z << ")";
  return s.str();
}


/*----------------------------------------------------------------------*
 *  to std::string                                                           *
 *----------------------------------------------------------------------*/
std::string CORE::GEO::DomainIntCell::toString() const
{
  std::ostringstream s;
  int numpoints = nodalpos_xi_domain_.numCols();
  s << "DomainIntCell:" << std::endl;
  s << " position in xi coordinates: " << std::endl;
  for (int inode = 0; inode < numpoints; ++inode)
    s << "   "
      << PosToString(nodalpos_xi_domain_(0, inode), nodalpos_xi_domain_(1, inode),
             nodalpos_xi_domain_(2, inode))
      << std::endl;
  s << " position in xyz coordinates: " << std::endl;
  for (int inode = 0; inode < numpoints; ++inode)
    s << "   "
      << PosToString(nodalpos_xyz_domain_(0, inode), nodalpos_xyz_domain_(1, inode),
             nodalpos_xyz_domain_(2, inode))
      << std::endl;

  s << std::endl
    << " Center : " << PosToString(phys_center_(0), phys_center_(1), phys_center_(2)) << std::endl;

  //  s << phys_center_ << std::endl;
  return s.str();
}



/*----------------------------------------------------------------------*
 * write Gmsh file with cell data in local coordinates                  *
 *----------------------------------------------------------------------*/
void CORE::GEO::DomainIntCell::xiToGmsh(const std::string& filename) const
{
  const CORE::LINALG::SerialDenseMatrix& cellpos = this->CellNodalPosXiDomain();

  std::ofstream f_system(filename.c_str());
  f_system << "View \" "
           << "Bad Cell \" {\n";
  f_system << IO::GMSH::cellWithScalarToString(this->Shape(), 0.0, cellpos);
  f_system << "};\n";
  f_system << "View[0].Axes = 3;\nView[0].AxesMikado = 1;\n";
  f_system.close();
}



/*----------------------------------------------------------------------*
 * write Gmsh file with cell data in global coordinates                 *
 *----------------------------------------------------------------------*/
void CORE::GEO::DomainIntCell::xToGmsh(const std::string& filename) const
{
  const CORE::LINALG::SerialDenseMatrix& cellpos = this->CellNodalPosXYZ();

  std::ofstream f_system(filename.c_str());
  f_system << "View \" "
           << "Bad Cell \" {\n";
  f_system << IO::GMSH::cellWithScalarToString(this->Shape(), 0.0, cellpos);
  f_system << "};\n";
  f_system << "View[0].Axes = 3;\nView[0].AxesMikado = 1;\n";
  f_system.close();
}



/*----------------------------------------------------------------------*
 * compute volume in XiDomain coordinates
 *----------------------------------------------------------------------*/
double CORE::GEO::DomainIntCell::VolumeInXiDomain(const ::DRT::Element& ele) const
{
  const double volume_cell = CORE::GEO::ElementVolume(this->Shape(), nodalpos_xi_domain_);
  if (volume_cell <= 0.0)
  {
    std::cout << this->toString() << std::endl;
    this->xiToGmsh("cell_with_negative_volume.pos");
    dserror("GLOBAL ELEMENT NO.%i\n NEGATIVE VOLUME OF INTEGRATION CELL: %20.16f", ele.Id(),
        volume_cell);
  }
  const double normed_cell_volume =
      volume_cell / CORE::DRT::UTILS::getSizeInLocalCoordinates(ele.Shape());
  return normed_cell_volume;
}

/*----------------------------------------------------------------------*
 * compute volume in physical domain
 *----------------------------------------------------------------------*/
double CORE::GEO::DomainIntCell::VolumeInPhysicalDomain() const
{
  double vol = CORE::GEO::ElementVolume(this->Shape(), nodalpos_xyz_domain_);

  if (vol < 0.0) dserror("cell with negative volume");

  return vol;
}

////////////// Boundary integration cell ////////////////////////////////

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::GEO::BoundaryIntCell* CORE::GEO::BoundaryIntCell::Create(
    const ::DRT::Element::DiscretizationType& distype, const int& surface_ele_gid,
    const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix* eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool& indomainplus)
{
  const int probdim = ::DRT::Problem::Instance()->NDim();
  if (probdim == 3)
  {
    return new BoundaryIntCell(distype, surface_ele_gid, xfemEleDomainCoordinates,
        *eleBoundaryCoordinates, physDomainCoordinates, indomainplus);
  }

  // do this only, if it is no 3-dimensional problem
  switch (distype)
  {
    case ::DRT::Element::point1:
      return CORE::GEO::CreateConcreteBoundaryIntCell<::DRT::Element::point1>(surface_ele_gid,
          xfemEleDomainCoordinates, eleBoundaryCoordinates, physDomainCoordinates, indomainplus,
          probdim);
    case ::DRT::Element::line2:
      return CORE::GEO::CreateConcreteBoundaryIntCell<::DRT::Element::line2>(surface_ele_gid,
          xfemEleDomainCoordinates, eleBoundaryCoordinates, physDomainCoordinates, indomainplus,
          probdim);
    default:
      dserror("Unsupported cell type = %s", ::DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }
  // you should never reach this point
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------*
 * constructor Boundary integration cells                               *
 *----------------------------------------------------------------------*/
CORE::GEO::BoundaryIntCell::BoundaryIntCell(const ::DRT::Element::DiscretizationType& distype,
    const int surface_ele_gid, const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix& eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates)
    : IntCell(distype),
      surface_ele_gid_(surface_ele_gid),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xi_boundary_(eleBoundaryCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = true;
}

/*----------------------------------------------------------------------*
 * constructor Boundary integration cells                               *
 *----------------------------------------------------------------------*/
CORE::GEO::BoundaryIntCell::BoundaryIntCell(const ::DRT::Element::DiscretizationType& distype,
    const int surface_ele_gid, const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix& eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool indomainplus)
    : IntCell(distype),
      surface_ele_gid_(surface_ele_gid),
      nodalpos_xi_domain_(xfemEleDomainCoordinates),
      nodalpos_xi_boundary_(eleBoundaryCoordinates),
      nodalpos_xyz_domain_(physDomainCoordinates),
      phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = indomainplus;
}


/*----------------------------------------------------------------------*
 * copy constructor Boundary integration cells                          *
 *----------------------------------------------------------------------*/
CORE::GEO::BoundaryIntCell::BoundaryIntCell(const BoundaryIntCell& old)
    : IntCell(old),
      surface_ele_gid_(old.surface_ele_gid_),
      nodalpos_xi_domain_(old.nodalpos_xi_domain_),
      nodalpos_xi_boundary_(old.nodalpos_xi_boundary_),
      nodalpos_xyz_domain_(old.nodalpos_xyz_domain_)
{
  indomainplus_ = old.indomainplus_;
}



/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
CORE::GEO::BoundaryIntCell::~BoundaryIntCell() {}

/*----------------------------------------------------------------------*
 * assignment operator                                                  *
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*
 * to std::string                                                            *
 *----------------------------------------------------------------------*/
std::string CORE::GEO::BoundaryIntCell::toString() const
{
  std::ostringstream s;
  s << "BoundaryIntCell" << std::endl;
  s << nodalpos_xi_domain_ << std::endl;
  return s.str();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::BoundaryIntCell::Print(std::ostream& stream) const
{
  stream << toString() << "\n" << std::flush;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::GEO::BoundaryIntCell::BoundaryIntCell(
    ::DRT::Element::DiscretizationType distype, const int& surface_ele_gid)
    : IntCell(distype), surface_ele_gid_(surface_ele_gid), phys_center_(true)
{
  /* intentionally left blank */
}

////////////// Concrete boundary integration cell //////////////////////////////

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <::DRT::Element::DiscretizationType cellType>
CORE::GEO::BoundaryIntCell* CORE::GEO::CreateConcreteBoundaryIntCell(const int& surface_ele_gid,
    const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix* eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool& indomainplus,
    const unsigned& probDim)
{
  switch (probDim)
  {
    case 2:
      return new CORE::GEO::ConcreteBoundaryIntCell<2, cellType>(surface_ele_gid,
          xfemEleDomainCoordinates, eleBoundaryCoordinates, physDomainCoordinates, indomainplus);
    default:
      dserror("Unsupported problem dimension! (probDim == %d", probDim);
      exit(EXIT_FAILURE);
  }
  // you should never reach this point
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probDim, ::DRT::Element::DiscretizationType cellType, unsigned dim,
    unsigned numNodePerEle>
CORE::GEO::ConcreteBoundaryIntCell<probDim, cellType, dim, numNodePerEle>::ConcreteBoundaryIntCell(
    const int& surface_ele_gid, const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix* eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool& indomainplus)
    : CORE::GEO::BoundaryIntCell(cellType, surface_ele_gid),
      xyz_center_(this->phys_center_.A(), true)
{
  nodalpos_xi_domain_ = xfemEleDomainCoordinates;
  nodalpos_xi_boundary_ = (eleBoundaryCoordinates != nullptr ? *eleBoundaryCoordinates
                                                             : CORE::LINALG::SerialDenseMatrix());
  nodalpos_xyz_domain_ = physDomainCoordinates;
  indomainplus_ = indomainplus;
  CORE::GEO::ComputePhysicalCenterPosition<probDim, cellType>(physDomainCoordinates, xyz_center_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <unsigned probDim, ::DRT::Element::DiscretizationType cellType, unsigned dim,
    unsigned numNodePerEle>
CORE::GEO::BoundaryIntCell&
CORE::GEO::ConcreteBoundaryIntCell<probDim, cellType, dim, numNodePerEle>::operator=(
    const CORE::GEO::BoundaryIntCell& boundaryintcell)
{
  return operator=(
      dynamic_cast<const ConcreteBoundaryIntCell<probDim, cellType, dim, numNodePerEle>&>(
          boundaryintcell));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <unsigned probDim, ::DRT::Element::DiscretizationType cellType, unsigned dim,
    unsigned numNodePerEle>
CORE::GEO::BoundaryIntCell&
CORE::GEO::ConcreteBoundaryIntCell<probDim, cellType, dim, numNodePerEle>::operator=(
    const CORE::GEO::ConcreteBoundaryIntCell<probDim, cellType, dim, numNodePerEle>&
        boundaryintcell)
{
  this->CORE::GEO::BoundaryIntCell::operator=(boundaryintcell);
  // here we set a view to the base class member
  this->xyz_center_ = CORE::LINALG::Matrix<probDim, 1>(this->phys_center_.A(), true);

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <> /* function specialization */
void CORE::GEO::ComputePhysicalCenterPosition<2, ::DRT::Element::line2>(
    const CORE::LINALG::SerialDenseMatrix& xyze, CORE::LINALG::Matrix<2, 1>& phys_center)
{
  const ::DRT::Element::DiscretizationType cellType = ::DRT::Element::line2;
  const unsigned dim = CORE::DRT::UTILS::DisTypeToDim<cellType>::dim;
  CORE::LINALG::Matrix<dim, 1> localcenterpos =
      CORE::DRT::UTILS::getLocalCenterPosition<1>(cellType);
  CORE::GEO::elementToCurrentCoordinates(cellType, xyze, localcenterpos, phys_center);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <> /* function specialization */
void CORE::GEO::ComputePhysicalCenterPosition<2, ::DRT::Element::point1>(
    const CORE::LINALG::SerialDenseMatrix& xyze, CORE::LINALG::Matrix<2, 1>& phys_center)
{
  phys_center.SetCopy(xyze.values());
}

template class CORE::GEO::ConcreteBoundaryIntCell<2, ::DRT::Element::point1>;
template class CORE::GEO::ConcreteBoundaryIntCell<2, ::DRT::Element::line2>;

template void CORE::GEO::ComputePhysicalCenterPosition<2, ::DRT::Element::point1>(
    const CORE::LINALG::SerialDenseMatrix& xyze, CORE::LINALG::Matrix<2, 1>& phys_center);
template void CORE::GEO::ComputePhysicalCenterPosition<2, ::DRT::Element::line2>(
    const CORE::LINALG::SerialDenseMatrix& xyze, CORE::LINALG::Matrix<2, 1>& phys_center);

template CORE::GEO::BoundaryIntCell*
CORE::GEO::CreateConcreteBoundaryIntCell<::DRT::Element::point1>(const int& surface_ele_gid,
    const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix* eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool& indomainplus,
    const unsigned& probDim);
template CORE::GEO::BoundaryIntCell*
CORE::GEO::CreateConcreteBoundaryIntCell<::DRT::Element::line2>(const int& surface_ele_gid,
    const CORE::LINALG::SerialDenseMatrix& xfemEleDomainCoordinates,
    const CORE::LINALG::SerialDenseMatrix* eleBoundaryCoordinates,
    const CORE::LINALG::SerialDenseMatrix& physDomainCoordinates, const bool& indomainplus,
    const unsigned& probDim);
