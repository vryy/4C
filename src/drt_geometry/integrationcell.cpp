/*!----------------------------------------------------------------------
\file integrationcell.cpp

\brief integration cell classes for domain and boundary integration

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

<pre>
\level 3
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*----------------------------------------------------------------------*/


#include "integrationcell.H"
#include "element_coordtrafo.H"
#include "element_volume.H"
#include "../drt_io/io_gmsh.H"


/*----------------------------------------------------------------------*
 * get center in physical coordinates
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> GEO::IntCell::ComputePhysicalCenterPosition(
  const DRT::Element::DiscretizationType&   distype,
  const LINALG::SerialDenseMatrix&          xyze) const
{
  // center in local coordinates
  const LINALG::Matrix<3,1> localcenterpos(DRT::UTILS::getLocalCenterPosition<3>(distype));
  // center in physical coordinates
  static LINALG::Matrix<3,1> pyhsicalcenterpos;
  GEO::elementToCurrentCoordinates(distype, xyze, localcenterpos, pyhsicalcenterpos);
  return pyhsicalcenterpos;
}


////////////// Integration cell ////////////////////////////////////////

/*----------------------------------------------------------------------*
 * constructor int cell                                                 *
 *----------------------------------------------------------------------*/
GEO::IntCell::IntCell(
        const DRT::Element::DiscretizationType& distype) :
            distype_(distype)
{}



/*----------------------------------------------------------------------*
 * copy constructor                                                     *
 *----------------------------------------------------------------------*/
GEO::IntCell::IntCell(
        const IntCell& old) :
            distype_(old.distype_)
{}



/*----------------------------------------------------------------------*
 * assignment operator                                                  *
 *----------------------------------------------------------------------*/
GEO::IntCell& GEO::IntCell::operator=(const IntCell& intcell)
{
  distype_ = intcell.distype_;
  return *this;
}



/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
GEO::IntCell::~IntCell()
{}


int GEO::IntCell::NumNode() const
{
  return DRT::UTILS::getNumberOfElementNodes(Shape());
}


/*----------------------------------------------------------------------*
 * to std::string                                                            *
 *----------------------------------------------------------------------*/
std::string GEO::IntCell::toString() const
{
  return "";
}



////////////// Domain integration cell //////////////////////////////////

/*----------------------------------------------------------------------*
 * Constructor Domain integration cell                                  *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType&     distype,
        const LINALG::SerialDenseMatrix&            xfemEleDomainCoordinates,
        const LINALG::SerialDenseMatrix&            physDomainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_(xfemEleDomainCoordinates),
            nodalpos_xyz_domain_(physDomainCoordinates),
            phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = false;
}



/*----------------------------------------------------------------------*
 * Constructor Domain integration cell                  rasthofer 07/09 *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType&     distype,
        const LINALG::SerialDenseMatrix&            xfemEleDomainCoordinates,
        const LINALG::SerialDenseMatrix&            physDomainCoordinates,
        const bool                                  indomainplus) :
            IntCell(distype),
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
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType&   distype,
        const LINALG::SerialDenseMatrix&          xyze_ele) :
            IntCell(distype),
            nodalpos_xi_domain_(DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype)),
            nodalpos_xyz_domain_(xyze_ele),
            phys_center_(ComputePhysicalCenterPosition(distype, xyze_ele))
{
  indomainplus_ = false;
}



/*----------------------------------------------------------------------*
 * Copy Constructor Domain integration cell                             *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
      const DomainIntCell& old) :
          IntCell(old),
          nodalpos_xi_domain_(old.nodalpos_xi_domain_),
          nodalpos_xyz_domain_(old.nodalpos_xyz_domain_),
          phys_center_(old.phys_center_)
{
  indomainplus_ = old.indomainplus_;
}


/*----------------------------------------------------------------------*
 * assignment operator                                                  *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell& GEO::DomainIntCell::operator=(
    const GEO::DomainIntCell& domainintcell)
{
  this->GEO::IntCell::operator = (domainintcell);
  nodalpos_xi_domain_ = domainintcell.nodalpos_xi_domain_;
  nodalpos_xyz_domain_ = domainintcell.nodalpos_xyz_domain_;
  phys_center_ = domainintcell.phys_center_;
  indomainplus_ = domainintcell.indomainplus_;
  return *this;
}



/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::~DomainIntCell()
{}



/*----------------------------------------------------------------------*
 * output                                                           *
 *----------------------------------------------------------------------*/
static std::string PosToString(double x, double y, double z)
{
  std::ostringstream s;
  s << "(" << std::setw(14) << std::scientific << x <<
       "," << std::setw(14) << std::scientific << y <<
       "," << std::setw(14) << std::scientific << z << ")";
  return s.str();
}


/*----------------------------------------------------------------------*
 *  to std::string                                                           *
 *----------------------------------------------------------------------*/
std::string GEO::DomainIntCell::toString() const
{
  std::ostringstream s;
  int numpoints = nodalpos_xi_domain_.N();
  s << "DomainIntCell:" << std::endl;
  s << " position in xi coordinates: " << std::endl;
  for (int inode = 0; inode < numpoints; ++inode)
    s << "   " << PosToString(nodalpos_xi_domain_(0,inode),nodalpos_xi_domain_(1,inode),nodalpos_xi_domain_(2,inode)) << std::endl;
  s << " position in xyz coordinates: " << std::endl;
  for (int inode = 0; inode < numpoints; ++inode)
    s << "   " << PosToString(nodalpos_xyz_domain_(0,inode),nodalpos_xyz_domain_(1,inode),nodalpos_xyz_domain_(2,inode)) << std::endl;

  s << std::endl << " Center : " << PosToString(phys_center_(0),phys_center_(1),phys_center_(2)) << std::endl;

//  s << phys_center_ << std::endl;
  return s.str();
}



/*----------------------------------------------------------------------*
 * write Gmsh file with cell data in local coordinates                  *
 *----------------------------------------------------------------------*/
void GEO::DomainIntCell::xiToGmsh(const std::string& filename) const
{
  const LINALG::SerialDenseMatrix& cellpos = this->CellNodalPosXiDomain();

  std::ofstream f_system(filename.c_str());
  f_system << "View \" " << "Bad Cell \" {\n";
  f_system << IO::GMSH::cellWithScalarToString(this->Shape(), 0.0, cellpos);
  f_system << "};\n";
  f_system << "View[0].Axes = 3;\nView[0].AxesMikado = 1;\n";
  f_system.close();
}



/*----------------------------------------------------------------------*
 * write Gmsh file with cell data in global coordinates                 *
 *----------------------------------------------------------------------*/
void GEO::DomainIntCell::xToGmsh(const std::string& filename) const
{
  const LINALG::SerialDenseMatrix& cellpos = this->CellNodalPosXYZ();

  std::ofstream f_system(filename.c_str());
  f_system << "View \" " << "Bad Cell \" {\n";
  f_system << IO::GMSH::cellWithScalarToString(this->Shape(), 0.0, cellpos);
  f_system << "};\n";
  f_system << "View[0].Axes = 3;\nView[0].AxesMikado = 1;\n";
  f_system.close();
}



/*----------------------------------------------------------------------*
 * compute volume in XiDomain coordinates
 *----------------------------------------------------------------------*/
double GEO::DomainIntCell::VolumeInXiDomain(
        const DRT::Element&           ele
        ) const
{
  const double volume_cell = GEO::ElementVolume(this->Shape(), nodalpos_xi_domain_);
  if(volume_cell <= 0.0)
  {
    std::cout << this->toString() << std::endl;
    this->xiToGmsh("cell_with_negative_volume.pos");
    dserror("GLOBAL ELEMENT NO.%i\n NEGATIVE VOLUME OF INTEGRATION CELL: %20.16f", ele.Id(), volume_cell);
  }
  const double normed_cell_volume = volume_cell/DRT::UTILS::getSizeInLocalCoordinates(ele.Shape());
  return normed_cell_volume;
}

/*----------------------------------------------------------------------*
 * compute volume in physical domain
 *----------------------------------------------------------------------*/
double GEO::DomainIntCell::VolumeInPhysicalDomain() const
{
  double vol = GEO::ElementVolume(this->Shape(), nodalpos_xyz_domain_);

  if(vol < 0.0)
    dserror("cell with negative volume");

  return vol;
}

////////////// Boundary integration cell ////////////////////////////////

/*----------------------------------------------------------------------*
 * constructor Boundary integration cells                               *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType&   distype,
        const int                                 surface_ele_gid,
        const LINALG::SerialDenseMatrix&          xfemEleDomainCoordinates,
        const LINALG::SerialDenseMatrix&          eleBoundaryCoordinates,
        const LINALG::SerialDenseMatrix&          physDomainCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_(    xfemEleDomainCoordinates),
            nodalpos_xi_boundary_(  eleBoundaryCoordinates),
            nodalpos_xyz_domain_(   physDomainCoordinates),
            phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = true;
}

/*----------------------------------------------------------------------*
 * constructor Boundary integration cells                               *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType&   distype,
        const int                                 surface_ele_gid,
        const LINALG::SerialDenseMatrix&          xfemEleDomainCoordinates,
        const LINALG::SerialDenseMatrix&          eleBoundaryCoordinates,
        const LINALG::SerialDenseMatrix&          physDomainCoordinates,
        const bool                                indomainplus) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_(    xfemEleDomainCoordinates),
            nodalpos_xi_boundary_(  eleBoundaryCoordinates),
            nodalpos_xyz_domain_(   physDomainCoordinates),
            phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates))
{
  indomainplus_ = indomainplus;
}


/*----------------------------------------------------------------------*
 * copy constructor Boundary integration cells                          *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            surface_ele_gid_(old.surface_ele_gid_),
            nodalpos_xi_domain_(    old.nodalpos_xi_domain_),
            nodalpos_xi_boundary_(  old.nodalpos_xi_boundary_),
            nodalpos_xyz_domain_(   old.nodalpos_xyz_domain_)
{
  indomainplus_ = old.indomainplus_;
}



/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::~BoundaryIntCell()
{}




/*----------------------------------------------------------------------*
 * assignment operator                                                  *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell& GEO::BoundaryIntCell::operator=(const GEO::BoundaryIntCell& boundaryintcell)
{
  this->GEO::IntCell::operator = (boundaryintcell);
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
std::string GEO::BoundaryIntCell::toString() const
{
  std::ostringstream s;
  s << "BoundaryIntCell" << std::endl;
  s << nodalpos_xi_domain_ << std::endl;
  return s.str();
}

