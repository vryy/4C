/*!----------------------------------------------------------------------
\file integrationcell.cpp

\brief integration cell classes for domain and boundary integration

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "integrationcell.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_io/io_gmsh.H"
//#include "../drt_io/io_gmsh_xfem_extension.H"



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
 * destructor                                                           *
 *----------------------------------------------------------------------*/
GEO::IntCell::~IntCell()
{}

 
/*----------------------------------------------------------------------*
 * to string                                                            *
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
            phys_center_(ComputePhysicalCenterPosition(distype, physDomainCoordinates)),
            label_(-1)
{}

         
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
            phys_center_(ComputePhysicalCenterPosition(distype, xyze_ele)),
            label_(-1)
{}
       
        
/*----------------------------------------------------------------------*
 * Constructor Domain integration cell                                  *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
      const DomainIntCell& old) :
          IntCell(old),
          nodalpos_xi_domain_(old.nodalpos_xi_domain_),
          nodalpos_xyz_domain_(old.nodalpos_xyz_domain_),
          phys_center_(old.phys_center_),
          label_(old.label_)
{}


/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::~DomainIntCell()
{}


static string PosToString(double x, double y, double z)
{
  std::stringstream s;
  s << "(" << std::setw(14) << scientific << x <<
       "," << std::setw(14) << scientific << y <<
       "," << std::setw(14) << scientific << z << ")";
  return s.str();
}

/*----------------------------------------------------------------------*
 *  to string                                                           *
 *----------------------------------------------------------------------*/
std::string GEO::DomainIntCell::toString() const
{
  std::stringstream s;
  s << "DomainIntCell:" << endl;
  s << " position in xi coordinates: " << endl;
  for (int inode = 0; inode < 4; ++inode)
    s << "   " << PosToString(nodalpos_xi_domain_(0,inode),nodalpos_xi_domain_(1,inode),nodalpos_xi_domain_(2,inode)) << endl;
  s << " position in xyz coordinates: " << endl;
  for (int inode = 0; inode < 4; ++inode)
    s << "   " << PosToString(nodalpos_xyz_domain_(0,inode),nodalpos_xyz_domain_(1,inode),nodalpos_xyz_domain_(2,inode)) << endl;
  
  s << endl << " Center : " << PosToString(phys_center_(0),phys_center_(1),phys_center_(2)) << endl;
  
//  s << phys_center_ << endl;
  return s.str();
}


/*----------------------------------------------------------------------*
 *  to string                                                           *
 *----------------------------------------------------------------------*/
void GEO::DomainIntCell::toGmsh(const std::string& filename) const
{
  const LINALG::SerialDenseMatrix& cellpos = this->CellNodalPosXiDomain();
  
  std::ofstream f_system(filename.c_str());
  f_system << "View \" " << "Bad Cell \" {\n";
  f_system << IO::GMSH::cellWithScalarToString(this->Shape(), 0.0, cellpos) << "\n";
  f_system << "};\n";
  f_system << "View[0].Axes = 3;\nView[0].AxesMikado = 1;\n";
  f_system.close();
}


/*----------------------------------------------------------------------*
 * get center in physical coordinates
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> GEO::DomainIntCell::ComputePhysicalCenterPosition(
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


/*----------------------------------------------------------------------*
 * set xfem label, if fluid label = 0; if solid label = solid id
 *----------------------------------------------------------------------*/
void GEO::DomainIntCell::setLabel(const int   label)
{
  label_ = label;
}


/*----------------------------------------------------------------------*
 * check for coplanar corner points
 *----------------------------------------------------------------------*/
bool GEO::DomainIntCell::CoplanarCornerPoints() const
{
  const bool tetcell = ((Shape() == DRT::Element::tet4) or (Shape() == DRT::Element::tet10));
  dsassert(tetcell,"test for coplanar points only for tetrahedral integration cells");
  
  // create plane of first 3 nodes
  LINALG::Matrix<3,3> plane;
  for (int ipoint = 0; ipoint < 3; ++ipoint)
    for (int isd = 0; isd < 3; ++isd)
      plane(isd,ipoint) = nodalpos_xi_domain_(isd,ipoint);
  
  // 4th node is the test point
  LINALG::Matrix<3,1> testpoint;
  for (int isd = 0; isd < 3; ++isd)
    testpoint(isd,0) = nodalpos_xi_domain_(isd,3);
  
  const bool coplanar_tet  = GEO::testForCoplanarTet(nodalpos_xi_domain_);
  const bool coplanar_tet2 = GEO::pointsInPlaneSurfaceElement(plane,testpoint);
  
  if (coplanar_tet != coplanar_tet2)
    dserror("which method is better?");
  
  // check for coplanar points
  return GEO::pointsInPlaneSurfaceElement(plane,testpoint);
}

/*----------------------------------------------------------------------*
 * compute volume in XiDomain coordinates
 *----------------------------------------------------------------------*/
double GEO::DomainIntCell::VolumeInXiDomain(
        const DRT::Element&           ele
        ) const
{
  if (CoplanarCornerPoints())
  {
    cout << "coplanar" << endl;
    return 0.0;
  }

  DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
  switch (this->Shape())
  {
    case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
    {
      gaussrule = DRT::UTILS::intrule_hex_8point;
      break;
    }
    case DRT::Element::tet4: case DRT::Element::tet10:
    {
      gaussrule = DRT::UTILS::intrule_tet_4point;
      break;
    }
    default:
      dserror("add your element type here...");
  }
  
  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  double volume_cell = 0.0;
  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point in cell coordinates \eta
    static LINALG::Matrix<3,1> pos_eta_domain;
    pos_eta_domain(0) = intpoints.qxg[iquad][0];
    pos_eta_domain(1) = intpoints.qxg[iquad][1];
    pos_eta_domain(2) = intpoints.qxg[iquad][2];
    
    const double detcell = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*this, pos_eta_domain);
    const double fac = intpoints.qwgt[iquad]*detcell;

    if(detcell < 0.0)
    {
      cout << scientific << detcell << endl;
      cout << this->toString() << endl;
      this->toGmsh("cell.pos");
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT OF INTEGRATION CELL: %20.16f", ele.Id(), detcell);
    }    
    volume_cell += fac;

  } // end loop over gauss points
  return volume_cell/DRT::UTILS::getSizeInLocalCoordinates(ele.Shape());
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
            nodalpos_xyz_domain_(   physDomainCoordinates)
{}
         
        
/*----------------------------------------------------------------------*
 * constructor Boundary integration cells                               *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            surface_ele_gid_(old.surface_ele_gid_),
            nodalpos_xi_domain_(    old.nodalpos_xi_domain_),
            nodalpos_xi_boundary_(  old.nodalpos_xi_boundary_),
            nodalpos_xyz_domain_(   old.nodalpos_xyz_domain_)
{}
 
 
/*----------------------------------------------------------------------*
 * destructor                                                           *
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::~BoundaryIntCell()
{}
 
 
/*----------------------------------------------------------------------*
 * to string                                                            *
 *----------------------------------------------------------------------*/    
std::string GEO::BoundaryIntCell::toString() const
{
  std::stringstream s;
  s << "BoundaryIntCell" << endl;
  s << nodalpos_xi_domain_ << endl;
  return s.str();
}


#endif  // #ifdef CCADISCRET
