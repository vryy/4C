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
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_fem_general/drt_utils_integration.H"



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


/*----------------------------------------------------------------------*
 *  to string                                                           *
 *----------------------------------------------------------------------*/
std::string GEO::DomainIntCell::toString() const
{
  std::stringstream s;
  s << "DomainIntCell" << endl;
  s << nodalpos_xi_domain_ << endl;
  return s.str();
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
 * compute size in XiDomain coordinates
 *----------------------------------------------------------------------*/
double GEO::DomainIntCell::SizeXiDomain(
        const DRT::Element&           ele
        ) const
{
  double volume_cell = 0.0;
  
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

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point in cell coordinates \eta
    static LINALG::Matrix<3,1> pos_eta_domain;
    pos_eta_domain(0) = intpoints.qxg[iquad][0];
    pos_eta_domain(1) = intpoints.qxg[iquad][1];
    pos_eta_domain(2) = intpoints.qxg[iquad][2];


    // coordinates of the current integration point in element coordinates \xi
    static LINALG::Matrix<3,1> posXiDomain;
    GEO::mapEtaToXi3D<XFEM::xfem_assembly>(*this, pos_eta_domain, posXiDomain);
    const double detcell = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*this, pos_eta_domain);

    const double fac = intpoints.qwgt[iquad]*detcell;

    if(detcell < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele.Id(), detcell);
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
