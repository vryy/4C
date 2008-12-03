/*!
\file interfacexfsi.cpp

\brief provides a class that represents an enriched physical scalar field for XFSI problems

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include "interfacexfsi.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_serialdensevector.H"

#include "xfem_condition.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_gmsh_xfem_extension.H"
#include "../drt_geometry/intersection.H"
#include "../drt_geometry/intersection_exp.H"
#include "coordinate_transformation.H"
#include "../drt_fem_general/drt_utils_integration.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandleXFSI::InterfaceHandleXFSI(
    const Teuchos::RCP<DRT::Discretization>&  xfemdis,
    const Teuchos::RCP<DRT::Discretization>&  cutterdis,
    const int step,
    const bool experimental_intersection
    ) : InterfaceHandle(xfemdis),
        cutterdis_(cutterdis)
{
  if (xfemdis->Comm().MyPID() == 0)
    std::cout << "Constructing InterfaceHandle" << std::endl;
      
  FillCurrentCutterPositionMap(cutterdis, *cutterdis->GetState("idispcolnp"), cutterposnp_);
  FillCurrentCutterPositionMap(cutterdis, *cutterdis->GetState("idispcoln") , cutterposn_ );
  
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = false;
  
  const int myrank = xfemdis_->Comm().MyPID();

  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_uncut_elements_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_uncut_elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis_);
    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, cutterposnp_);
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }

  CheckXFEMElementSize();
  
  elementalDomainIntCells_.clear();
  elementalBoundaryIntCells_.clear();

  
  if (experimental_intersection)
  {
    GEO::IntersectionExp is;
    is.computeIntersection(xfemdis, cutterdis, cutterposnp_, elementalDomainIntCells_, elementalBoundaryIntCells_);
  }
  else
  {
    GEO::Intersection is;
    is.computeIntersection(xfemdis, cutterdis, cutterposnp_, elementalDomainIntCells_, elementalBoundaryIntCells_); 
  }
  
  SanityChecks();
  
  PrintStatistics();
  
  elementsByLabel_.clear();
  CollectElementsByXFEMCouplingLabel(*cutterdis, elementsByLabel_);
  
  // invert collection
  labelPerElementId_.clear();
  XFEM::InvertElementsByLabel(elementsByLabel_, labelPerElementId_);

  const LINALG::Matrix<3,2> cutterAABB = GEO::getXAABBofDis(*cutterdis,cutterposnp_);
  const LINALG::Matrix<3,2> xfemAABB =GEO::getXAABBofDis(*xfemdis);
  const LINALG::Matrix<3,2> AABB = GEO::mergeAABB(cutterAABB, xfemAABB);
//  octTreenp_ = rcp( new GEO::SearchTree(5));
  octTreenp_->initializeTree(AABB, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  //octTreen_ = rcp( new GEO::SearchTree(5));
  octTreen_->initializeTree(AABB, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  
  // find malicious entries
  const std::set<int> ele_to_delete = FindDoubleCountedIntersectedElements();
  
  // remove malicious entries from both maps
  for (set<int>::const_iterator eleid = ele_to_delete.begin(); eleid != ele_to_delete.end(); ++eleid)
  {
    elementalDomainIntCells_.erase(*eleid);
    elementalBoundaryIntCells_.erase(*eleid);
  }
  
  GenerateSpaceTimeLayer(cutterdis_, cutterposnp_, cutterposn_);
  
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandleXFSI::~InterfaceHandleXFSI()
{
    return;
}


void XFEM::InterfaceHandleXFSI::FillCurrentCutterPositionMap(
    const Teuchos::RCP<DRT::Discretization>& cutterdis,
    const Epetra_Vector&                     idispcol,
    std::map<int,LINALG::Matrix<3,1> >&      currentcutterpositions
    ) const
{
  currentcutterpositions.clear();
  
  vector<int> lm(3);
  for (int lid = 0; lid < cutterdis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = cutterdis->lColNode(lid);
    cutterdis->Dof(node, 0, lm);
    vector<double> mydisp(3);
    DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);
    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    currentcutterpositions[node->Id()] = currpos;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::InterfaceHandleXFSI::FindDoubleCountedIntersectedElements() const
{
  // clean up double counted intersections
  std::set<int> ele_to_delete;

  // find unintersected elements and put their Id them in aboves set
  std::map<int,GEO::DomainIntCells >::const_iterator entry;
  for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
  {
    const GEO::DomainIntCells cells = entry->second;
    DRT::Element* xfemele = xfemdis_->gElement(entry->first);
    std::set<int> labelset;
    bool one_cell_is_fluid = false;
    for (GEO::DomainIntCells::const_iterator cell = cells.begin(); cell != cells.end(); ++cell)
    {
      const LINALG::Matrix<3,1> cellcenter(cell->GetPhysicalCenterPosition(*xfemele));
      const int current_label = PositionWithinConditionNP(cellcenter);
      if (current_label == 0)
      {
        one_cell_is_fluid = true;
      }
      labelset.insert(current_label);
    }
    const bool all_cells_in_same_domain = (labelset.size() == 1);
    if (all_cells_in_same_domain and not one_cell_is_fluid)
    {
      ele_to_delete.insert(xfemele->Id());
    }
  }
  return ele_to_delete;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::toGmsh(const int step) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
  
  const bool screen_out = false;
  
  const bool gmsh_tree_output = false;
  
  const int myrank = xfemdis_->Comm().MyPID();
  
  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_elements_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::XdisToString("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_);
    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, cutterposnp_);
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }
  
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_domains_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_domains_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
    {
      // stringstream for domains
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;
      
      for (int i=0; i<xfemdis_->NumMyRowElements(); ++i)
      {
        DRT::Element* actele = xfemdis_->lRowElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele->Id(), actele->Shape());
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          LINALG::SerialDenseMatrix cellpos(3,cell->NumNode()); 
          cell->NodalPosXYZ(*actele, cellpos);
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition(*actele));
          const int domain_id = PositionWithinConditionNP(cellcenterpos);
          //const double color = domain_id*100000+(closestElementId);
          const double color = domain_id;
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), color, cellpos) << endl;
        };
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }
  
  if (gmshdebugout) // print space time layer
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_spacetime_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_spacetime_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
    {
      // stringstream for domains
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "SpaceTime cells \" {" << endl;
      LINALG::SerialDenseVector vals(8);
      vals(0) = 0.0;vals(1) = 0.0;vals(2) = 0.0;vals(3) = 0.0;
      vals(4) = 1.0;vals(5) = 1.0;vals(6) = 1.0;vals(7) = 1.0;
      for (std::map<int,XFEM::SpaceTimeBoundaryCell>::const_iterator slabiter = stlayer_.begin(); slabiter != stlayer_.end(); ++slabiter)
      {
        const XFEM::SpaceTimeBoundaryCell& slabitem = slabiter->second;
        
        gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(DRT::Element::hex8, vals, slabitem.get_xyzt()) << endl;
      }
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }
  
  
  if (gmsh_tree_output)
  {
    // debug: write information about which structure we are in
    std::stringstream filenameP;
    std::stringstream filenamePdel;
    filenameP    << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_points_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamePdel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_points_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamePdel.str().c_str());

    std::cout << "writing " << left << std::setw(50) <<filenameP.str()<<"...";
    std::ofstream f_systemP(filenameP.str().c_str());
    {
      // stringstream for cellcenter points
      stringstream gmshfilecontentP;
      gmshfilecontentP << "View \" " << "CellCenter of Elements and Integration Cells \" {" << endl;
     
      for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
      {
        DRT::Element* actele = xfemdis_->lColElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele->Id(), actele->Shape());
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          LINALG::SerialDenseMatrix cellpos(3,cell->NumNode()); 
          cell->NodalPosXYZ(*actele, cellpos);
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition(*actele));
          
          //const int domain_id = PositionWithinConditionNP(cellcenterpos);
          
          LINALG::SerialDenseMatrix point(3,1);
          point(0,0)=cellcenterpos(0);
          point(1,0)=cellcenterpos(1);
          point(2,0)=cellcenterpos(2);

          gmshfilecontentP << IO::GMSH::cellWithScalarToString(DRT::Element::point1, (actele->Id()), point) << endl;              
        };
      };
      gmshfilecontentP << "};" << endl;
      f_systemP << gmshfilecontentP.str();
    }
    f_systemP.close();
    cout << " done" << endl;
    
    octTreenp_->printTree(DRT::Problem::Instance()->OutputControlFile()->FileName(), step);
    octTreenp_->evaluateTreeMetrics(step);
  }
  
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::InterfaceHandleXFSI::FindSpaceTimeLayerCell(
    const LINALG::Matrix<3,1>&        querypos,
    XFEM::SpaceTimeBoundaryCell&      stcell,
    LINALG::Matrix<3,1>&              rst
) const
{
  bool in_spacetimecell = false;
  // loop over space time slab cells until one is found - brute force
  for (std::map<int,XFEM::SpaceTimeBoundaryCell>::const_iterator slabiter = stlayer_.begin(); slabiter != stlayer_.end(); ++slabiter)
  {
    const XFEM::SpaceTimeBoundaryCell slabitem = slabiter->second;
    LINALG::Matrix<3,1> xsi(true);

    in_spacetimecell = GEO::currentToVolumeElementCoordinates(DRT::Element::hex8, slabitem.get_xyzt(), querypos, xsi);
    in_spacetimecell = GEO::checkPositionWithinElementParameterSpace(xsi, DRT::Element::hex8);
    if (in_spacetimecell)
    {
      stcell = XFEM::SpaceTimeBoundaryCell(slabitem);
//      cout << "slabitem " << slabitem.toString() << endl;
      rst = xsi;    
      return true;
    }
  }
  
//  if (not in_spacetimecell)
//  {
//    dserror("should be in one space time cell");
//  }
//  cout << stcell.toString() << endl;
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionNP(
    const LINALG::Matrix<3,1>&        x_in
) const
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionNP");
  return octTreenp_->queryXFEMFSIPointType(*(cutterdis_), cutterposnp_, x_in);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionN(
    const LINALG::Matrix<3,1>&        x_in
) const
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionN");
  return octTreen_->queryXFEMFSIPointType(*(cutterdis_), cutterposn_, x_in);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionNP(
    const LINALG::Matrix<3,1>&        x_in,
    GEO::NearestObject&               nearestobject
) const
{
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionNP");
  return octTreenp_->queryFSINearestObject(*(cutterdis_), cutterposnp_, x_in, nearestobject);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionN(
    const LINALG::Matrix<3,1>&        x_in,
    GEO::NearestObject&               nearestobject
) const
{
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionN");
  return octTreen_->queryFSINearestObject(*(cutterdis_), cutterposn_, x_in, nearestobject);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::GenerateSpaceTimeLayer(
    const Teuchos::RCP<DRT::Discretization>&            cutterdis,
    const std::map<int,LINALG::Matrix<3,1> >&           cutterposnp,
    const std::map<int,LINALG::Matrix<3,1> >&           cutterposn
)
{
 
  for (int i=0; i<cutterdis->NumMyColElements(); ++i)
  {
    const DRT::Element* cutterele = cutterdis->lColElement(i);
    const int* nodeids = cutterele->NodeIds();
    LINALG::SerialDenseMatrix posnp(3,4);
    for (int inode = 0; inode != 4; ++inode) // fill n+1 position
    {
      const int nodeid = nodeids[inode];
      const LINALG::Matrix<3,1> nodexyz = cutterposnp.find(nodeid)->second;
      for (int isd = 0; isd != 3; ++isd)
      {
        posnp(isd,inode) = nodexyz(isd);
      }
    }
    LINALG::SerialDenseMatrix posn(3,4);
    for (int inode = 0; inode != 4; ++inode) // fill n   position
    {
      const int nodeid = nodeids[inode];
      const LINALG::Matrix<3,1> nodexyz = cutterposn.find(nodeid)->second;
      for (int isd = 0; isd != 3; ++isd)
      {
        posn(isd,inode) = nodexyz(isd);
      }
    }
    //XFEM::SpaceTimeBoundaryCell slab(cutterele->Id(),posnp,posn);
    stlayer_.insert(make_pair(cutterele->Id(),XFEM::SpaceTimeBoundaryCell(cutterele->Id(),posnp,posn)));
    //cout << "XFEM::SpaceTimeBoundaryCell" << slab.getBeleId() << endl;
  }
  
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::PrintStatistics() const
{

  // loop intersected elements and count intcells
  const unsigned numintersectedele = elementalDomainIntCells_.size();
  
  if (numintersectedele > 0)
  {
    unsigned numcells = 0;
    std::map<int,GEO::DomainIntCells >::const_iterator entry;
    for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
    {
      const GEO::DomainIntCells cells = entry->second;
      numcells += cells.size();
    }
    const unsigned avgnumcellperele = numcells/numintersectedele;
    cout << "Avg. Number of DomainIntCells per intersected xfem element: " << avgnumcellperele << endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::CheckSurfaceElementSize() const
{
  cout << "Warning: Not implemented, yet! No checking performed!" << endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::CheckXFEMElementSize() const
{
  
  int small_ele_counter = 0;
  double smallest_edge_length = 1.0e12;
  
  const double smallest_allowed_edgelength = 1.0e-2;
  
  for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
  {
    const DRT::Element* xfemele = xfemdis_->lColElement(i);
    const double elevol = XFEM::ElementVolume(*xfemele);
    
    // assume that vol = h^3 for a cube
    const double characteristic_edge_length = std::pow(elevol,1.0/3.0);
//    cout << "characteristic_edge_length: " << characteristic_edge_length << endl;
    
    smallest_edge_length = std::min(smallest_edge_length,characteristic_edge_length);
    // count to small edges
    if (characteristic_edge_length < smallest_allowed_edgelength)
    {
      small_ele_counter++;
    }
  }
  
  if (small_ele_counter > 0)
  {
    cout << "Warning: Smallest edge length is " << smallest_edge_length << ", ";
    cout << small_ele_counter << " elements are too small!" << endl;
  }
}

/*! 
 * \brief calculates the volume of an element in initial configuration
 * 
 * \return physical volume of the element in initial configuration (without any displacement)
 */
template <DRT::Element::DiscretizationType DISTYPE>
double ElementVolumeT(
        const DRT::Element&           ele
        )
{
  if (ele.Shape() != DISTYPE)
  {
    dserror("mismatch between element shape and template parameter DISTYPE! This is a bug!");
  }

  // number of nodes for element
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  
  // dimension for 3d element
  const int nsd = 3;
  
  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(&ele, xyze);
  
  // physical element volume 
  double vol = 0.0;
  
  DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
  switch (DISTYPE)
  {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      gaussrule = DRT::UTILS::intrule_hex_8point;
      break;
    }
    case DRT::Element::tet4:
    case DRT::Element::tet10:
    {
      gaussrule = DRT::UTILS::intrule_tet_4point;
      break;
    }
    default:
      dserror("add your element type here...");
  }
  
  // integration points
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point in element coordinates \xi
    static LINALG::Matrix<3,1> posXiDomain;
    posXiDomain(0) = intpoints.qxg[iquad][0];
    posXiDomain(1) = intpoints.qxg[iquad][1];
    posXiDomain(2) = intpoints.qxg[iquad][2];
    
    // shape functions and their first derivatives
    static LINALG::Matrix<nsd,numnode> deriv;
    DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    // get transposed of the jacobian matrix d x / d \xi
    static LINALG::Matrix<3,3> xjm;
    xjm.MultiplyNT(deriv,xyze);

    const double det = xjm.Determinant();
    const double fac = intpoints.qwgt[iquad]*det;

    if (det <= 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele.Id(), det);
    }
   
    vol += fac;
      
  }
  return vol;
}

double XFEM::ElementVolume(
        const DRT::Element&           ele
        )
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return ElementVolumeT<DRT::Element::hex8>(ele);
    case DRT::Element::hex20:
      return ElementVolumeT<DRT::Element::hex20>(ele);
    case DRT::Element::hex27:
      return ElementVolumeT<DRT::Element::hex27>(ele);
    case DRT::Element::tet4:
      return ElementVolumeT<DRT::Element::tet4>(ele);
    case DRT::Element::tet10:
      return ElementVolumeT<DRT::Element::tet10>(ele);
    default:
      dserror("add you distype here...");
      exit(1);
  }
}

#endif  // #ifdef CCADISCRET
