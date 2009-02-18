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
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_gmsh_xfem_extension.H"
#include "../drt_geometry/intersection.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "enrichment_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"



/*----------------------------------------------------------------------*
 * standard constructor
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
  currentXAABBs_ = GEO::getCurrentXAABBs(*cutterdis, cutterposnp_);
  
  static std::map<int,LINALG::Matrix<3,2> > currentXAABBs; 
  currentXAABBs = GEO::getCurrentXAABBs(*cutterdis, cutterposnp_);
  
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = false;
  
  const int myrank = xfemdis_->Comm().MyPID();

  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".uncut_elements_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".uncut_elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis_);
    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, cutterposnp_);
    f_system.close();
    if (screen_out) cout << " done" << endl;
  }

  GEO::Intersection is;
  is.computeIntersection(xfemdis, cutterdis, cutterposnp_, currentXAABBs, elementalDomainIntCells_, elementalBoundaryIntCells_);  

  xfemdis->Comm().Barrier();
  
  EraseTinyDomainIntCells(elementalDomainIntCells_);
  
  xfemdis->Comm().Barrier();
   
  PrintStatistics();
  
  // collect elements by xfem label for tree and invert
  DRT::UTILS::CollectElementsByConditionLabel(*cutterdis, elementsByLabel_, "XFEMCoupling");
  InvertElementsPerLabel();
  
  const LINALG::Matrix<3,2> cutterAABB = GEO::getXAABBofDis(*cutterdis,cutterposnp_);
  const LINALG::Matrix<3,2> xfemAABB =GEO::getXAABBofDis(*xfemdis);
  const LINALG::Matrix<3,2> AABB = GEO::mergeAABB(cutterAABB, xfemAABB);

  octTreenp_->initializeTree(AABB, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  octTreen_->initializeTree(AABB, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  
  ClassifyIntegrationCells();
 
  GenerateSpaceTimeLayer(cutterdis_, cutterposnp_, cutterposn_);

  xfemdis->Comm().Barrier();
}

    

/*----------------------------------------------------------------------*
 * destructor
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandleXFSI::~InterfaceHandleXFSI()
{
    return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
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
void XFEM::InterfaceHandleXFSI::ClassifyIntegrationCells()
{
  // clean up double counted intersections
  std::set<int> ele_to_delete;

  // find unintersected elements and put their Id them in aboves set
  std::map<int,GEO::DomainIntCells >::const_iterator entry;
  for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
  {
    GEO::DomainIntCells cells = entry->second;
    const DRT::Element* xfemele = xfemdis_->gElement(entry->first);
    std::set<int> labelset;
    bool one_cell_is_fluid = false;
    for (GEO::DomainIntCells::iterator cell = cells.begin(); cell != cells.end(); ++cell)
    {
      const LINALG::Matrix<3,1> cellcenter(cell->GetPhysicalCenterPosition());
      const int current_label = PositionWithinConditionNP(cellcenter);
      cell->setLabel(current_label);
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
  
  // remove malicious entries from both maps
  for (std::set<int>::const_iterator eleid = ele_to_delete.begin(); eleid != ele_to_delete.end(); ++eleid)
  {
    elementalDomainIntCells_.erase(*eleid);
    elementalBoundaryIntCells_.erase(*eleid);
  }
 
  return;
}

/*----------------------------------------------------------------------*
 * delete tiny cells, if they are only a small fraction
 * of the parent element
 * 
 * this is purely for speedup, conditioning of the tangent stiffness
 * is not affected 
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::EraseTinyDomainIntCells(
    std::map<int,GEO::DomainIntCells >& elementalDomainIntCells) const
{
  
  unsigned small_cell_count = 0;
  const double small_cell_treshold = 1.0e-12;
  
  // check size of each cell and delete, if to small
  // implementation copies each cell ones!!!
  // not very efficient - improve such that no copy is necessary
  // (can't find command to erase from vector without killing the iterator)
  std::map<int,GEO::DomainIntCells >::iterator entry;
  for (entry = elementalDomainIntCells.begin(); entry != elementalDomainIntCells.end(); ++entry)
  {
    const GEO::DomainIntCells cells_old = entry->second;
    GEO::DomainIntCells& cells = entry->second;
    cells.clear();
    
    const DRT::Element* xfemele = xfemdis_->gElement(entry->first);
    for (GEO::DomainIntCells::const_iterator cell = cells_old.begin(); cell != cells_old.end(); ++cell)
    {
      const double size = cell->SizeXiDomain(*xfemele);
      if (size < small_cell_treshold)
      {
//        cout << RED << size << END_COLOR << endl;
        small_cell_count++;
      }
      else
      {
        cells.push_back(*cell);
//        cout << GREEN << size << END_COLOR << endl;
      }
    }
  }
  cout << " " << small_cell_count << " small cells ( v_ele / v_cell < " << small_cell_treshold << " % ) deleted." << endl;
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionNP(
    const LINALG::Matrix<3,1>&        x_in) const
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionNP");
  return octTreenp_->queryXFEMFSIPointType(*(cutterdis_), cutterposnp_, currentXAABBs_, x_in);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionN(
    const LINALG::Matrix<3,1>&        x_in) const
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionN");
  return octTreen_->queryXFEMFSIPointType(*(cutterdis_), cutterposn_, currentXAABBs_, x_in);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionNP(
    const LINALG::Matrix<3,1>&        x_in,
    GEO::NearestObject&               nearestobject) const
{
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionNP");
  return octTreenp_->queryFSINearestObject(*(cutterdis_), cutterposnp_, currentXAABBs_, x_in, nearestobject);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithinConditionN(
    const LINALG::Matrix<3,1>&        x_in,
    GEO::NearestObject&               nearestobject) const
{
  TEUCHOS_FUNC_TIME_MONITOR(" - search - InterfaceHandle::PositionWithinConditionN");
  return octTreen_->queryFSINearestObject(*(cutterdis_), cutterposn_, currentXAABBs_, x_in, nearestobject);
}


/*----------------------------------------------------------------------*
 * generate space time layer
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::GenerateSpaceTimeLayer(
    const Teuchos::RCP<DRT::Discretization>&            cutterdis,
    const std::map<int,LINALG::Matrix<3,1> >&           cutterposnp,
    const std::map<int,LINALG::Matrix<3,1> >&           cutterposn)
{
 
  for (int i=0; i<cutterdis->NumMyColElements(); ++i)
  {
    const DRT::Element* cutterele = cutterdis->lColElement(i);
    const int* nodeids = cutterele->NodeIds();
    const int numNode = cutterele->NumNode();
    LINALG::SerialDenseMatrix posnp(3, numNode);
    for (int inode = 0; inode != numNode; ++inode) // fill n+1 position
    {
      const int nodeid = nodeids[inode];
      const LINALG::Matrix<3,1> nodexyz = cutterposnp.find(nodeid)->second;
      for (int isd = 0; isd != 3; ++isd)
      {
        posnp(isd,inode) = nodexyz(isd);
      }
    }
    LINALG::SerialDenseMatrix posn(3,numNode);
    for (int inode = 0; inode != numNode; ++inode) // fill n   position
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
 * find space time layer
 *----------------------------------------------------------------------*/
bool XFEM::InterfaceHandleXFSI::FindSpaceTimeLayerCell(
    const LINALG::Matrix<3,1>&        querypos,
    XFEM::SpaceTimeBoundaryCell&      stcell,
    LINALG::Matrix<3,1>&              rst) const
{
  bool in_spacetimecell = false;
  // loop over space time slab cells until one is found - brute force
  for (std::map<int,XFEM::SpaceTimeBoundaryCell>::const_iterator slabiter = stlayer_.begin(); slabiter != stlayer_.end(); ++slabiter)
  {
    const XFEM::SpaceTimeBoundaryCell slabitem = slabiter->second;
    LINALG::Matrix<3,1> xsi(true);

    // TODO check if hardcoding hex8 is ok (needs more theoretical work)
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
 * Debug only
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
 * Debug only
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
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
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
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".domains_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
    {
      // stringstream for domains
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;
      
      for (int i=0; i<xfemdis_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = xfemdis_->lRowElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::SerialDenseMatrix& cellpos = cell->CellNodalPosXYZ();
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
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
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".spacetime_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".spacetime_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
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
        const DRT::Element* actele = xfemdis_->lColElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
          
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



#endif  // #ifdef CCADISCRET
