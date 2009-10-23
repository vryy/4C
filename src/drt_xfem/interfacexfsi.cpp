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
#include "../drt_geometry/intersection_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "enrichment_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"



/*----------------------------------------------------------------------*
 * standard constructor
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandleXFSI::InterfaceHandleXFSI(
    const Teuchos::RCP<DRT::Discretization>  xfemdis,
    const Teuchos::RCP<DRT::Discretization>  cutterdis
    ) : InterfaceHandle(xfemdis),
        cutterdis_(cutterdis)
{
  if (cutterdis == Teuchos::null)
    dserror("We need a real boundary discretization here!");

  if (xfemdis->Comm().MyPID() == 0)
    std::cout << "Constructing InterfaceHandle" << std::endl;

  FillCurrentCutterPositionMap(cutterdis, *cutterdis->GetState("idispcolnp"), cutterposnp_);
  FillCurrentCutterPositionMap(cutterdis, *cutterdis->GetState("idispcoln") , cutterposn_ );
  currentXAABBs_ = GEO::getCurrentXAABBs(*cutterdis, cutterposnp_);

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = false;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("uncut_elements_coupled_system", 999999, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
    IO::GMSH::disToStream("Fluid", 0.0, xfemdis_, gmshfilecontent);
    IO::GMSH::disToStream("Boundary", 1.0, cutterdis_, cutterposnp_, gmshfilecontent);
    gmshfilecontent.close();
    if (screen_out) cout << " done" << endl;
  }

  // call before intersection
  // collect elements by xfem label for tree and invert
  DRT::UTILS::CollectElementsByConditionLabel(*cutterdis, boundaryElementsByLabel_, "XFEMCoupling");
  InvertElementsPerLabel();

  if(cutterdis_->NumMyColElements()!=0)
  {
#ifdef QHULL
    Teuchos::RCP<GEO::Intersection> is = rcp(new GEO::Intersection());
    is->computeIntersection(xfemdis, cutterdis, cutterposnp_, currentXAABBs_, elementalDomainIntCells_, elementalBoundaryIntCells_, labelPerBoundaryElementId_);
#else
    dserror("you have to compile with the QHULL flag!");
#endif
  }


  xfemdis->Comm().Barrier();

  PrintStatistics();

  const LINALG::Matrix<3,2> cutterAABBnp = GEO::getXAABBofDis(*cutterdis,cutterposnp_);
  const LINALG::Matrix<3,2> cutterAABBn = GEO::getXAABBofDis(*cutterdis,cutterposn_);
  const LINALG::Matrix<3,2> cutterAABB = GEO::mergeAABB(cutterAABBn, cutterAABBnp);
  const LINALG::Matrix<3,2> xfemAABB = GEO::getXAABBofDis(*xfemdis);
  const LINALG::Matrix<3,2> AABB = GEO::mergeAABB(cutterAABB, xfemAABB);

  octTreenp_->initializeTree(AABB, boundaryElementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  octTreen_->initializeTree(AABB, boundaryElementsByLabel_, GEO::TreeType(GEO::OCTTREE));

  for (std::map<int,std::set<int> >::const_iterator entry = boundaryElementsByLabel_.begin();
      entry != boundaryElementsByLabel_.end();
      ++entry)
  {
    const int label = entry->first;
    map<int,set<int> > onelabelmap;
    onelabelmap[label] = entry->second;
    octTreePerLabelnp_.insert(make_pair(label,rcp( new GEO::SearchTree(20))));
    octTreePerLabelnp_[label]->initializeTree(AABB, onelabelmap, GEO::TreeType(GEO::OCTTREE));
    octTreePerLabeln_.insert(make_pair(label,rcp( new GEO::SearchTree(20))));
    octTreePerLabeln_[label]->initializeTree(AABB, onelabelmap, GEO::TreeType(GEO::OCTTREE));
  }

//  ClassifyIntegrationCells();

  GenerateSpaceTimeLayer(cutterdis_, cutterposnp_, cutterposn_);

  xfemdis->Comm().Barrier();
  if (xfemdis->Comm().MyPID() == 0)
    cout << "Interfacehandle constructed" << endl;
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

  for (int lid = 0; lid < cutterdis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = cutterdis->lColNode(lid);
    vector<int> lm;
    cutterdis->Dof(node, lm);
    vector<double> mydisp;
    DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);
    if (mydisp.size() != 3)
      dserror("we need 3 displacements here");

    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    currentcutterpositions.insert(make_pair(node->Id(),currpos));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::ClassifyIntegrationCells()
{
  // clean up double counted intersections
  std::set<int> ele_to_delete;

  // find uncut elements and put their Id them in aboves set
  std::map<int,GEO::DomainIntCells >::iterator entry;
  for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
  {
    GEO::DomainIntCells& cells = entry->second;
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
 * check if intcell fill their xfem ele                     u.may 04/09 *
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandleXFSI::TestDomainIntCells() const
{
  std::map<int,GEO::DomainIntCells >::const_iterator entry;
  for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
  {
    const GEO::DomainIntCells cells = entry->second;
    const DRT::Element* xfemele = xfemdis_->gElement(entry->first);
    double cellFillFactor = 0.0;

    for (GEO::DomainIntCells::const_iterator cell = cells.begin(); cell != cells.end(); ++cell)
      cellFillFactor += cell->VolumeInXiDomain(*xfemele);

    if (std::abs(1.0 - cellFillFactor) > 1.0e-7)
    {
      cout << "unit element volume integrated by using integration cells does not sum up to 1" << endl;
      cout << scientific << (1.0-cellFillFactor) << endl;
      dserror("");
    }
  }
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
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithRespectToInterfaceNP(
    const LINALG::Matrix<3,1>&        x_in,
    const int label) const
{
  return octTreePerLabelnp_.find(label)->second->queryXFEMFSIPointType(*(cutterdis_), cutterposnp_, currentXAABBs_, x_in);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::InterfaceHandleXFSI::PositionWithRespectToInterfaceN(
    const LINALG::Matrix<3,1>&        x_in,
    const int label) const
{
  return octTreePerLabeln_.find(label)->second->queryXFEMFSIPointType(*(cutterdis_), cutterposn_, currentXAABBs_, x_in);
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
    stlayer_.insert(make_pair(cutterele->Id(),XFEM::SpaceTimeBoundaryCell(cutterele->Id(),cutterele->Shape(),posnp,posn)));
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
  const std::size_t numintersectedele = elementalDomainIntCells_.size();

  if (numintersectedele > 0)
  {
    std::size_t numcells = 0;
    std::map<int,GEO::DomainIntCells >::const_iterator entry;
    for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
    {
      const GEO::DomainIntCells& cells = entry->second;
      numcells += cells.size();
    }
    const std::size_t avgnumcellperele = numcells/numintersectedele;
    cout << " Avg. Number of DomainIntCells per intersected xfem element: " << avgnumcellperele << endl;
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

  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("elements_coupled_system", step, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
    IO::GMSH::XdisToStream("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_, gmshfilecontent);
    IO::GMSH::disToStream("Solid", 1.0, cutterdis_, cutterposnp_, gmshfilecontent);
    gmshfilecontent.close();
    if (screen_out) cout << " done" << endl;
  }

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("interface_patches", step, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
//    IO::GMSH::XdisToStream("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_, gmshfilecontent);

    const std::string s("Interface patches");
    gmshfilecontent << "View \" " << s << " Elements \" {\n";
    for (int i=0; i<cutterdis_->NumMyColElements(); ++i)
    {
      const DRT::Element* actele = cutterdis_->lColElement(i);
      const int scalar = this->GetLabelPerBoundaryElementId(actele->Id());
      IO::GMSH::cellWithScalarToStream(actele->Shape(),
          scalar, GEO::getCurrentNodalPositions(actele,cutterposnp_), gmshfilecontent);
    };
    gmshfilecontent << "};\n";

    gmshfilecontent.close();
    if (screen_out) cout << " done" << endl;
  }

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("domains", step, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
    {
      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;

      for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = xfemdis_->lColElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::SerialDenseMatrix& cellpos = cell->CellNodalPosXYZ();
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
          const int domain_id = PositionWithinConditionNP(cellcenterpos);
          const double color = domain_id;
          IO::GMSH::cellWithScalarToStream(cell->Shape(), color, cellpos, gmshfilecontent);
        };
      };
      gmshfilecontent << "};" << endl;
    }
    gmshfilecontent.close();
    if (screen_out) cout << " done" << endl;
  }

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("domains_for_patches", step, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
    for (map<int, const Teuchos::RCP<GEO::SearchTree> >::const_iterator entry = octTreePerLabelnp_.begin();
        entry != octTreePerLabelnp_.end();
        ++entry)
    {
      const int label = entry->first;
      const Teuchos::RCP<GEO::SearchTree> tree = entry->second;


      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;

      for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = xfemdis_->lColElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::SerialDenseMatrix& cellpos = cell->CellNodalPosXYZ();
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());
          const int domain_id = octTreePerLabelnp_.find(label)->second->queryXFEMFSIPointType(*(cutterdis_), cutterposnp_, currentXAABBs_, cellcenterpos);
          //const double color = domain_id*100000+(closestElementId);
          const double color = domain_id;
          IO::GMSH::cellWithScalarToStream(cell->Shape(), color, cellpos, gmshfilecontent);
        };
      };
      gmshfilecontent << "};" << endl;
    }
    gmshfilecontent.close();
    if (screen_out) cout << " done" << endl;
  }

  if (gmshdebugout) // print space time layer
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("spacetime", step, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
    {
      gmshfilecontent << "View \" " << "SpaceTime cells \" {" << endl;
      LINALG::SerialDenseVector vals(8);
      vals(0) = 0.0;vals(1) = 0.0;vals(2) = 0.0;vals(3) = 0.0;
      vals(4) = 1.0;vals(5) = 1.0;vals(6) = 1.0;vals(7) = 1.0;
      for (std::map<int,XFEM::SpaceTimeBoundaryCell>::const_iterator slabiter = stlayer_.begin(); slabiter != stlayer_.end(); ++slabiter)
      {
        const XFEM::SpaceTimeBoundaryCell& slabitem = slabiter->second;

        IO::GMSH::cellWithScalarFieldToStream(DRT::Element::hex8, vals, slabitem.get_xyzt(), gmshfilecontent);
      }
      gmshfilecontent << "};" << endl;
    }
    gmshfilecontent.close();
    if (screen_out) cout << " done" << endl;
  }


  if (gmsh_tree_output)
  {
    // debug: write information about which structure we are in
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("points", step, 5, screen_out);
    std::ofstream gmshfilecontent(filename.c_str());
    {
      gmshfilecontent << "View \" " << "CellCenter of Elements and Integration Cells \" {" << endl;

      for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = xfemdis_->lColElement(i);
        const GEO::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele);
        GEO::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          const LINALG::Matrix<3,1> cellcenterpos(cell->GetPhysicalCenterPosition());

          LINALG::SerialDenseMatrix point(3,1);
          point(0,0)=cellcenterpos(0);
          point(1,0)=cellcenterpos(1);
          point(2,0)=cellcenterpos(2);

          IO::GMSH::cellWithScalarToStream(DRT::Element::point1, (actele->Id()), point, gmshfilecontent);
        };
      };
      gmshfilecontent << "};" << endl;
    }
    gmshfilecontent.close();
    cout << " done" << endl;

    octTreenp_->printTree(DRT::Problem::Instance()->OutputControlFile()->FileName(), step);
    octTreenp_->evaluateTreeMetrics(step);
  }

  return;
}



#endif  // #ifdef CCADISCRET
