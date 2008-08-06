/*!
\file interface.cpp

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include "interface.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "xfsi_searchtree.H"

#include "xfem_condition.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_gmsh_xfem_extension.H"
#include "integrationcell.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::InterfaceHandle(
    const Teuchos::RCP<DRT::Discretization>  xfemdis, 
    const Teuchos::RCP<DRT::Discretization>  cutterdis,
    const Epetra_Vector&                     idispcol) :
      xfemdis_(xfemdis),
      cutterdis_(cutterdis)
{
  FillCurrentCutterPositionMap(cutterdis, idispcol, currentcutterpositions_);
  
  elementalDomainIntCells_.clear();
  elementalBoundaryIntCells_.clear();
  XFEM::Intersection is;
  is.computeIntersection(xfemdis, cutterdis, currentcutterpositions_,elementalDomainIntCells_, elementalBoundaryIntCells_);
  
//  std::cout << "numcuttedelements (elementalDomainIntCells_)   = " << elementalDomainIntCells_.size() << endl;
//  std::cout << "numcuttedelements (elementalBoundaryIntCells_) = " << elementalBoundaryIntCells_.size() << endl;
  if (elementalDomainIntCells_.size() != elementalBoundaryIntCells_.size())
  {
    dserror("mismatch in cutted elements maps!");  
  }
  
  // sanity check, whether, we really have integration cells in the map
  for (std::map<int,DomainIntCells>::const_iterator 
      tmp = elementalDomainIntCells_.begin();
      tmp != elementalDomainIntCells_.end();
      ++tmp)
  {
    dsassert(tmp->second.empty() == false, "this is a bug!");
  }
  
  // sanity check, whether, we really have integration cells in the map
  for (std::map<int,BoundaryIntCells>::const_iterator 
      tmp = elementalBoundaryIntCells_.begin();
      tmp != elementalBoundaryIntCells_.end();
      ++tmp)
  {
    dsassert(tmp->second.empty() == false, "this is a bug!");
  }
  
  elementsByLabel_.clear();
  CollectElementsByXFEMCouplingLabel(*cutterdis, elementsByLabel_);

  
  const BlitzMat3x2 cutterAABB = XFEM::getXAABBofDis(*cutterdis,currentcutterpositions_);
  const BlitzMat3x2 xfemAABB = XFEM::getXAABBofDis(*xfemdis);
  const BlitzMat3x2 AABB = XFEM::mergeAABB(cutterAABB, xfemAABB);
  const int max_treedepth = 5;
  octTree_ = rcp( new GEO::OctTree(max_treedepth, AABB ) );
  octTree_->initializeTree(elementsByLabel_); 
 
  // find malicious entries
  const std::set<int> ele_to_delete = FindDoubleCountedIntersectedElements();
  
  // remove malicious entries from both maps
  for (set<int>::const_iterator eleid = ele_to_delete.begin(); eleid != ele_to_delete.end(); ++eleid)
  {
    elementalDomainIntCells_.erase(*eleid);
    elementalBoundaryIntCells_.erase(*eleid);
  }
  
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::~InterfaceHandle()
{
    return;
}

void XFEM::InterfaceHandle::FillCurrentCutterPositionMap(
    const Teuchos::RCP<DRT::Discretization>  cutterdis,
    const Epetra_Vector&                     idispcol,
    std::map<int,BlitzVec3>&                 currentcutterpositions
    ) const
{
  currentcutterpositions.clear();
  
  for (int lid = 0; lid < cutterdis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = cutterdis->lColNode(lid);
    vector<int> lm;
    lm.reserve(3);
    cutterdis->Dof(node, lm);
    vector<double> mydisp(3);
    DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);
    BlitzVec3 currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    currentcutterpositions[node->Id()] = currpos;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::InterfaceHandle::FindDoubleCountedIntersectedElements() const
{
  // clean up double counted intersections
  std::set<int> ele_to_delete;

  // find unintersected elements and put their Id them in aboves set
  std::map<int, DomainIntCells >::const_iterator entry;
  for (entry = elementalDomainIntCells_.begin(); entry != elementalDomainIntCells_.end(); ++entry)
  {
    const DomainIntCells cells = entry->second;
    DRT::Element* xfemele = xfemdis_->gElement(entry->first);
    DomainIntCells::const_iterator cell;
    std::set<int> labelset;
    bool one_cell_is_fluid = false;
    for (cell = cells.begin(); cell != cells.end(); ++cell)
    {
      const BlitzVec3 cellcenter(cell->GetPhysicalCenterPosition(*xfemele));
      const int current_label = PositionWithinCondition(cellcenter, *this);
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
std::string XFEM::InterfaceHandle::toString() const
{
	std::stringstream s(" ");
	return s.str();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandle::toGmsh(const int step) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
  
  const bool gmsh_tree_output = false;
  
  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    std::stringstream filenamedel;
    filename << allfiles.outputfile_kenner << "_elements_coupled_system_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamedel << allfiles.outputfile_kenner << "_elements_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::XdisToString("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_);
    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, currentcutterpositions_);
    f_system.close();
    cout << " done" << endl;
  }
  
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename << allfiles.outputfile_kenner << "_domains_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamedel << allfiles.outputfile_kenner << "_domains_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
    {
      // stringstream for domains
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;
      
      for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
      {
        DRT::Element* actele = xfemdis_->lColElement(i);
        const XFEM::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele->Id(), actele->Shape());
        XFEM::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          
          BlitzMat cellpos(3,cell->NumNode()); 
          cell->NodalPosXYZ(*actele, cellpos);
          const BlitzVec3 cellcenterpos(cell->GetPhysicalCenterPosition(*actele));
          const int domain_id = PositionWithinCondition(cellcenterpos, *this);
          //const double color = domain_id*100000+(closestElementId);
          const double color = domain_id;
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), color, cellpos) << endl;
        };
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    cout << " done" << endl;
  }
  if (gmsh_tree_output)
  {
    // debug: write information about which structure we are in
    std::stringstream filenameP;
    std::stringstream filenamePdel;
    filenameP << allfiles.outputfile_kenner << "_points_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamePdel << allfiles.outputfile_kenner << "_points_" << std::setw(5) << setfill('0') << step-5 << ".pos";
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
        const XFEM::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele->Id(), actele->Shape());
        XFEM::DomainIntCells::const_iterator cell;
        for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
        {
          
          BlitzMat cellpos(3,cell->NumNode()); 
          cell->NodalPosXYZ(*actele, cellpos);
          const BlitzVec3 cellcenterpos(cell->GetPhysicalCenterPosition(*actele));
          
          //const int domain_id = PositionWithinCondition(cellcenterpos, *this);
          
          BlitzMat point(3,1);
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
    
    octTree_->printTree(allfiles.outputfile_kenner, step);
  }
  // TODO implement for octtree xTree_->printTreeMetrics(step);
  // TODO implement for octtree xTree_->printTreeMetricsFile(step);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::DomainIntCells XFEM::InterfaceHandle::GetDomainIntCells(
    const int gid,
    const DRT::Element::DiscretizationType distype
) const
{
  std::map<int,DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(gid);
  if (tmp == elementalDomainIntCells_.end())
  {   
    // create default set with one dummy DomainIntCell of proper size
    XFEM::DomainIntCells cells;
    cells.push_back(XFEM::DomainIntCell(distype));
    return cells;
  }
  return tmp->second;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::BoundaryIntCells XFEM::InterfaceHandle::GetBoundaryIntCells(
    const int gid
) const
{
  std::map<int,XFEM::BoundaryIntCells>::const_iterator tmp = elementalBoundaryIntCells_.find(gid);
  if (tmp == elementalBoundaryIntCells_.end())
  {   
    // return empty list
    return XFEM::BoundaryIntCells();
  }
  return tmp->second;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::InterfaceHandle::ElementIntersected(
    const int element_gid
) const
{
  if (elementalDomainIntCells_.find(element_gid) == elementalDomainIntCells_.end())
  {   
    return false;
  }
  else
  {
    return true;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinCondition(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih
)
{
  return PositionWithinConditionOctTree(x_in, ih);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinConditionOctTree(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - search - PositionWithinConditionOctTree");
  Teuchos::RCP<GEO::OctTree> xt = ih.getOctTree(); // pointer is constant, object is not!
  const int XFEMlabel = xt->queryXFEMFSIPointType(*ih.cutterdis() , *ih.currentcutterpositions(), x_in);

  // for parallel : there is nothing to be extended for parallel execution

  return XFEMlabel;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinConditionBruteForce(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih
)
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - PositionWithinCondition");
  //init
  std::map<int,bool>  posInCondition; // not really needed, but this method could go away soon, so it won't be cleaned
  const std::map<int,std::set<int> >& elementsByLabel = *(ih.elementsByLabel());
  
  /////////////////
  // loop labels
  /////////////////
  int label = 0;
  for(std::map<int,std::set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    label = conditer->first;
    posInCondition[label] = false; 
    
    // point lies opposite to a element (basis point within element parameter space)
    // works only, if I can loop over ALL surface elements
    // MUST be modified, if only a subset of the surface is used
    bool in_element = false;
    double min_ele_distance = 1.0e12;
    const DRT::Element* closest_element;
    for (set<int>::const_iterator elegid = conditer->second.begin(); elegid != conditer->second.end(); ++elegid)
    {
      const DRT::Element* cutterele = ih.cutterdis()->gElement(*elegid);
      const BlitzMat xyze_cutter(DRT::UTILS::getCurrentNodalPositions(cutterele, *ih.currentcutterpositions()));
      double distance = 0.0;
      BlitzVec2 eleCoord;
      BlitzVec3 normal;
      in_element = searchForNearestPointOnSurface(cutterele,xyze_cutter,x_in,eleCoord,normal,distance);
      if (in_element)
      {
        if (abs(distance) < abs(min_ele_distance))
        {
          closest_element = cutterele;
          min_ele_distance = distance;
        }
      }
    }
    
    if (in_element)
    {
      if (min_ele_distance < 0.0)
      {
        posInCondition[label] = true;
        break;
      }
    }
    
  } // end loop label
  
  // TODO: in parallel, we have to ask all processors, whether there is any match!!!!
#ifdef PARALLEL
  dserror("not implemented, yet");
#endif
  return label;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::PositionWithinAnyInfluencingCondition(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih,
    const std::set<int>&              xlabelset
)
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - PositionWithinAnyInfluencingCondition");
 
  const int label = PositionWithinCondition(x_in, ih);
  
  bool compute = false;
  if (label == 0) // fluid
  {
    compute = true;
  }
  else if (xlabelset.size() > 1)
  {
    compute = true;
  }
  else if (xlabelset.find(label) == xlabelset.end())
  {
    compute = true;
  }
  //compute = true;
  // TODO: in parallel, we have to ask all processors, whether there is any match!!!!
#ifdef PARALLEL
  dserror("not implemented, yet");
#endif
  return compute;
}



#endif  // #ifdef CCADISCRET
