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
#include "../drt_lib/drt_globalproblem.H"
#include "xfsi_searchtree.H"
#include "xfem_condition.H"
#include "../drt_io/io_gmsh.H"
#include "integrationcell.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::InterfaceHandle(
    const RCP<DRT::Discretization>        xfemdis, 
    const RCP<DRT::Discretization>        cutterdis,
    const Epetra_Vector&                  idispcol) :
      xfemdis_(xfemdis),
      cutterdis_(cutterdis)
{
  currentcutterpositions_.clear();
  {
    for (int lid = 0; lid < cutterdis->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = cutterdis->lColNode(lid);
      vector<int> lm;
      lm.reserve(3);
      cutterdis->Dof(node, lm);
      vector<double> mydisp(3);
      DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);
      static BlitzVec3 currpos;
      currpos(0) = node->X()[0] + mydisp[0];
      currpos(1) = node->X()[1] + mydisp[1];
      currpos(2) = node->X()[2] + mydisp[2];
      currentcutterpositions_[node->Id()] = currpos;
    }
  }
  
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
  
  // sanity check, whether, we realy have integration cells in the map
  for (std::map<int,DomainIntCells>::const_iterator 
      tmp = elementalDomainIntCells_.begin();
      tmp != elementalDomainIntCells_.end();
      ++tmp)
  {
    dsassert(tmp->second.empty() == false, "this is a bug!");
  }
  
  // sanity check, whether, we realy have integration cells in the map
  for (std::map<int,BoundaryIntCells>::const_iterator 
      tmp = elementalBoundaryIntCells_.begin();
      tmp != elementalBoundaryIntCells_.end();
      ++tmp)
  {
    dsassert(tmp->second.empty() == false, "this is a bug!");
  }
  
  elementsByLabel_.clear();
  CollectElementsByXFEMCouplingLabel(*cutterdis, elementsByLabel_);

  //cout << "create new xTree_ object" << endl;
  const BlitzMat3x2 cutterAABB = XFEM::getXAABBofDis(*cutterdis,currentcutterpositions_);
  const BlitzMat3x2 xfemAABB = XFEM::getXAABBofDis(*xfemdis);
  const BlitzMat3x2 AABB = XFEM::mergeAABB(cutterAABB, xfemAABB);
  xTree_ = rcp(new XSearchTree(AABB));
}
		
/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::~InterfaceHandle()
{
    return;
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
std::string XFEM::InterfaceHandle::toString() const
{
	std::stringstream s(" ");
	return s.str();
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::InterfaceHandle::toGmsh(const int step) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (xfemparams.get<std::string>("GMSH_DEBUG_OUT") == "Yes");
  
  const bool gmsh_tree_output = false;
  
  if (gmshdebugout)
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    filename << allfiles.outputfile_kenner << "_elements_coupled_system_" << std::setw(5) << setfill('0') << step << ".pos";
    std::cout << "writing '"<<filename.str()<<"'...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_);
    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, currentcutterpositions_);
    //f_system << IO::GMSH::getConfigString(3);
    f_system.close();
    cout << " done" << endl;
  }
  
  if (gmshdebugout)
  {
    std::stringstream filename;
    filename << allfiles.outputfile_kenner << "_domains_" << std::setw(5) << setfill('0') << step << ".pos";
    cout << "writing '"<<filename.str()<<"...";

    std::ofstream f_system(filename.str().c_str());
//    f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_);
//    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, currentcutterpositions_);
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
          int closestElementId;
          double distance;
          const int domain_id = PositionWithinCondition(cellcenterpos, *this, closestElementId, distance);
          stringstream text;
          text.precision(3);
          text << "(<"<< (actele->Id()) << "," << closestElementId << ">,\n"<< fixed << distance << ")";
          gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), domain_id*100000+(closestElementId), cellpos) << endl;
          BlitzMat point(3,1);
          point(0,0)=cellcenterpos(0);
          point(1,0)=cellcenterpos(1);
          point(2,0)=cellcenterpos(2);
             
        };
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    //f_system << IO::GMSH::getConfigString(3);
    f_system.close();
    cout << " done" << endl;
  }
  if (gmsh_tree_output)
  {
    // debug: write information about which structure we are in
    std::stringstream filenameP;
    filenameP << allfiles.outputfile_kenner << "_points_" << std::setw(5) << setfill('0') << step << ".pos";
    std::stringstream filenameAnnotations;
    filenameAnnotations << allfiles.outputfile_kenner << "_annotations_" << std::setw(5) << setfill('0') << step << ".pos";

    cout <<endl << "writing '"<<filenameP.str()<<"...";
    std::ofstream f_systemP(filenameP.str().c_str());
    cout <<endl << "writing '"<<filenameAnnotations.str()<<"...";
    std::ofstream f_systemAnnotations(filenameAnnotations.str().c_str());
    {
      // stringstream for cellcenter points
      stringstream gmshfilecontentP;
      gmshfilecontentP << "View \" " << "CellCenter of Elements and Integration Cells \" {" << endl;
      stringstream gmshfilecontentAnnotations;
      gmshfilecontentAnnotations << "View \" " << "Annotations for cell centers \" {" << endl;
      
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
          
          int closestElementId;
          double distance;
 
          const int domain_id = PositionWithinCondition(cellcenterpos, *this, closestElementId, distance);
          
          BlitzMat point(3,1);
          point(0,0)=cellcenterpos(0);
          point(1,0)=cellcenterpos(1);
          point(2,0)=cellcenterpos(2);

          stringstream text;
          text.precision(3);
          text << "(<"<< (actele->Id()) << "," << closestElementId << ">,\n"<< fixed << distance << ")";
       
          gmshfilecontentAnnotations << IO::GMSH::text3dToString(cellcenterpos, text.str(), 2) << endl;              
          gmshfilecontentP << IO::GMSH::cellWithScalarToString(DRT::Element::point1, (actele->Id()), point) << endl;              
        };
      };
      gmshfilecontentP << "};" << endl;
      gmshfilecontentAnnotations << "};" << endl;
      f_systemP << gmshfilecontentP.str();
      f_systemAnnotations << gmshfilecontentAnnotations.str();
    }
    //f_system << IO::GMSH::getConfigString(3);
    f_systemP.close();
    f_systemAnnotations.close();
    cout << " done" << endl;
    xTree_->printTree(allfiles.outputfile_kenner, step);
  }
  return;
}

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
 |  CLI:    checks if a position is within condition-enclosed region      p.ede 05/08|   
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinCondition(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih
)
{
  int closestElementId;
  double distance; 
  return PositionWithinCondition(x_in,ih,closestElementId, distance);
}

/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within condition-enclosed region      p.ede 05/08|   
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinCondition(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih,
    int& closestElementId, 
    double& distance
)
{
  //PositionWithinConditionBruteForce(x_in, ih, posInCondition);
  const int label = PositionWithinConditionTree(x_in, ih, closestElementId,distance);
//  const std::map<int,set<int> >& elementsByLabel = *(ih.elementsByLabel());
//  for(std::map<int,set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
//   {
//     const int label = conditer->first;
//     if (posInCondition1[label]!=posInCondition2[label])
//     {
////       cout << " bruteforce posInCondition[" << label <<"] = "<< posInCondition1[label] << endl;    
////       cout << "xsearchtree posInCondition[" << label <<"] = "<< posInCondition2[label] << endl;
//       //cout <<  posInCondition1[label] << " " << posInCondition2[label] << endl;
//       flush(cout);  
//       dserror("results for searchtree and brute force do not match");
//     }
//   }
  //posInCondition = posInCondition1;
    
  // TODO: in parallel, we have to ask all processors, whether there is any match!!!!
#ifdef PARALLEL
  dserror("not implemented, yet");
#endif
  return label;
}

/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within condition-enclosed region      a.ger 12/07|   
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinConditionBruteForce(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih
)
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - PositionWithinCondition");
  //init
  std::map<int,bool>  posInCondition; // not really needed, but this method could go away soon, so it won't be cleaned
  const std::map<int,set<int> >& elementsByLabel = *(ih.elementsByLabel());
  
  /////////////////
  // loop labels
  /////////////////
  int label = 0;
  for(std::map<int,set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
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
      const BlitzMat xyze_cutter(getCurrentNodalPositions(cutterele, *ih.currentcutterpositions()));
      double distance = 0.0;
      static BlitzVec2 eleCoord;
      static BlitzVec3 normal;
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
 |  CLI:    checks if a position is within condition-enclosed region      p.ede 05/08|   
 *----------------------------------------------------------------------*/
int XFEM::PositionWithinConditionTree(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih,
    int& closestElementId, 
    double& distance
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - search - PositionWithinConditionTree");
  Teuchos::RCP<XSearchTree> xt = ih.getSearchTree(); // pointer is constant, object is not!
  int l = xt->queryPointType(*ih.cutterdis() , *ih.currentcutterpositions(), x_in, closestElementId, distance);

  #ifdef PARALLEL
  dserror("not implemented, yet");
#endif

  return l;
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within condition-enclosed region      a.ger 12/07|   
 *----------------------------------------------------------------------*/
bool XFEM::PositionWithinAnyInfluencingCondition(
    const BlitzVec3&                  x_in,
    const XFEM::InterfaceHandle&      ih,
    const std::set<int>&              xlabelset
)
{
  
  TEUCHOS_FUNC_TIME_MONITOR(" - search - PositionWithinAnyInfluencingCondition");
  
  if (xlabelset.size() > 1)
  {
    dserror("this function has to be extended to allow contact problems!");
  }
  int closestElementId; 
  double distance;
  const int label = PositionWithinCondition(x_in, ih,closestElementId, distance);
  bool compute = false;
  if (xlabelset.find(label) == xlabelset.end())
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


//const XFEM::InterfaceHandle::emptyBoundaryIntCells_ = XFEM::BoundaryIntCells(0);

#endif  // #ifdef CCADISCRET
