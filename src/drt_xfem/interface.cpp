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
#include "xfem_condition.H"

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::InterfaceHandle(
		const RCP<DRT::Discretization>        xfemdis, 
		const RCP<DRT::Discretization>        cutterdis,
		const RCP<Epetra_Vector>              idispcol) :
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
      DRT::UTILS::ExtractMyValues(*idispcol,mydisp,lm);
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
	is.computeIntersection(
	        xfemdis,
	        cutterdis,
	        currentcutterpositions_,
	        elementalDomainIntCells_,
	        elementalBoundaryIntCells_);
	
	std::cout << "numcuttedelements (elementalDomainIntCells_)   = " << elementalDomainIntCells_.size() << endl;
	std::cout << "numcuttedelements (elementalBoundaryIntCells_) = " << elementalBoundaryIntCells_.size() << endl;
	dsassert(elementalDomainIntCells_.size() == elementalBoundaryIntCells_.size(), "mismatch in cutted elements maps!");
	  
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
  {
    // debug: write both meshes to file in Gmsh format
    std::stringstream filename;
    filename << "elements_coupled_system_" << std::setw(5) << setfill('0') << step << ".pos";
    std::cout << "writing '"<<filename.str()<<"'...";
    std::ofstream f_system(filename.str().c_str());
    f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis_, elementalDomainIntCells_, elementalBoundaryIntCells_);
    f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, currentcutterpositions_);
//    {
//        stringstream gmshfilecontent;
//        gmshfilecontent << "View \" " << "CellCenter" << " Elements and Integration Cells \" {" << endl;
//        for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
//        {
//            const DRT::Element* xfemele = xfemdis_->lColElement(i);
//            const XFEM::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(xfemele->Id(), xfemele->Shape());
//            XFEM::DomainIntCells::const_iterator cell;
//            for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
//            {
//                const BlitzVec3 cellcenterpos(cell->GetPhysicalCenterPosition(*xfemele));
//                gmshfilecontent << "SP(";
//                gmshfilecontent << scientific << cellcenterpos(0) << ",";
//                gmshfilecontent << scientific << cellcenterpos(1) << ",";
//                gmshfilecontent << scientific << cellcenterpos(2);
//                gmshfilecontent << "){";
//                gmshfilecontent << "0.0};" << endl;
//            };
//        };
//        gmshfilecontent << "};" << endl;
//        f_system << gmshfilecontent.str();
//    }
    f_system << IO::GMSH::getConfigString(3);
    f_system.close();
    std::cout << " done" << endl;
  }
  
  // debug: write information about which structure we are in
  {
    std::stringstream filename;
    filename << "domains_" << std::setw(5) << setfill('0') << step << ".pos";
    std::cout << "writing '"<<filename.str()<<"'...";
    std::ofstream f_system(filename.str().c_str());
    //f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis, elementalDomainIntCells_, elementalBoundaryIntCells_);
    //f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis_, currentcutterpositions_);
    {
        BlitzMat cellpos(3,27);
      
        map<int,bool> posInCondition;
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "Domains using CellCenter of Elements and Integration Cells \" {" << endl;
        for (int i=0; i<xfemdis_->NumMyColElements(); ++i)
        {
            DRT::Element* actele = xfemdis_->lColElement(i);
            const XFEM::DomainIntCells& elementDomainIntCells = this->GetDomainIntCells(actele->Id(), actele->Shape());
            XFEM::DomainIntCells::const_iterator cell;
            for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
            {
              cell->NodalPosXYZ(*actele, cellpos);
              const BlitzVec3 cellcenterpos(cell->GetPhysicalCenterPosition(*actele));
              //cout << cellcenterpos << endl;
              PositionWithinCondition(cellcenterpos, *this, posInCondition);
              int domain_id = 0;
              // loop conditions
              for (map<int,bool>::const_iterator entry = posInCondition.begin(); entry != posInCondition.end(); ++entry)
              {
                const int actdomain_id = entry->first;
                const bool inside_condition = entry->second;
                //cout << "Domain: " << actdomain_id << " -> inside? " << inside_condition << endl;
                if (inside_condition)
                {
                  //cout << "inside" << endl;
                  domain_id = actdomain_id;
                  break;
                }
              }
              gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(), domain_id, cellpos) << endl;
            };
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    //f_system << IO::GMSH::getConfigString(3);
    f_system.close();
    std::cout << " done" << endl;
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
    std::map<int,DomainIntCells>::const_iterator tmp = elementalDomainIntCells_.find(element_gid);
    if (tmp == elementalDomainIntCells_.end())
    {   
        return false;
    }
    else
    {
        return true;
    }
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within condition-enclosed region      a.ger 12/07|   
 *----------------------------------------------------------------------*/
void XFEM::PositionWithinCondition(
        const BlitzVec3&                  x_in,
        const XFEM::InterfaceHandle&      ih,
        std::map<int,bool>&               posInCondition
    )
{
  posInCondition.clear();

  const std::map<int,set<int> >& elementsByLabel = *(ih.elementsByLabel());
  
  /////////////////
  // loop labels
  /////////////////
  for(std::map<int,set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    const int label = conditer->first;
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
      }
    }
    
  } // end loop label

    // TODO: in parallel, we have to ask all processors, whether there is any match!!!!
#ifdef PARALLEL
    dserror("not implemented, yet");
#endif
    //exit(0);
    return;
}

//const XFEM::InterfaceHandle::emptyBoundaryIntCells_ = XFEM::BoundaryIntCells(0);

#endif  // #ifdef CCADISCRET
