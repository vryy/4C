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

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::InterfaceHandle::InterfaceHandle(
		const RefCountPtr<DRT::Discretization>        xfemdis, 
		const RefCountPtr<DRT::Discretization>        cutterdis) :
			xfemdis_(xfemdis),
			cutterdis_(cutterdis)
{
	  elementalDomainIntCells_.clear();
	  elementalBoundaryIntCells_.clear();
	  map< int, vector< DRT::Element* > >            cutterElementMap;
      map< int, RefCountPtr<DRT::Node> >             cutterNodeMap;
	  XFEM::Intersection is;
	  is.computeIntersection(
	          xfemdis,
	          cutterdis,
	          elementalDomainIntCells_,
	          elementalBoundaryIntCells_,
	          cutterElementMap,
	          cutterNodeMap);
	  std::cout << "numcuttedelements = " << elementalDomainIntCells_.size() << endl;
	  
	  
	  
	  
#if 0
	  // debug: write both meshes to file in Gmsh format
	  std::ofstream f_system("elements_coupled_system.pos");
	  f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis, elementalDomainIntCells_);
	  f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis);
	  {
	      stringstream gmshfilecontent;
	      gmshfilecontent << "View \" " << "CellCenter" << " Elements and Integration Cells \" {" << endl;
	      for (int i=0; i<xfemdis->NumMyColElements(); ++i)
	      {
	          DRT::Element* actele = xfemdis->lColElement(i);
	          XFEM::DomainIntCells elementDomainIntCells = this->GetDomainIntCells(actele->Id(), actele->Shape());
	          XFEM::DomainIntCells::const_iterator cell;
	          for(cell = elementDomainIntCells.begin(); cell != elementDomainIntCells.end(); ++cell )
	          {
	              const blitz::Array<double,1> cellcenterpos(cell->GetPhysicalCenterPosition(*actele));
	              gmshfilecontent << "SP(";
	              gmshfilecontent << scientific << cellcenterpos(0) << ",";
	              gmshfilecontent << scientific << cellcenterpos(1) << ",";
	              gmshfilecontent << scientific << cellcenterpos(2);
	              gmshfilecontent << "){";
	              gmshfilecontent << "0.0};" << endl;
	          };
	      };
	      gmshfilecontent << "};" << endl;
	      f_system << gmshfilecontent.str();
	  }
	  
	  f_system << IO::GMSH::getConfigString(3);
	  f_system.close();
#endif
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

XFEM::DomainIntCells XFEM::InterfaceHandle::GetDomainIntCells(
        const int gid,
        const DRT::Element::DiscretizationType distype
        )const
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
        const XFEM::BoundaryIntCells cells;
        return cells;
    }
    return tmp->second;
}

bool XFEM::InterfaceHandle::PositionWithinSpecificClosedRegion(
        const blitz::Array<double,1>& actpos,
        const int xfemcondition_label
        ) const
{
    return PositionWithinCondition(actpos, xfemcondition_label,cutterdis_);
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


#endif  // #ifdef CCADISCRET
