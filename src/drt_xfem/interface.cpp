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
#include "intersection.H"
#include "../io/gmsh.H"

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
	  XFEM::Intersection is;
	  is.computeIntersection(xfemdis,cutterdis,elementalDomainIntCells_,elementalBoundaryIntCells_);
	  cout << "numcuttedelements = " << elementalDomainIntCells_.size() << endl;
	  // debug: write both meshes to file in Gmsh format
	  ofstream f_system("elements_coupled_system.pos");
	  f_system << IO::GMSH::disToString("Fluid", 0.0, xfemdis, elementalDomainIntCells_);
	  f_system << IO::GMSH::disToString("Solid", 1.0, cutterdis);
	  f_system << IO::GMSH::getConfigString(2);
	  f_system.close();
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
string XFEM::InterfaceHandle::toString() const
{
	stringstream s(" ");
	return s.str();
}

#endif  // #ifdef CCADISCRET
