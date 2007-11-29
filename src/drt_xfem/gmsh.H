/*!
\file gmsh.H

\brief simple element print library for Gmsh (mostly for debuging, only sequential)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifndef GMSH_H
#define GMSH_H
#ifdef CCADISCRET

#include <blitz/array.h>
#include <string>
#include <vector>

#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

namespace GMSH
{

std::string elementToString(const double scalar, DRT::Element* ele);
std::string elementToString(const vector<double> scalarfield, DRT::Element* ele);
std::string elementToString(const vector<vector<double> > vectorfield, DRT::Element* ele);

std::string cellToString(const double scalar, const XFEM::DomainIntCell& ele, DRT::Element* ele);

std::string cellToString(const double scalar, const XFEM::BoundaryIntCell& ele, DRT::Element* ele);

std::string disToString(const std::string s, const double scalar, RefCountPtr<DRT::Discretization> dis);

std::string disToString(const std::string s, const double scalar, RefCountPtr<DRT::Discretization> dis,
        map<int, XFEM::DomainIntCells >&  elementDomainIntCellsMap);

std::string getConfigString(const int numview);

std::string distypeToGmshElementHeader(const DRT::Element::DiscretizationType distype);

int distypeToGmshNumNode(const DRT::Element::DiscretizationType distype);

};

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
