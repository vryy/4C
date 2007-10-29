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

using namespace XFEM;

namespace GMSH
{

std::string elementToGmshString(DRT::Element* ele);

std::string intCellToGmshString(DRT::Element* ele, const IntCell ele);

std::string meshToGmshString(std::string s, RefCountPtr<DRT::Discretization> dis);

std::string meshToGmshString(std::string s, RefCountPtr<DRT::Discretization> dis,
        map<int, vector<IntCell> >&  elementIntCellMap);

};

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
