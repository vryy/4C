/*!
\file gmsh.H

\brief simple element print library for Gmsh (debuging only)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifndef GMSH_H_
#define GMSH_H_
#ifdef CCADISCRET

#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include <string>
#include <vector>

using namespace XFEM;

namespace GMSH
{

std::string elementToGmshString(DRT::Element* ele);

std::string intCellToGmshString(DRT::Element* ele, const Integrationcell ele);

};

#endif
#endif /*GMSH_H_*/
