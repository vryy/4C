/*!
 * \file io_gmsh_xfem_extension.H
 * \brief simple element print library for Gmsh (mostly for debuging, only sequential)
 * 
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 * 
 * For a complete description, see the section about post processing formations in the Gmsh documentation 
 */
#ifndef IO_GMSH_XFEM_EXTENSION_H
#define IO_GMSH_XFEM_EXTENSION_H
#ifdef CCADISCRET

#include <blitz/array.h>
#include <string>
#include <vector>

#include "../drt_geometry/integrationcell.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

namespace IO
{
  //! routines to transform results (elements, nodes, integrationcells, whole discretizations...) to Gmsh output
  namespace GMSH
  {
    //! print discretization including intersection information in initial configuration (t = 0)
    std::string XdisToString(
        const std::string& s,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        const std::map<int,XFEM::DomainIntCells >& elementDomainIntCellsMap,
        const std::map<int,XFEM::BoundaryIntCells >& elementBoundaryIntCellsMap);

  }
}

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
