/*!
 * \file io_gmsh_xfem_extension.H
 * \brief simple element print library for Gmsh (mostly for debugging, only sequential)
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

#include "../drt_geometry/integrationcell.H"
#include "../drt_lib/drt_discret.H"

namespace IO
{
  namespace GMSH
  {
    //! print discretization including intersection information in initial configuration (t = 0)
    void XdisToStream(
        const std::string& text,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        const std::map<int,GEO::DomainIntCells >& elementDomainIntCellsMap,
        const std::map<int,GEO::BoundaryIntCells >& elementBoundaryIntCellsMap,
        std::ostream& s);

    //! print discretization including intersection information in initial configuration (t = 0)
    std::string XdisToString(
        const std::string& text,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        const std::map<int,GEO::DomainIntCells >& elementDomainIntCellsMap,
        const std::map<int,GEO::BoundaryIntCells >& elementBoundaryIntCellsMap);

  }
}

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
