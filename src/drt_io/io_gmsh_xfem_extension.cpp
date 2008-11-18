/*!
\file io_gmsh_xfem_extension.cpp

\brief simple element print library for Gmsh (debugging only)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "io_gmsh.H"
#include "io_gmsh_xfem_extension.H"

std::string IO::GMSH::XdisToString(
    const std::string& s,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis,
    const std::map<int, GEO::DomainIntCells >& elementDomainIntCellsMap,
    const std::map<int, GEO::BoundaryIntCells >& elementBoundaryIntCellsMap)
{
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" " << s << " Elements and Integration Cells \" {\n";
  
  for (int i=0; i<dis->NumMyRowElements(); ++i)
  {
    const DRT::Element* actele = dis->lRowElement(i);
    const int id = actele->Id();
    // print integration cells, if available
    if (elementDomainIntCellsMap.count(id) > 0)
    {
      const GEO::DomainIntCells cells = elementDomainIntCellsMap.find(id)->second;
      for (GEO::DomainIntCells::const_iterator cell = cells.begin(); cell
          != cells.end(); ++cell)
      {
        static LINALG::SerialDenseMatrix dxyz_ele(3, 27);
        cell->NodalPosXYZ(*actele, dxyz_ele);
        gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(),
            scalar, dxyz_ele) << "\n";
      }
      const GEO::BoundaryIntCells bcells = elementBoundaryIntCellsMap.find(id)->second;
      for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell
          != bcells.end(); ++bcell)
      {
        static LINALG::SerialDenseMatrix bxyz_ele(3, 9);
        bcell->NodalPosXYZ(*actele, bxyz_ele);
        gmshfilecontent << IO::GMSH::cellWithScalarToString(bcell->Shape(),
            scalar, bxyz_ele) << "\n";
      }
    }
    // print element
    gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(scalar, actele) << "\n";
  };
  gmshfilecontent << "};\n";
  return gmshfilecontent.str();
}


#endif // #ifdef CCADISCRET
