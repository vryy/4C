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

void IO::GMSH::XdisToStream(
    const std::string& text,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis,
    const std::map<int, GEO::DomainIntCells >& elementDomainIntCellsMap,
    const std::map<int, GEO::BoundaryIntCells >& elementBoundaryIntCellsMap,
    std::ostream& s)
{
  s << "View \" " << text << " Elements and Integration Cells \" {\n";

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
        const LINALG::SerialDenseMatrix& dxyz_ele = cell->CellNodalPosXYZ();
        IO::GMSH::cellWithScalarToStream(cell->Shape(), scalar, dxyz_ele, s);
      }
    }
    if (elementBoundaryIntCellsMap.count(id) > 0)
    {
      const GEO::BoundaryIntCells bcells = elementBoundaryIntCellsMap.find(id)->second;
      for (GEO::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell
          != bcells.end(); ++bcell)
      {
        const LINALG::SerialDenseMatrix& bxyz_ele = bcell->CellNodalPosXYZ();
        IO::GMSH::cellWithScalarToStream(bcell->Shape(), scalar, bxyz_ele, s);
      }
    }
    // print element
    IO::GMSH::elementAtInitialPositionToStream(scalar, actele, s);
  };
  s << "};\n";
}

std::string IO::GMSH::XdisToString(
    const std::string& text,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis,
    const std::map<int, GEO::DomainIntCells >& elementDomainIntCellsMap,
    const std::map<int, GEO::BoundaryIntCells >& elementBoundaryIntCellsMap)
{
  std::ostringstream s;
  XdisToStream(text, scalar, dis, elementDomainIntCellsMap, elementBoundaryIntCellsMap, s);
  return s.str();
}

#endif // #ifdef CCADISCRET
