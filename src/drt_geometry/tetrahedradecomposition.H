/*!-----------------------------------------------------------------------------------------------*
\file tetrahedradecomposition.H

  \brief tetrahedralization procedure based on level set function

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

<pre>
\level 3
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TETRAHEDRADECOMPOSITION_H_
#define TETRAHEDRADECOMPOSITION_H_


#include "../drt_lib/drt_element.H"
#include "geo_utils.H"


namespace COMBUST
{
  class RefinementCell;
}

namespace GEO
{

class TetrahedraDecomposition
{

public:

  TetrahedraDecomposition(
      const COMBUST::RefinementCell* cell,
      GEO::BoundaryIntCells& listBoundaryIntCellsperEle,
      GEO::DomainIntCells & listDomainIntCellsperEle
  );

private:

  void DecomposeIntoTetrahedra(
      std::vector<GEO::DomainIntCell>& tets,
      LINALG::SerialDenseMatrix& tets_gvalues
  );

  void IntersectTetrahdera(
      std::vector<GEO::DomainIntCell>& tets,
      GEO::BoundaryIntCells& tris,
      LINALG::SerialDenseMatrix& tets_gvalues
  );

  int  FindIntersectionPointsTet4(
      GEO::DomainIntCell& tet,
      std::vector<double>& tet_gvalues,
      LINALG::SerialDenseMatrix& vertices,
      LINALG::SerialDenseMatrix& nodes,
      int& sign_node_four
  );

  void checkTriangleOrientation(
      LINALG::SerialDenseMatrix& coords,
      LINALG::SerialDenseMatrix& pcoords,
      LINALG::Matrix<3,1>& node_four,
      int sign_node_four
  );

  void createTetrahedra(
      unsigned int num_tetrahedra,
      LINALG::SerialDenseMatrix & tet_coords,
      LINALG::SerialDenseMatrix & vertices, const int* tetindex,
      LINALG::SerialDenseMatrix & tet_pcoords, LINALG::SerialDenseMatrix & nodes,
      bool & intersected, std::vector<GEO::DomainIntCell> & tets,
      int & sign_node_four,
      bool tetgplus[2]
  );

  void createTriangles(
      unsigned int num_triangles,
      LINALG::SerialDenseMatrix trianglecoord,
      LINALG::SerialDenseMatrix vertices,
      const int* triindex,
      LINALG::SerialDenseMatrix phystrianglecoord,
      LINALG::SerialDenseMatrix nodes,
      LINALG::Matrix<3,1> temp,
      int comparison_vertex,
      int sign_node_four,
      GEO::BoundaryIntCells & tris,
      bool intersected
  );

  //! compute average GfuncValue for integration cell
  bool GetIntCellDomainInElement(
      const LINALG::SerialDenseMatrix&       IntCellCoord,
      const std::vector<double>&             gfuncvalues_ele,
      const DRT::Element::DiscretizationType xfem_distype,
      const DRT::Element::DiscretizationType cell_distype
  );

  // refinement cell to be decomposed
  const COMBUST::RefinementCell* cell_;

  // number of degenerated tetrahedra
  int num_degeneratedTets_;

};
}

#endif /* TetrahedraDecomposition_H_ */
