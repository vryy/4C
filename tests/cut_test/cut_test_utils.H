#ifndef GEOMETRY_TEST_UTILS_H
#define GEOMETRY_TEST_UTILS_H

#include <vector>
#include <Epetra_SerialDenseMatrix.h>

namespace GEO
{
  namespace CUT
  {
    class MeshIntersection;
    class Mesh;
    class Element;
    class Side;
  }
}

GEO::CUT::Element* create_hex8( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
// GEO::CUT::Element* create_hex20( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
// GEO::CUT::Element* create_hex27( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
GEO::CUT::Element* create_tet4( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
GEO::CUT::Element* create_wedge6( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
GEO::CUT::Element* create_pyramid5( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );

GEO::CUT::Side* create_quad4( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
// GEO::CUT::Side* create_quad8( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
// GEO::CUT::Side* create_quad9( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );

GEO::CUT::Element* create_hex8( GEO::CUT::Mesh & mesh, double dx=0, double dy=0, double dz=0 );

void create_hex8_mesh( GEO::CUT::Mesh & mesh, int rows, int cols, int depth );

void create_quad4_mesh( GEO::CUT::Mesh & mesh, int rows, int cols, std::vector<GEO::CUT::Side*> & sides );

void create_quad4_cylinder_mesh( GEO::CUT::MeshIntersection & intersection, double x, double y, int rows, int cols );

void cutmesh( GEO::CUT::Mesh & mesh );

void cutelement( GEO::CUT::Mesh & mesh, GEO::CUT::Element * element );

#endif
