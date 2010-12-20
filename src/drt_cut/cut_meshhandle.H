#ifndef CUT_MESHHANDLE_H
#define CUT_MESHHANDLE_H

#include <list>
#include <map>
#include <set>
#include <vector>

#include <Teuchos_RCP.hpp>

#include "cut_mesh.H"
#include "cut_sidehandle.H"
#include "cut_elementhandle.H"

namespace GEO
{
namespace CUT
{
  class Element;
  class LinearElement;
  class QuadraticElement;
  class Mesh;

  class SideHandle;
  class LinearSideHandle;
  class QuadraticSideHandle;

  class ElementHandle;
  class LinearElementHandle;
  class QuadraticElementHandle;

class MeshHandle
{
public:

  MeshHandle( double norm=1, Teuchos::RCP<PointPool> pp=Teuchos::null, bool cutmesh=false )
    : mesh_( norm, pp, cutmesh )
  {
  }

  SideHandle * CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype );

  ElementHandle * CreateElement( int eid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype );

  SideHandle * GetSide( int eid );

  ElementHandle * GetElement( int eid );

  Mesh & LinearMesh() { return mesh_; }

private:

  Mesh mesh_;
  std::map<int, LinearElementHandle> linearelements_;
  std::map<int, Teuchos::RCP<QuadraticElementHandle> > quadraticelements_;
  std::map<int, LinearSideHandle> linearsides_;
  std::map<int, Teuchos::RCP<QuadraticSideHandle> > quadraticsides_;
};

}
}

#endif
