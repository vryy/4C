#ifndef CUT_SIDEHANDLE_H
#define CUT_SIDEHANDLE_H

#include <vector>

#include "cut_side.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Side;
  class Mesh;

class SideHandle
{
public:

  virtual ~SideHandle() {}
};

class LinearSideHandle : public SideHandle
{
public:

  LinearSideHandle()
    : side_( NULL )
  {
  }

  explicit LinearSideHandle( Side * s )
    : side_( s )
  {
  }

private:
  Side * side_;
};

class QuadraticSideHandle : public SideHandle
{
public:

protected:
  std::vector<Side*> subsides_;
  std::vector<Node*> nodes_;
};

class Tri6SideHandle : public QuadraticSideHandle
{
public:

  Tri6SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

};

class Quad8SideHandle : public QuadraticSideHandle
{
public:

  Quad8SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

};

class Quad9SideHandle : public QuadraticSideHandle
{
public:

  Quad9SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

};

}
}

#endif
