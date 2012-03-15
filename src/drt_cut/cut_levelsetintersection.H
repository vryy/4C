#ifndef CUT_LEVELSETINTERSECTION_H
#define CUT_LEVELSETINTERSECTION_H

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

#include "cut_meshhandle.H"
#include "cut_pointpool.H"
#include "cut_levelsetside.H"
#include "cut_options.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;

  /*!
  \brief Interface class for the levelset cut.
  */
class LevelSetIntersection
{
public:

  LevelSetIntersection()
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( options_, 1, pp_, false ),
      side_( 1 )
  {
  }

  void AddElement( int eid,
                   const std::vector<int> & nids,
                   const Epetra_SerialDenseMatrix & xyz,
                   const double * lsv,
                   DRT::Element::DiscretizationType distype );

  /*!
  \brief Performs all the levelset cut operations.
   */
  void Cut( bool include_inner=true );

  ElementHandle * GetElement( int eid );

  Mesh & NormalMesh() { return mesh_.LinearMesh(); }

  Options & CreateOptions() { return options_; }

private:

  Teuchos::RCP<PointPool> pp_;

  Options options_;
  MeshHandle mesh_;
  LevelSetSide side_;
};

}
}

#endif
