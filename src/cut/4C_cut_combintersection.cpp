/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief provides the basic functionality for cutting a mesh with a level set function and/or a
       mesh


\level 2
*/
/*------------------------------------------------------------------------------------------------*/
#include "4C_cut_combintersection.hpp"

#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_levelsetside.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_tolerance.hpp"  // for EXTENDED_CUT_DEBUG_OUTPUT

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------*
 * constructur for Combined Intersection class (Levelset and Mesh intersection in one class)
 *-----------------------------------------------------------------------------------------*/

CORE::GEO::CUT::CombIntersection::CombIntersection(int myrank)
    : ParentIntersection(myrank), LevelSetIntersection(myrank, false), MeshIntersection(1, myrank)
{
  // call also the ParentIntersection-constructor first, otherwise according to the public virtual
  // inheritance of diamond shape the standard constructor of ParentIntersection() with default
  // value of myrank is called
}


void CORE::GEO::CUT::CombIntersection::Cut(bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::GEO::CUT --- 4/6 --- Cut_Intersection");

  if (myrank_ == 0 and screenoutput) IO::cout << "\t * 4/6 Cut_Intersection ...";

  Mesh& m = NormalMesh();

  // Remark: we assume that there is no overlap between levelset-isocontour and mesh

  // find cut points with levelset-side
  if (side_ != Teuchos::null)
  {
    m.Cut(*side_);
  }
//
// find cut points with cut mesh
#if EXTENDED_CUT_DEBUG_OUTPUT
  if (myrank_ == 0 and screenoutput) IO::cout << "\t\t...Finding cut points";
  double t_start = Teuchos::Time::wallTime();
#endif
  m.FindCutPoints();
#if EXTENDED_CUT_DEBUG_OUTPUT
  double t_diff = Teuchos::Time::wallTime() - t_start;
  if (myrank_ == 0 and screenoutput)
    IO::cout << "\t\t\t... Success (" << t_diff << " secs)" << IO::endl;
  if (myrank_ == 0 and screenoutput) IO::cout << "\t\t...Finding cut lines" << IO::endl;
  t_start = Teuchos::Time::wallTime();
#endif
  m.MakeCutLines();
#if EXTENDED_CUT_DEBUG_OUTPUT
  t_diff = Teuchos::Time::wallTime() - t_start;
  if (myrank_ == 0 and screenoutput)
    IO::cout << "\t\t\t... Success (" << t_diff << " secs)" << IO::endl;
  if (myrank_ == 0 and screenoutput) IO::cout << "\t\t...Finding cut facets" << IO::endl;
  t_start = Teuchos::Time::wallTime();
#endif
  m.MakeFacets();
#if EXTENDED_CUT_DEBUG_OUTPUT
  t_diff = Teuchos::Time::wallTime() - t_start;
  if (myrank_ == 0 and screenoutput)
    IO::cout << "\t\t\t... Success (" << t_diff << " secs)" << IO::endl;
  if (myrank_ == 0 and screenoutput) IO::cout << "\t\t...Finding cut cells" << IO::endl;
  t_start = Teuchos::Time::wallTime();
#endif
  m.MakeVolumeCells();
#if EXTENDED_CUT_DEBUG_OUTPUT
  t_diff = Teuchos::Time::wallTime() - t_start;
  if (myrank_ == 0 and screenoutput)
    IO::cout << "\t\t\t... Success (" << t_diff << " secs)" << IO::endl;
#endif
}


void CORE::GEO::CUT::CombIntersection::FindNodePositions()
{
  // TODO: this function and the overall inside-outside position strategy still has to be adapted
  // for more complex cases

  // NOTE: this will only work if mesh-cut area and level-set cut area are not overlapping

  Mesh& m = NormalMesh();

  // first, set the position for the mesh cut
  m.FindNodePositions();

  // second, set the position for the level-set cut
  m.FindLSNodePositions();
}


void CORE::GEO::CUT::CombIntersection::AddElement(int eid, const std::vector<int>& nids,
    const CORE::LINALG::SerialDenseMatrix& xyz, CORE::FE::CellType distype, const double* lsv,
    const bool lsv_only_plus_domain)
{
  CORE::GEO::CUT::ElementHandle* e = nullptr;

  // consider level-set values to decide whether the element has to be added or not
  if (lsv != nullptr)
  {
    // NOTE: dependent on whether one or two phases are used for the computation, the number of
    // degrees of freedom is determined via the cut status of elements, if both fluid phases have to
    // be considered, we have to add only cut elements, as uncut elements always carry physical
    // degrees of freedom if only the plus domain is a physical field, we have to add also elements
    // with pure negative level-set values (nodes in the ghost-domain) such that the
    // Dofset-Management does not produce degrees of freedom for such nodes

    e = LevelSetIntersection::AddElement(eid, nids, xyz, distype, lsv, lsv_only_plus_domain);
  }

  // no check necessary if element lies within bounding box of cut surface
  if (e != nullptr) return;

  MeshIntersection::AddElement(eid, nids, xyz, distype, lsv);
}

void CORE::GEO::CUT::CombIntersection::AddLevelSetSide(int levelset_side)
{
  LevelSetIntersection::AddCutSide(levelset_side);
}

void CORE::GEO::CUT::CombIntersection::AddMeshCuttingSide(int sid, const std::vector<int>& nids,
    const CORE::LINALG::SerialDenseMatrix& xyz, CORE::FE::CellType distype, int mi)
{
  MeshIntersection::AddCutSide(sid, nids, xyz, distype, mi);
}

FOUR_C_NAMESPACE_CLOSE
