/*---------------------------------------------------------------------*/
/*!

\brief Tessellate the element by means of a Delauney triangulation

\level 3

\maintainer Christoph Ager
*/
/*---------------------------------------------------------------------*/

#include <stack>

#include "cut_tetmesh.H"
#include "cut_tetmeshintersection.H"
#include "cut_levelsetside.H"
#include "cut_volumecell.H"

extern "C"
{
#include <qhull_a.h>
}

namespace GEO
{
  namespace CUT
  {
    namespace
    {
      class NullFile
      {
       public:
        NullFile() : f_(NULL) {}

        ~NullFile()
        {
          if (f_ != NULL) fclose(f_);
        }

        operator FILE*()
        {
          if (f_ == NULL) f_ = fopen("/dev/null", "w");
          return f_;
        }

       private:
        FILE* f_;
      };
    }  // namespace
  }    // namespace CUT
}  // namespace GEO


GEO::CUT::TetMesh::TetMesh(
    const std::vector<Point*>& points, const plain_facet_set& facets, bool project)
    : points_(points), facets_(facets)
{
  std::vector<std::vector<int>> original_tets;

  CallQHull(points_, original_tets, project);

  tets_.reserve(original_tets.size());
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
  int counter = 0;
  int counter_saved = 0;
#endif
  for (std::vector<std::vector<int>>::iterator i = original_tets.begin(); i != original_tets.end();
       ++i)
  {
    std::vector<int>& t = *i;
    std::vector<Point*> tet;
    tet.reserve(t.size());
    for (std::vector<int>::iterator i = t.begin(); i != t.end(); ++i)
    {
      tet.push_back(points_[*i]);
    }
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << "counter: " << counter << ", ";
#endif
    if (IsValidTet(tet))
    {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      std::cout << " counter_saved: " << counter_saved;
      std::cout << ", IsValidTet( tet ): "
                << "TRUE";
      counter_saved++;
#endif
      // tets_.push_back( std::vector<int>() );
      // std::swap( tets_.back(), t );
      tets_.push_back(t);
    }
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    else
    {
      std::cout << "IsValidTet, threw away tet with nodes: " << t[0] << ", " << t[1] << ", " << t[2]
                << ", " << t[3] << std::endl;
      std::cout << "With volume: " << CalcVolumeOfTet(tet) << std::endl;
    }
#endif
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << std::endl;
    counter++;
#endif
  }

#if 0
  // if this is not the first cut, it might be fine not to have all points
  TestUsedPoints( tets_ );
#endif

  Init();

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
  GmshWriteCells();
#endif
}

/* Loop over the facets of the tetmesh. Can every facet (triangulated at this stage) be associated
   with an unique TET from the delauney-cell decomposition?

   @Magnus:
   WARNING: For a LevelSet-cut with many cut-points this might lead to a failure to fill the
   FacetMesh as no info about the triangulated point is provided. This is however taken care of in
   the TetMeshIntersection, (the TETs are not cut) It could however be dealt with if the
   triangulated point is added to the cut information. Might have to look into this....

*/
bool GEO::CUT::TetMesh::FillFacetMesh()
{
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    if (not FillFacet(f))  // Possible to FillFacet?
    {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      std::cout << "Facet unable to be filled! ";
      std::cout << "f->OnCutSide(): " << f->OnCutSide() << std::endl;
      f->Print();
      std::cout << std::endl;
#endif
      return false;
    }
    if (f->HasHoles())  // Does Facet contain holes?
    {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      std::cout << "f->HasHoles()";
#endif
      const plain_facet_set& holes = f->Holes();
      for (plain_facet_set::const_iterator i = holes.begin(); i != holes.end(); ++i)
      {
        Facet* h = *i;
        if (not FillFacet(h))
        {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
          std::cout << "Facet unable to be filled!";
#endif
          return false;
        }
      }
    }
  }
  return true;
}


/* From the input create TETs.
   1) First remove the TETs which are too small for arithmetic precision from the "accepted_tet_set"
   of TETs. However they might be needed later (for boundary-cell creation) so don't throw them
   away. 2) Try to "FillFacetMesh(...)", this means that for the given given element (either a tet
   or a hex8) find the facets on the boundaries and find if they are associated with an unique TET.

    2.1) If NOT assume the tesselation (i.e. on TET or more) is cut, and call the Tesselation
   recursively. 2.2) If YES fill the inner part of the element with accepted tets and finish the
   tesselation for this element.

*/
void GEO::CUT::TetMesh::CreateElementTets(Mesh& mesh, Element* element,
    const plain_volumecell_set& cells,
    const plain_side_set& cut_sides,  //<- cut_facets_ of parent ele.
    int count, bool tetcellsonly)
{
  FixBrokenTets();

#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
  int sum = 0;
  for (unsigned k = 0; k < accept_tets_.size(); k++)
  {
    sum += accept_tets_[k];
  }
  std::cout << "accept_tets_.sum(): " << sum << std::endl;
#endif

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
  std::cout << "Entered GEO::CUT::TetMesh::CreateElementTets(...), count: " << count << std::endl;
#endif
  if (count == 1)
  {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << "  OUTPUT ActiveCells and SurfaceCells..." << std::endl;
#endif
    GmshWriteActiveCells();
    GmshWriteSurfaceCells();
  }
#endif

  // Can the mesh be "filled from the facets"
  if (FillFacetMesh())
  {
    for (plain_volumecell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
    {
      VolumeCell* vc = *i;

      // Describes the "domain" of the volume cell (i.e. its tets and borders)
      Domain<4> cell_domain;
      // The tets which are to be added in this volume cell
      PlainEntitySet<4>& cell_members = cell_domain.Members();
      // The border (i.e. surface-tris) between the tets in the cell.
      PlainEntitySet<3>& cell_border = cell_domain.Border();

#ifdef DEBUGCUTLIBRARY
      std::vector<Side*> facet_sides;
#endif

      const plain_facet_set& facets = vc->Facets();
      for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
      {
        Facet* f = *i;
        FacetMesh& fm = facet_mesh_[f];
        const PlainEntitySet<3>& tris = fm.SurfaceTris();

#ifdef DEBUGCUTLIBRARY
        facet_sides.push_back(f->ParentSide());
#endif

        // All tris are added to the cell_border, i.e. (border_ in class Domain).
        //  These are later used to create "done_border_" later in SeedDomain.
        std::copy(tris.begin(), tris.end(), std::inserter(cell_border, cell_border.begin()));
      }

      for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
      {
        Facet* f = *i;
        SeedDomain(cell_domain, f);
      }

      if (cell_domain.Empty())
      {
        // Emergency. If this is the only volume cell within the element (an
        // element is surrounded by cut surfaces), try and force any available
        // tet into the volume cell. (This is a special case that turns up in
        // debug situations only.)

        if (cells.size() == 1)
        {
          for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
          {
            Facet* f = *i;
            SeedDomain(cell_domain, f, true);
          }
        }

        if (cell_domain.Empty())
        {
          // Assume the volume cell in question is degenerated and does not
          // contain any tets.

          // Test that the VC is completely on a cut-side. (added by Magnus)
#ifdef TEST_EMPTY_VC_IS_COMPLETELY_ON_CUTSIDE
#ifdef DEBUGCUTLIBRARY
          std::cout << "Warning a volume cell has no integration-cells!" << std::endl;
#endif
#ifdef DEBUGCUTLIBRARY
          double num_facet_counter = 0;
#endif
          bool volume_cell_on_facet = true;
          bool does_not_share_same_cutside = true;
          for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
          {
            Facet* f = *i;

            std::cout << "f->ShareSameCutSide(*facets.begin()): "
                      << f->ShareSameCutSide(*facets.begin()) << std::endl;

#ifdef DEBUGCUTLIBRARY
            std::cout << "num_facet_counter: " << num_facet_counter << std::endl;
            std::cout << "f->OnCutSide(): " << f->OnCutSide() << std::endl;
            std::cout << "f->ParentSide()->Id(): " << f->ParentSide()->Id() << std::endl;
            std::cout << "f->IsPlanar(mesh,f->CornerPoints()): "
                      << f->IsPlanar(mesh, f->CornerPoints()) << std::endl;
            f->PrintPointIds();
            num_facet_counter++;
#endif

            if (not f->ShareSameCutSide(*facets.begin()))
            {
              std::cout << "NOT SAME ParentSide!" << std::endl;
              does_not_share_same_cutside = false;
              continue;
            }

            if (not f->IsPlanar(mesh, f->CornerPoints()))
            {
              std::cout << "Warning: This volume cell has no integration cells AND not all on same "
                           "cut-side."
                        << std::endl;
              volume_cell_on_facet = false;
              continue;
            }
          }

          if (not volume_cell_on_facet and not does_not_share_same_cutside)
            throw std::runtime_error(
                "cell_domain.Empty() and facets not on same cut-side and not planar");
          else if (not volume_cell_on_facet)
            throw std::runtime_error("cell_domain.Empty() and facets not planar");
          else if (not does_not_share_same_cutside)
            throw std::runtime_error("cell_domain.Empty() and facets not on same cut-side");
#endif

          continue;
        }
      }

      cell_domain.Fill();

      //#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      //      std::cout << "cell_domain.Members().size(): " << cell_domain.Members().size();
      //      std::cout << ", cell_domain.Border().size(): " << cell_domain.Border().size() <<
      //      std::endl; std::cout << "Tet-IDs: "; for(unsigned j=0; j<cell_domain.Members().size();
      //      j++)
      //      {
      //        std::cout << cell_domain.Members()[j]->Id() << ", ";
      //      }
      //      std::cout << "Border-IDs: ";
      //      for(unsigned j=0; j<cell_domain.Border().size(); j++)
      //      {
      //        std::cout << std::endl << cell_domain.Border()[j]->Id() << ": ";
      //        for(unsigned l=0; l<3; l++)
      //        {
      //          std::cout << cell_domain.Border()[j]->GetHandle()[l] << ", ";
      //        }
      //      }
      //      std::cout << std::endl;
      //#endif


#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      GmshWriteTriSet("cell_border", cell_border);
      GmshWriteTetSet("cell_members", cell_members);
#endif

      std::vector<std::vector<Point*>> tets;
      tets.reserve(cell_members.size());

      for (PlainEntitySet<4>::iterator i = cell_members.begin(); i != cell_members.end(); ++i)
      {
        Entity<4>& t = **i;
        if (accept_tets_[t.Id()])
        {
          std::vector<int>& fixedtet = tets_[t.Id()];
          if (fixedtet.size() != 4) throw std::runtime_error("confused");
          tets.push_back(std::vector<Point*>(4));
          std::vector<Point*>& tet = tets.back();
          for (int i = 0; i < 4; ++i)
          {
            tet[i] = points_[fixedtet[i]];
          }
        }
      }

      // all facets (whose facet-mesh are just tri3s here) which are on cut surface obtain the
      // coordinates to create boundary integration cells then
      std::map<Facet*, std::vector<Point*>> sides_xyz;

      for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
      {
        Facet* f = *i;
        if (f->OnBoundaryCellSide())  // This is entered for both sides of the volume cell.
        {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
          std::cout << "facet on cut side: " << std::endl;
          f->PrintPointIds();
#endif

          FacetMesh& fm = facet_mesh_[f];
          const PlainEntitySet<3>& tris = fm.SurfaceTris();  // tris from the triangulated surface
          std::vector<Point*>& side_coords =
              sides_xyz[f];  // create entry for facet and get a reference to the facets side
                             // coordinates
          std::vector<std::vector<int>> sides;
          FindProperSides(tris, sides, &cell_members);
          CollectCoordinates(sides, side_coords);  // fill the side coordinates, if all the side's
                                                   // coordinates are on cut surface
        }
      }

      vc->CreateTet4IntegrationCells(mesh, tets, sides_xyz);

      if (vc->Empty())
      {
        throw std::runtime_error("empty volume cell detected");
      }
    }
  }
  else
  {
    // if ( count <= 3 )
    {
      TetMeshIntersection intersection(
          mesh.CreateOptions(), element, tets_, accept_tets_, points_, cut_sides);
      intersection.Cut(
          mesh, element, cells, count, tetcellsonly);  // writes tetmesh output in debug mode!!!
    }
  }
}

/* Initialize the data-structures needed for the TetMesh. The call to Qhull has to be performed
 * before.
 */
void GEO::CUT::TetMesh::Init()
{
  unsigned numtets = tets_.size();

  for (unsigned i = 0; i < numtets; ++i)
  {
    const std::vector<int>& t = tets_[i];
    tet_entities_.push_back(Entity<4>(i, Handle<4>(&t[0])));
  }

  for (std::vector<Entity<4>>::iterator i = tet_entities_.begin(); i != tet_entities_.end(); ++i)
  {
    Entity<4>& e = *i;
    e.CreateChildren(tet_surfaces_);
  }

  for (std::map<Handle<3>, Entity<3>>::iterator i = tet_surfaces_.begin(); i != tet_surfaces_.end();
       ++i)
  {
    Entity<3>& e = i->second;
    e.CreateChildren(tet_lines_);
  }

  accept_tets_.resize(tets_.size());
  std::fill(accept_tets_.begin(), accept_tets_.end(), true);

  //   active_surface_tris_.resize( tet_surfaces_.size() );
  //   std::fill( active_surface_tris_.begin(), active_surface_tris_.end(), false );
}


/* Call up QHull and let it create the delauney-cells. There are different options for how QHull can
   be called. 3 options are tested in order, if the previous option failed. Fill the std::vector
   tets with tetrahedrons provided from the delauney triangulization of the "convex"-hull.

   One might want to take a look at the call to Qhull to improve it for non-local coordinate cuts.
*/
void GEO::CUT::TetMesh::CallQHull(
    const std::vector<Point*>& points, std::vector<std::vector<int>>& tets, bool project)
{
  const int dim = 3;
  const int n = points.size();

  if (n < 4)
  {
    dserror("where coming from?");
    run_time_error("illegal element topology");
  }
  //   if ( n == 4 )
  //   {
  //     throw std::runtime_error( "no need to triangulate" );
  //   }

  std::vector<double> coordinates(dim * n);

  if (project)  // Never used... Not working properly either, I think...
  {
    LINALG::Matrix<3, 1> m;
    m = 0;
    double scale = 1. / n;
    for (int i = 0; i < n; ++i)  // Find mid-point (m)
    {
      Point* p = points[i];
      LINALG::Matrix<3, 1> x(p->X());
      m.Update(scale, x, 1);
    }
    double length = 0;
    LINALG::Matrix<3, 1> l;
    for (int i = 0; i < n;
         ++i)  // Find the distance to the point furthest away from the mid-point (length)
    {
      Point* p = points[i];
      LINALG::Matrix<3, 1> x(p->X());
      l = m;
      l.Update(1, x, -1);
      double n = l.Norm2();
      length = std::max(n, length);
    }
#ifdef DEBUGCUTLIBRARY
    std::ofstream pointfile("points.plot");
    pointfile << m(0) << " " << m(1) << " " << m(2) << "\n";
#endif
    for (int i = 0; i < n; ++i)
    {
      Point* p = points[i];
      LINALG::Matrix<3, 1> x(p->X());
      l = m;
      l.Update(1, x, -1);
      double n = l.Norm2();
      l.Scale(length / n);
      l.Update(1, m, 1);
      std::copy(l.A(), l.A() + 3, &coordinates[dim * i]);
#ifdef DEBUGCUTLIBRARY
      pointfile << l(0) << " " << l(1) << " " << l(2) << "\n";
#endif
    }
  }
  else
  {
    for (int i = 0; i < n; ++i)
    {
      Point* p = points[i];
      p->Coordinates(&coordinates[dim * i]);
    }
  }

  boolT ismalloc = false;

  // a set of option we try to process the input with
  // Qz seems to be required for rotational symmetric input
  std::vector<std::string> options;

  // Qhull options:
  //      Qbb - scale last coordinate to [0,m]. Preferable for integer input
  //      QbB - scales input to unit cube [-0.5,0.5]^3. Scaling can reduce precision errors if coord
  //      values wary widely.
#ifdef DIFF_QHULL_CALL
  options.push_back("qhull d Qt QbB Qc Pp");
#endif

  options.push_back("qhull d Qt Qbb Qc Pp");
  options.push_back("qhull d Qt Qbb Qc Qz Pp");
  options.push_back("qhull d Qt Qbb Qc QJ Pp");

  // If you want some debugging information replace the 0 pointer
  // with stdout or some other file open for writing.

#ifdef QHULL_EXTENDED_DEBUG_OUTPUT
  FILE* outfile = stdout;
#else
  FILE* outfile = 0;
#endif
//#define QHULL_DEBUG_OUTPUT
#ifdef QHULL_DEBUG_OUTPUT
#if 0
  static FILE * errfile;
  if ( errfile==NULL )
    errfile = fopen( "qhull_error.log", "w" );
#else
  FILE* errfile = stderr;
#endif
#else
  static NullFile errfile;
#endif

#if 0
  static FILE * debug_f;
  if ( debug_f==NULL )
    debug_f = fopen( "qhull.debug", "w" );
  for ( int i=0; i<n; ++i )
  {
    Point * p = points[i];
    const double * x = p->X();
    fprintf( debug_f, "% .20f % .20f % .20f ", x[0], x[1], x[2] );
  }
  fprintf( debug_f, "\n");
  fflush( debug_f );
#endif
#ifdef QHULL_EXTENDED_DEBUG_OUTPUT
  int counter_qhull = 0;
#endif
  for (std::vector<std::string>::iterator i = options.begin(); i != options.end(); ++i)
  {
#ifdef QHULL_EXTENDED_DEBUG_OUTPUT
    std::cout << "counter_qhull: " << counter_qhull << std::endl;
#endif
    std::string& ostr = *i;
    if (not qh_new_qhull(
            dim, n, &coordinates[0], ismalloc, const_cast<char*>(ostr.c_str()), outfile, errfile))
    {
      // triangulate non-simplicial facets
      qh_triangulate();

      facetT* facet;
      int nf = 0;

      FORALLfacets
      {
        if (not facet->upperdelaunay) nf++;

        // Double check
        if (not facet->simplicial)
        {
          throw std::runtime_error(
              "Qhull returned non-simplicial facets -- try delaunayn with different options");
        }
      }

      tets.reserve(nf);

      FORALLfacets
      {
        if (not facet->upperdelaunay)
        {
          if (facet->vertices)
          {
            std::vector<int> ids;
            ids.reserve(dim + 1);
            vertexT* vertex = NULL;

            // FOREACHvertex_(facet->vertices)
            for (void** vertexp = &facet->vertices->e[0].p;
                 (vertex = static_cast<vertexT*>(*vertexp++));)
            {
              // if delaunayn crashes, enable this check
#if 0
              if (j > dim)
              {
                std::runtime_error("internal error. Qhull returned non-tetsicial facets");
              }
#endif

              int p = qh_pointid(vertex->point);
              if (p >= n)
              {
                throw std::runtime_error("new node in delaunay");
              }
              ids.push_back(p);
            }

            tets.push_back(ids);
          }
        }
      }
    }

    // free long memory
    qh_freeqhull(not qh_ALL);

    // free short memory and memory allocator
    int curlong, totlong;
    qh_memfreeshort(&curlong, &totlong);

    if (curlong or totlong)
    {
      std::stringstream str;
      str << "did not free " << totlong << " bytes of long memory (" << curlong << " pieces)";
      throw std::runtime_error(str.str());
    }

    if (tets.size() > 0)
    {
      plain_int_set used_points;
      for (std::vector<std::vector<int>>::iterator i = tets.begin(); i != tets.end(); ++i)
      {
        std::vector<int>& t = *i;
        std::copy(t.begin(), t.end(), std::inserter(used_points, used_points.begin()));
      }
      if (used_points.size() == points.size())
      {
        return;
      }
      // throw std::runtime_error( "failed to triangulate all points" );

      // failed! start a new iteration.
      tets.clear();
    }
#ifdef QHULL_EXTENDED_DEBUG_OUTPUT
    counter_qhull++;
#endif
  }

#if 1
  // debug output to be read by qhull_test programm
  FILE* debug_f = fopen("qhull.debug", "w");
  for (int i = 0; i < n; ++i)
  {
    Point* p = points[i];
    const double* x = p->X();
    fprintf(debug_f, "% .20f % .20f % .20f ", x[0], x[1], x[2]);
  }
  fprintf(debug_f, "\n");
  fclose(debug_f);
#endif

#ifdef QHULL_DEBUG_OUTPUT
  fflush(errfile);
#endif

  throw std::runtime_error(
      "qhull failed: Maybe the wrong version is used. Check your installation.");
}

/* First check if all the points of a tet share a cut-side, i.e. does the tet lie on a cut-side?
   If not, check if the points share a facet, if not ->         TET IS VALID.
   If yes, check if the cut_side is a LevelSetSide, if yes ->   TET IS VALID.
           if not LevelSetSide,
              check if the points share a facet, if no ->       TET IS VALID.
                   if yes,
                      check if facet is triangulated, if yes -> TET IS VALID.
                                                      if no  -> TET IS INVALID.

  + This was the way envisioned by Kuettler. However, it makes little sense accepting a TET
 completely on a cut-side. A simplification by me (Magnus) is to remove TETS whos points are all on
 a cut-side or more than one cut-side (A case not uncommon for MeshIntersection, however this works?
 Would it imply a TET is on a line?).
  + The change proposed does not change any existing test-case or any cut-test. However, the
 algorithm implemented, by Kuettler, is probably not there unnecessarily (hopefully?). Thus keep the
 code for reference in the future.
 !+ For now use the old way of Kuettler as after further investigation it seems to be somewhat more
 exact...

*/
bool GEO::CUT::TetMesh::IsValidTet(const std::vector<Point*>& t)
{
  plain_side_set sides;
  // Find if the points of the tet share a common side.
  FindCommonSides(t, sides);
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
  std::cout << "sides.size(): " << sides.size() << std::endl;
  for (unsigned k2 = 0; k2 < sides.size(); k2++)
  {
    std::cout << "#" << k2 << ": ";
    std::cout << "IsLevelSetSide: " << sides[k2]->IsLevelSetSide() << std::endl;
    sides[k2]->Print();
    std::cout << std::endl;
  }
#endif

  if (sides.size() == 0)
  {
    plain_facet_set facets;
    // Find if the points of the tet share a common facet.
    FindCommonFacets(t, facets);
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << "facets.size(): " << facets.size() << std::endl;
#endif
    if (facets.size() == 0)
    {
      return true;
    }
    // Why do we have to enter here? Shouldn't it be clear already if the side=0, the points of the
    // tet
    //  can't share a common facet? A cut-side can be outside of the element, thus a case can occur,
    //  when a tet is not on a cut-side. BUT shares a facet. However this is already tested in
    //  FindCommonSides, as the points know what sides it cuts.
#ifdef DEBUGCUTLIBRARY
    std::cout << "t.size(): " << t.size() << std::endl;
    for (unsigned k = 0; k < t.size(); k++)
    {
      std::cout << "t[" << k << "]->Id()" << t[k]->Id() << std::endl;
      std::cout << "t[k]->CutSides().size(): " << t[k]->CutSides().size() << std::endl;
      for (unsigned k2 = 0; k2 < t[k]->CutSides().size(); k2++)
      {
        std::cout << "#" << k2 << ": ";
        std::cout << "IsLevelSetSide: " << t[k]->CutSides()[k2]->IsLevelSetSide() << std::endl;
        t[k]->CutSides()[k2]->Print();
        std::cout << std::endl;
      }
      for (unsigned k2 = 0; k2 < t[k]->Facets().size(); k2++)
      {
        std::cout << "#" << k2 << ": ";
        // std::cout << "IsLevelSetSide: " << t[k]->CutSides()[k2]->IsLevelSetSide() << std::endl;
        t[k]->Facets()[k2]->PrintPointIds();
        std::cout << std::endl;
      }
    }
#endif
    dserror(
        "You have encountered a case where the a tet is not on a cut-side BUT is on a facet!!! "
        "CHECK IT!");

    for (plain_facet_set::iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (not f->IsTriangulated())
      {
        return false;
      }
    }
    return true;
  }  // end if(sides.size()==0)

  // Why would we want to accept a tet completely on a LevelSet-side?
  // POSSIBLE ANS: Might be a cut with more than 3nodes (i.e. 4,5,6) nodes NOT in a plane
  //               and tet defined by these 4 points?
  if (sides.size() == 1)
  {
    // Is this side a LevelSetSide, i.e. does the cut lie on a LevelSetSide?
    // The points of a tet can all be on a cut-side that is a LevelSetSide, but not share a common
    // facet. This occurs when there exists a degenerate cut (i.e. a cut we can't cut well). One VC
    // is created and that's about all.... Might have to investigate...
    if ((*sides.begin())->IsLevelSetSide())
    {
      // This part is done mainly for debugging purposes. It is probably not necessary,
      //  but could play an important role when/if the LevelSetSide is remodeled.
      plain_facet_set facets;
      // Find if the points of the tet share a common facet.
      FindCommonFacets(t, facets);

      // Eventhough all nodes share a level set side.
      if (facets.size() != 1)
      {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
        std::cout << "side.size()==1 and LEVELSetSide with facets.size(): " << facets.size() << "!"
                  << std::endl;
        std::cout << "This is most likely a degenerate cut case for the level-set. Something "
                     "should have been done before this stage...."
                  << std::endl;
#endif
        //        std::cout << "facets.size(): " << facets.size();
        //        std::cout << "t.size(): " << t.size() << std::endl;
        //        for(unsigned k=0; k<t.size(); k++)
        //        {
        //          std::cout << "t[" << k << "]->Id()" << t[k]->Id() << std::endl;
        //          std::cout << "t[k]->CutSides().size(): " << t[k]->CutSides().size() <<
        //          std::endl; for(unsigned k2=0; k2<t[k]->CutSides().size(); k2++)
        //          {
        //            std::cout << "#"<< k2 << ": ";
        //            std::cout << "IsLevelSetSide: " << t[k]->CutSides()[k2]->IsLevelSetSide() <<
        //            std::endl; t[k]->CutSides()[k2]->Print(); std::cout << std::endl;
        //          }
        //          for(unsigned k2=0; k2<t[k]->Facets().size(); k2++)
        //          {
        //            std::cout << "#"<< k2 << ": ";
        //            t[k]->Facets()[k2]->PrintPointIds();
        //            std::cout << std::endl;
        //          }
        //        }
        // throw std::runtime_error("A LevelSetSide should BE associated to ONE facet. CHECK
        // THIS!");
      }

      for (plain_facet_set::iterator i = facets.begin(); i != facets.end(); ++i)
      {
        Facet* f = *i;
        if (f->IsTriangulated())
        {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
          std::cout << "sides.size()==1 and LEVELSetSide with Triangulated facets!" << std::endl;
#endif
          return true;
        }
        else
          dserror("Tet completely on LevelSetSide, and NOT Triangulated.");
      }
      return true;
    }
  }  // end if(sides.size()==1)

  // If a TET lies completely on a cut-side, remove this tet.
  // This is not the way of Kuettler!!! It is unclear what his vision for this was...
#ifdef REMOVE_ALL_TETS_ON_CUTSIDE
  if (sides.size() > 0)
  {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << "One TET lies on MORE THAN one cut-side. CUT-SIDES OVERLAPPING! Discard this tet."
              << std::endl;
#endif
    return false;
  }
#else
  plain_facet_set facets;
  FindCommonFacets(t, facets);
  if (facets.size() == 0)
  {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << "sides.size()==" << sides.size() << ", and facets.size()==0" << std::endl;
#endif
    return true;
  }
  for (plain_facet_set::iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    if (f->IsTriangulated())  //(i.e. is this facet not a tri?)
    {
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      std::cout << "sides.size()==" << sides.size() << ", and f->IsTriangulated()" << std::endl;
#endif
      return true;
    }
  }
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
  std::cout << "sides.size()==" << sides.size() << ", and facets.size(): " << facets.size()
            << std::endl;
  for (unsigned k2 = 0; k2 < facets.size(); k2++)
  {
    std::cout << "#" << k2 << ": ";
    facets[k2]->PrintPointIds();
    std::cout << std::endl;
  }
#endif
#endif
  return false;
}

#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
double GEO::CUT::TetMesh::CalcVolumeOfTet(const std::vector<Point*>& t)
{
  if (t.size() != 4) dserror("Expected a tet. Size of vector is not 4.");
  // create planes consisting of 3 nodes each
  LINALG::Matrix<3, 1> p0(t[0]->X());
  LINALG::Matrix<3, 1> p1(t[1]->X());
  LINALG::Matrix<3, 1> p2(t[2]->X());
  LINALG::Matrix<3, 1> p3(t[3]->X());

  LINALG::Matrix<3, 1> v01;
  LINALG::Matrix<3, 1> v02;
  LINALG::Matrix<3, 1> v03;

  v01.Update(1, p1, -1, p0, 0);
  v02.Update(1, p2, -1, p0, 0);
  v03.Update(1, p3, -1, p0, 0);

  // create 4 normal vectors to each tet surface plane
  LINALG::Matrix<3, 1> nplane012;

  // cross product
  nplane012(0) = v01(1) * v02(2) - v01(2) * v02(1);
  nplane012(1) = v01(2) * v02(0) - v01(0) * v02(2);
  nplane012(2) = v01(0) * v02(1) - v01(1) * v02(0);

  // compute normal distance of point to plane of the three remaining points
  double distance = nplane012.Dot(v03);

  double vol_tet = distance / 6.0;

  return vol_tet;
}
#endif

/* This function is unused....
 */
void GEO::CUT::TetMesh::TestUsedPoints(const std::vector<std::vector<int>>& tets)
{
  plain_int_set used_points;
  for (std::vector<std::vector<int>>::const_iterator i = tets.begin(); i != tets.end(); ++i)
  {
    const std::vector<int>& t = *i;
    std::copy(t.begin(), t.end(), std::inserter(used_points, used_points.begin()));
  }
  if (used_points.size() != points_.size())
  {
    throw std::runtime_error("failed to triangulate all points");
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::TetMesh::FixBrokenTets()
{
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
  int counter = 0;
#endif
  for (std::vector<std::vector<int>>::iterator i = tets_.begin(); i != tets_.end(); ++i)
  {
    std::vector<int>& t = *i;

    // create planes consisting of 3 nodes each
    LINALG::Matrix<3, 1> p0(points_[t[0]]->X());
    LINALG::Matrix<3, 1> p1(points_[t[1]]->X());
    LINALG::Matrix<3, 1> p2(points_[t[2]]->X());
    LINALG::Matrix<3, 1> p3(points_[t[3]]->X());

    LINALG::Matrix<3, 1> v01;
    LINALG::Matrix<3, 1> v02;
    LINALG::Matrix<3, 1> v03;

    v01.Update(1, p1, -1, p0, 0);
    v02.Update(1, p2, -1, p0, 0);
    v03.Update(1, p3, -1, p0, 0);

    // create 4 normal vectors to each tet surface plane
    LINALG::Matrix<3, 1> nplane012;

    // cross product
    nplane012(0) = v01(1) * v02(2) - v01(2) * v02(1);
    nplane012(1) = v01(2) * v02(0) - v01(0) * v02(2);
    nplane012(2) = v01(0) * v02(1) - v01(1) * v02(0);

    // compute normal distance of point to plane of the three remaining points
    double distance = nplane012.Dot(v03);

    // compute norm (area) of plane
    // double norm012 = nplane012.Norm2();

    LINALG::Matrix<4, 1> temp(true);
    temp(0, 0) = p0.Norm2();  // Distance of points to origin
    temp(1, 0) = p1.Norm2();
    temp(2, 0) = p2.Norm2();
    temp(3, 0) = p3.Norm2();

    // This is to scale the tolerance, it determines our maximum precision (i.e. machine precision)
    double max_dist_to_orgin = temp.NormInf();

    LINALG::Matrix<3, 1> v04;
    v04.Update(1, p1, -1, p2, 0);

    temp(0, 0) = v01.Norm2();  // Distance between points in "base" triangle
    temp(1, 0) = v02.Norm2();
    temp(2, 0) = v04.Norm2();
    temp(3, 0) = fabs(distance);  //"Height" of tetrahedral

#ifdef NEW_POSTOL_TET
    // This is the smallest distance between the base points in the tet and the height.
    double min_dist_in_tet = temp.MinValue();
    // We want to test with this one I think... But might lead to problems.
    double tolerance = LINSOLVETOL * max_dist_to_orgin;
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    double vol_tet = distance / 6.0;
#endif
#else
    double vol_tet = distance / 6.0;
#ifdef DEBUGCUTLIBRARY
    // This is the smallest distance between the base points in the tet and the height.
    double min_dist_in_tet = temp.MinValue();
    // We want to test with this one I think... But might lead to problems.
    double tolerance = LINSOLVETOL * max_dist_to_orgin;
#endif
#endif

#if TETMESH_EXTENDED_DEBUG_OUTPUT
    std::cout << "=================================" << std::endl;
    std::cout << "tolerance: " << tolerance << std::endl;
    std::cout << "min_dist_in_tet: " << min_dist_in_tet << std::endl;
#endif

    // Deactivate all tets that are too small. We might still need the tet to
    // create a cut surface tri. Afterwards we will discard it.
#ifdef NEW_POSTOL_TET
    if (min_dist_in_tet < tolerance)
    {
      accept_tets_[i - tets_.begin()] = false;
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      std::cout << "Removed tet(" << counter << ") , with volume: " << vol_tet << std::endl;
#endif
    }
#else
    if (fabs(vol_tet) < VOLUMETOL)
    {
      accept_tets_[i - tets_.begin()] = false;
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
      std::cout << "Removed tet(" << counter << ") , with volume: " << vol_tet << std::endl;
#endif
    }
    //     else if ( fabs( distance / norm012 ) < 1e-7 )
    //     {
    //       accept_tets_[i - tets_.begin()] = false;
    //     }
#endif


#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    // Is the TET accepted/rejected in the new framework but rejected/accepted in the old?
    if ((min_dist_in_tet < tolerance) != (fabs(vol_tet) < VOLUMETOL))
    {
      std::cout << "=================================" << std::endl;
      std::cout << "tolerance: " << tolerance << std::endl;
      std::cout << "min_dist_in_tet: " << min_dist_in_tet << std::endl;
      std::cout << "VOLUMETOL: " << VOLUMETOL << std::endl;
      std::cout << "volume: " << vol_tet << std::endl;

      if (min_dist_in_tet < tolerance)
        std::cout << " TET is NOT ACCEPTED in new tolerance." << std::endl;
      else
        std::cout << " TET is ACCEPTED in new tolerance." << std::endl;
    }
#endif

    // tet numbering wrong exchange 1 with 3
    if (distance < 0)
    {
      std::swap(t[1], t[3]);
    }
#ifdef TETMESH_EXTENDED_DEBUG_OUTPUT
    counter++;
#endif
  }
}

/* Take the tri from the facet and test whether it belongs to more than one tet,
   if not add the tri as a side with an error check.
 */
void GEO::CUT::TetMesh::FindProperSides(const PlainEntitySet<3>& tris,
    std::vector<std::vector<int>>& sides, const PlainEntitySet<4>* members)
{
  sides.reserve(tris.size());
  for (PlainEntitySet<3>::const_iterator i = tris.begin(); i != tris.end(); ++i)
  {
    Entity<3>* tri = *i;

    bool done = false;
    const std::vector<Entity<4>*>& tets = tri->Parents();
    for (std::vector<Entity<4>*>::const_iterator i = tets.begin(); i != tets.end(); ++i)
    {
      Entity<4>* tet = *i;

      // Test if the tet exist as member in the volume-cell, if not skip this tet.
      if (members != NULL and members->count(tet) == 0)
      {
        continue;
      }

      // Get the tet from the set of original tets.
      //  QUESTION: Why do we do this?
      //  Possible ans: Want access to the node-IDs.
      std::vector<int>& original_tet = tets_[tet->Id()];

      if (original_tet.size() > 0)
      {
        if (done)
        {
          throw std::runtime_error("double tets at cut surface");
        }

        done = true;

        bool found = false;
        for (int i = 0; i < 4; ++i)
        {
          std::vector<int> side(3);
          for (int j = 0; j < 3; ++j)
          {
            side[j] = original_tet[DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j]];
          }
          if (tri->Equals(side))
          {
            found = true;
            sides.push_back(std::vector<int>());
            std::swap(sides.back(), side);
            break;
          }
        }
        if (not found)
        {
          throw std::runtime_error("failed to find side");
        }
      }
    }
    if (not done)
    {
      throw std::runtime_error("failed to find tet");
    }
  }
}

/// Collects the coordinates for the tri3 sides of the facet if all its points are on cut surface
void GEO::CUT::TetMesh::CollectCoordinates(
    const std::vector<std::vector<int>>& sides, std::vector<Point*>& side_coords)
{
  for (std::vector<std::vector<int>>::const_iterator i = sides.begin(); i != sides.end(); ++i)
  {
    const std::vector<int>& side = *i;

    Point* p1 = points_[side[0]];
    Point* p2 = points_[side[1]];
    Point* p3 = points_[side[2]];

    // Plausible: Might not need this check.
    //  However, left as is, as Tetmeshintersection is still not well understood.
    if (p1->HasAssociatedBoundaryCellFacet() and p2->HasAssociatedBoundaryCellFacet() and
        p3->HasAssociatedBoundaryCellFacet())
    {
      side_coords.push_back(p1);
      side_coords.push_back(p2);
      side_coords.push_back(p3);

#ifdef TETMESH_GMSH_DEBUG_OUTPUT
      surface_tris_.push_back(side);
#endif
    }

#ifdef DEBUGCUTLIBRARY
    else
    {
      throw std::runtime_error("Side not on cut or marked surface!!! Shouldn't it be?");
    }
#endif
  }
}

#ifdef TETMESH_GMSH_DEBUG_OUTPUT

void GEO::CUT::TetMesh::GmshWriteCells()
{
  std::ofstream file("delaunaycells.pos");
  file << "View \"delaunaycells\" {\n";
  for (std::vector<std::vector<int>>::iterator i = tets_.begin(); i < tets_.end(); ++i)
  {
    std::vector<int>& t = *i;
    GmshWriteTet(file, i - tets_.begin(), t);
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteActiveCells()
{
  std::ofstream file("activecells.pos");
  file << "View \"activecells\" {\n";
  for (std::vector<std::vector<int>>::iterator i = tets_.begin(); i < tets_.end(); ++i)
  {
    unsigned pos = i - tets_.begin();
    if (accept_tets_[pos])
    {
      std::vector<int>& t = *i;
      GmshWriteTet(file, i - tets_.begin(), t);
    }
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteSurfaceCells()
{
  std::ofstream file("surfacecells.pos");
  file << "View \"surfacecells\" {\n";
  for (std::map<Handle<3>, Entity<3>>::iterator i = tet_surfaces_.begin(); i != tet_surfaces_.end();
       ++i)
  {
    Entity<3>& tri = i->second;
    std::vector<int> t(tri(), tri() + 3);
    GmshWriteTri(file, tri.Id(), t);
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteSurfaceTris()
{
  std::ofstream file("surfacetris.pos");
  file << "View \"surfacetris\" {\n";
  for (std::vector<std::vector<int>>::iterator i = surface_tris_.begin(); i < surface_tris_.end();
       ++i)
  {
    std::vector<int>& t = *i;
    GmshWriteTri(file, i - surface_tris_.begin(), t);
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteTriSet(const std::string& name, const PlainEntitySet<3>& tris)
{
  std::string filename = name + ".pos";
  std::ofstream file(filename.c_str());
  file << "View \"" << name << "\" {\n";
  for (PlainEntitySet<3>::const_iterator i = tris.begin(); i != tris.end(); ++i)
  {
    Entity<3>& tri = **i;
    std::vector<int> t(tri(), tri() + 3);
    GmshWriteTri(file, tri.Id(), t);
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteTetSet(const std::string& name, const PlainEntitySet<4>& tets)
{
  std::string filename = name + ".pos";
  std::ofstream file(filename.c_str());
  file << "View \"" << name << "\" {\n";
  for (PlainEntitySet<4>::const_iterator i = tets.begin(); i != tets.end(); ++i)
  {
    Entity<4>& tet = **i;
    std::vector<int> t(tet(), tet() + 4);
    GmshWriteTet(file, tet.Id(), t);
  }
  file << "};\n";
}

void GEO::CUT::TetMesh::GmshWriteTri(std::ostream& file, int eid, const std::vector<int>& t)
{
  GmshWriteConnect(file, "ST", t);
  GmshWritePosition(file, eid, t);
}

void GEO::CUT::TetMesh::GmshWriteTet(std::ostream& file, int eid, const std::vector<int>& t)
{
  GmshWriteConnect(file, "SS", t);
  GmshWritePosition(file, eid, t);
}

void GEO::CUT::TetMesh::GmshWriteConnect(
    std::ostream& file, std::string name, const std::vector<int>& t)
{
  file << name << "(";
  for (std::vector<int>::const_iterator j = t.begin(); j != t.end(); ++j)
  {
    Point* p = points_[*j];
    if (j != t.begin()) file << ",";
    file << p->X()[0] << "," << p->X()[1] << "," << p->X()[2];
  }
  file << "){";
}

void GEO::CUT::TetMesh::GmshWritePosition(std::ostream& file, int eid, const std::vector<int>& t)
{
  for (std::vector<int>::const_iterator j = t.begin(); j != t.end(); ++j)
  {
    Point* p = points_[*j];
    if (j != t.begin()) file << ",";
    file << p->Position();
  }
  file << "};  // " << eid << "\n";
}

#endif
