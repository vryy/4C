/*----------------------------------------------------------------------------*/
/*! \file

\brief Create a simple facet graph for 1-D and 2-D elements ( embedded in a
       higher dimensional space )

\date Nov 14, 2016

\level 2

 *------------------------------------------------------------------------------------------------*/
#include "cut_facetgraph_simple.H"

#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_side.H"
#include "cut_pointgraph_simple.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::SimpleFacetGraph_1D::SimpleFacetGraph_1D(
    const std::vector<Side*>& sides, const plain_facet_set& facets)
    : FacetGraph() /* call empty base class constructor */
{
  // loop over all facets
  for (plain_facet_set::const_iterator cit = facets.begin(); cit != facets.end(); ++cit)
  {
    Facet* f = *cit;

    if (f->HasHoles()) run_time_error("Holes are not yet considered for the simple case!");
  }

  // first store all facets
  all_facets_.reserve(facets.size());
  all_facets_.assign(facets.begin(), facets.end());

#if 0
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  for ( std::vector<Facet * >::const_iterator cit = all_facets_.begin();
        cit != all_facets_.end(); ++cit )
    (*cit)->Print();
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::SimpleFacetGraph_1D::CreateVolumeCells(
    Mesh& mesh, Element* element, plain_volumecell_set& cells)
{
  if (all_facets_.size() < 1) return;

  // sort the facets by their local line coordinate
  std::map<double, Facet*> sorted_facets;
  SortFacets(element, sorted_facets);

  // a pair of two consecutive facets forms a pseudo volume cell in 1D
  std::vector<plain_facet_set> volumes;
  volumes.reserve(all_facets_.size() - 1);
  CombineFacetsToLineVolumes(sorted_facets, volumes);

#if 0
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  unsigned vcount = 0;
  for ( std::vector<plain_facet_set>::const_iterator c_vec = volumes.begin();
      c_vec != volumes.end(); ++c_vec )
  {
    std::cout << "--- Volume " << vcount++ << std::endl;
    for ( plain_facet_set::const_iterator c_f = c_vec->begin();
        c_f != c_vec->end(); ++c_f )
    {
      (*c_f)->Print();
    }
  }
#endif

  // no volume lines in 1-D
  std::map<std::pair<Point*, Point*>, plain_facet_set> volume_lines;

  for (std::vector<plain_facet_set>::iterator it = volumes.begin(); it != volumes.end(); ++it)
  {
    plain_facet_set& collected_facets = *it;

    // add the pseudo volume cell
    cells.insert(mesh.NewVolumeCell(collected_facets, volume_lines, element));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::SimpleFacetGraph_1D::SortFacets(
    const Element* element, std::map<double, Facet*>& sorted_facets) const
{
  Edge* edge = element->Sides()[0]->Edges()[0];

  for (std::vector<Facet*>::const_iterator cit = all_facets_.begin(); cit != all_facets_.end();
       ++cit)
  {
    Facet* f = *cit;
    if (not f->Equals(DRT::Element::point1)) dserror("The given facets are supposed to be points!");

    Point* p = f->Points()[0];
    double t = p->t(edge);

    sorted_facets[t] = f;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::SimpleFacetGraph_1D::CombineFacetsToLineVolumes(
    const std::map<double, Facet*>& sorted_facets, std::vector<plain_facet_set>& volumes) const
{
  std::map<double, Facet*>::const_iterator end_cit = sorted_facets.end();
  --end_cit;
  for (std::map<double, Facet*>::const_iterator cit = sorted_facets.begin(); cit != end_cit;)
  {
    // get the line end points
    Facet* f1 = cit->second;
    Facet* f2 = (++cit)->second;

    // insert the line end points into the line volumes
    volumes.push_back(plain_facet_set());
    plain_facet_set& line = volumes.back();

    line.insert(f1);
    line.insert(f2);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::SimpleFacetGraph_2D::SimpleFacetGraph_2D(
    const std::vector<Side*>& sides, const plain_facet_set& facets)
    : FacetGraph() /* call empty base class constructor */
{
  // loop over all facets
  for (plain_facet_set::const_iterator cit = facets.begin(); cit != facets.end(); ++cit)
  {
    Facet* f = *cit;

    if (f->HasHoles()) run_time_error("Holes are not yet considered for the simple case!");
  }

  // first store all facets
  all_facets_.reserve(facets.size());
  all_facets_.assign(facets.begin(), facets.end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::SimpleFacetGraph_2D::CreateVolumeCells(
    Mesh& mesh, Element* element, plain_volumecell_set& cells)
{
  std::vector<plain_facet_set> volumes;

  const std::vector<Side*> sides = element->Sides();
  if (sides.size() != 1) dserror("A 2-D element is supposed to contain exactly one side!");

  Side* side = sides[0];

  IMPL::SimplePointGraph_2D pg_2d(
      mesh, element, side, IMPL::PointGraph::element_side, IMPL::PointGraph::all_lines);

  volumes.reserve(pg_2d.NumSurfaces());

  for (IMPL::SimplePointGraph_2D::surface_const_iterator sit = pg_2d.sbegin(); sit != pg_2d.send();
       ++sit)
  {
    Cycle vol_point_cycle = *sit;
    const unsigned cycle_length = vol_point_cycle().size();

    volumes.push_back(plain_facet_set());
    plain_facet_set& collected_facets = volumes.back();

    std::vector<Point*> line_points(2, NULL);
    for (unsigned i = 0; i < cycle_length; ++i)
    {
      line_points[0] = vol_point_cycle()[i % cycle_length];
      line_points[1] = vol_point_cycle()[(i + 1) % cycle_length];

      Facet* f = FindFacet(all_facets_, line_points);
      if (not f)
      {
        // debug output
        std::cout << "=====================================================" << std::endl;
        std::cout << "SimpleFacetGraph_2D::FindFacet failed! ( line = " << __LINE__ << " )\n";
        line_points[0]->Print(std::cout);
        line_points[1]->Print(std::cout);
        std::cout << "=====================================================" << std::endl;

        dserror("Could not find the correct line facet!");
      }

      collected_facets.insert(f);
    }
  }

  // finally add the volumes to volume cells set
  AddToVolumeCells(mesh, element, volumes, cells);
}
