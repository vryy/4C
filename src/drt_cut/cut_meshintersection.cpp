
#include "cut_meshintersection.H"
#include "cut_tetcutgenerator.H"

void GEO::CUT::MeshIntersection::AddElement( int eid,
                                             const std::vector<int> & nids,
                                             const Epetra_SerialDenseMatrix & xyz,
                                             DRT::Element::DiscretizationType distype )
{
  for ( std::vector<Teuchos::RCP<Mesh> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    Mesh & cut_mesh = **i;
    if ( cut_mesh.WithinBB( xyz ) )
    {
      int numnode = nids.size();
      if ( numnode != xyz.N() )
      {
        throw std::runtime_error( "node coordiante number mismatch" );
      }

      // make sure all nodes are there
      for ( int i=0; i<numnode; ++i )
      {
        mesh_.GetNode( nids[i], &xyz( 0, i ) );
      }

      // create element
      mesh_.CreateElement( eid, nids, distype );

      return;
    }
  }
}

void GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                             const std::vector<int> & nids,
                                             DRT::Element::DiscretizationType distype,
                                             int mi )
{
  // create side
  Mesh & cut_mesh = CutMesh( mi );
  cut_mesh.CreateSide( sid, nids, distype );
}

void GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                             const std::vector<int> & nids,
                                             const Epetra_SerialDenseMatrix & xyz,
                                             DRT::Element::DiscretizationType distype,
                                             int mi )
{
  Mesh & cut_mesh = CutMesh( mi );

  int numnode = nids.size();
  if ( numnode != xyz.N() )
  {
    throw std::runtime_error( "node coordiante number mismatch" );
  }

  // make sure all nodes are there
  for ( int i=0; i<numnode; ++i )
  {
    cut_mesh.GetNode( nids[i], &xyz( 0, i ) );
  }

  // create side
  cut_mesh.CreateSide( sid, nids, distype );
}

void GEO::CUT::MeshIntersection::Cut( CellGenerator * generator )
{
  for ( std::vector<Teuchos::RCP<Mesh> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    Mesh & cut_mesh = **i;
    cut_mesh.FillComplete();
  }
  mesh_.FillComplete();

  //Status();

  std::vector<Teuchos::RCP<CellGenerator> > cutgens;
  for ( int i=cut_mesh_.size()-1; i>0; --i )
  {
    cutgens.push_back( Teuchos::rcp( new TetCutGenerator( *this, generator, CutMesh( i ), pp_ ) ) );
    generator = &*cutgens.back();
  }

  // loop cut sides and cut against elements at the same position in space
  CutMesh().Cut( mesh_ );

  mesh_.MakeFacets();
  mesh_.FindNodePositions();

  //Status();

  mesh_.GenerateTetgen( generator );
}

void GEO::CUT::MeshIntersection::SelfCut()
{
  CutMesh().SelfCut();
}

GEO::CUT::Side * GEO::CUT::MeshIntersection::GetCutSides( int sid, int mi )
{
  const std::vector<GEO::CUT::Side*> & cut_sides = CutMesh( mi ).GetSides( sid );
  if ( cut_sides.size()==1 )
  {
    return cut_sides[0];
  }
  if ( cut_sides.size()>1 )
  {
    // the first entry has to be the quadratic element, the following entries
    // its linear parts
    return cut_sides[0];
  }
  throw std::runtime_error( "no such side" );
}

void GEO::CUT::MeshIntersection::Status()
{
  mesh_.Status();
  for ( std::vector<Teuchos::RCP<Mesh> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    Mesh & cut_mesh = **i;
    cut_mesh.Status();
  }

  mesh_.DumpGmsh( "mesh" );
  int count = 0;
  for ( std::vector<Teuchos::RCP<Mesh> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    Mesh & cut_mesh = **i;
    std::stringstream str;
    str << "cut_mesh" << count;
    cut_mesh.DumpGmsh( str.str().c_str() );
    count++;
  }
}
