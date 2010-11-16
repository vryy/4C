
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "cut_tetcutgenerator.H"
#include "cut_boundingbox.H"
#include "cut_element.H"

GEO::CUT::TetCutGenerator::TetCutGenerator( GEO::CUT::MeshIntersection & intersection,
                                            CellGenerator * parent,
                                            Mesh & cut_mesh,
                                            Teuchos::RCP<PointPool> pp )
  : intersection_( intersection ),
    parent_( parent ),
    cut_mesh_( cut_mesh ),
    pp_( pp )
{
}

void GEO::CUT::TetCutGenerator::NonCut( Element* element )
{
  std::cout << "uncut element!\n";
}

#ifdef QHULL
void GEO::CUT::TetCutGenerator::Generate( Element* element, const tetgenio & out )
{
  if ( cut_mesh_.WithinBB( *element ) )
  {
    const_cast<tetgenio &>( out ).save_nodes( const_cast<char*>( "tetcut" ) );
    const_cast<tetgenio &>( out ).save_elements( const_cast<char*>( "tetcut" ) );
    const_cast<tetgenio &>( out ).save_faces( const_cast<char*>( "tetcut" ) );

    Mesh mesh( pp_, false );

    // need to copy surface markers via side id
    for ( int i=0; i<out.numberoftrifaces; ++i )
    {
      if ( out.trifacemarkerlist[i] > -1 )
      {
        //GEO::CUT::Side * cut_side = intersection_.GetCutSides( out.trifacemarkerlist[i] );
        std::vector<int> nids;
        nids.reserve( 3 );
        for ( int j=0; j<3; ++j )
        {
          int pointidx = out.trifacelist[i*3+j] * 3;
          Node* n = mesh.GetNode( pointidx, &out.pointlist[pointidx] );
          nids.push_back( pointidx );
        }
        Side * side = mesh.CreateTri3( out.trifacemarkerlist[i], nids );
      }
    }

    //const int numTetNodes = DRT::UTILS::getNumberOfElementNodes(
    //DRT::Element::tet4 );
    const int numTetNodes = 4;

    for ( int i=0; i<out.numberoftetrahedra; ++i )
    {
      std::vector<int> nids;
      nids.reserve( numTetNodes );
      for ( int j=0; j<numTetNodes; ++j )
      {
        int pointidx = out.tetrahedronlist[i*out.numberofcorners+j] * 3;
        Node* n = mesh.GetNode( pointidx, &out.pointlist[pointidx] );
        nids.push_back( pointidx );
      }
      Element * e = mesh.CreateTet4( i, nids );
    }

    cut_mesh_.Cut( mesh );

    mesh.MakeFacets();
    mesh.FindNodePositions();

    mesh.Status();

    if ( parent_ != NULL )
    {
      TetCutConverter converter( parent_, element );
      mesh.GenerateTetgen( &converter );
    }
  }
  else
  {
    if ( parent_ != NULL )
    {
      parent_->Generate( element, out );
    }
  }
}
#endif

GEO::CUT::TetCutConverter::TetCutConverter( CellGenerator * parent, Element * element )
  : parent_( parent ),
    element_( element )
{
}

void GEO::CUT::TetCutConverter::NonCut( Element* element )
{
#ifdef QHULL
  tetgenio out;
  element->FillTetgen( out );
  Generate( element, out );
#endif
}

#ifdef QHULL
void GEO::CUT::TetCutConverter::Generate( Element* e, const tetgenio & out )
{
  // set original element
  parent_->Generate( element_, out );
}
#endif

