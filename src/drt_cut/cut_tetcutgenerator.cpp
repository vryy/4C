
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

    mesh.AddTetgen( out );

    cut_mesh_.Cut( mesh );

    mesh.MakeFacets();
    mesh.FindNodePositions();

    mesh.Status();

    if ( parent_ != NULL )
    {
      Mesh resultmesh( pp_, false );
      TetCutConverter converter( resultmesh );
      mesh.GenerateTetgen( &converter );

      tetgenio result;
      resultmesh.ExtractTetgen( result );
      parent_->Generate( element, result );
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

GEO::CUT::TetCutConverter::TetCutConverter( Mesh & mesh )
  : mesh_( mesh )
{
}

void GEO::CUT::TetCutConverter::NonCut( Element* element )
{
#ifdef QHULL
  tetgenio out;
  element->FillTetgen( out );
  mesh_.AddTetgen( out );
#endif
}

#ifdef QHULL
void GEO::CUT::TetCutConverter::Generate( Element* e, const tetgenio & out )
{
  mesh_.AddTetgen( out );
}
#endif

