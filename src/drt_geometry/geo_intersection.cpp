/*!-----------------------------------------------------------------------------------------------*
\file geo_intersection.cpp

\brief class that provides to set up a mesh cut based on a level set field or on further surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_Time.hpp>

#include "../drt_fluid/xfluid_defines.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_facet.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_side.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_parallel.H"

#include "geo_utils.H"

#include "geo_intersection.H"


GEO::CutWizard::CutWizard( DRT::Discretization & dis, bool include_inner, int numcutmesh )
  : dis_( dis ),
    include_inner_( include_inner )
{
  mesh_ = Teuchos::rcp( new GEO::CUT::MeshIntersection( numcutmesh ) );
}

void GEO::CutWizard::SetFindPositions( bool positions )
{
  mesh_->SetFindPositions( positions );
}

void GEO::CutWizard::AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze )
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );
  mesh_->AddCutSide( ele->Id(), nids, xyze, ele->Shape(), mi );
}

void GEO::CutWizard::AddElement( DRT::Element * ele )
{
  const int numnode = ele->NumNode();
  const DRT::Node * const * nodes = ele->Nodes();
  const int * nodeids = ele->NodeIds();

  Epetra_SerialDenseMatrix xyze( 3, numnode );

  for ( int i=0; i < numnode; ++i )
  {
    const DRT::Node & node = *nodes[i];
    std::copy( node.X(), node.X()+3, &xyze( 0, i ) );
  }

  std::vector<int> nids( nodeids, nodeids+numnode );
  mesh_->AddElement( ele->Id(), nids, xyze, ele->Shape() );
}

GEO::CUT::Side * GEO::CutWizard::GetSide( std::vector<int>& nodeids )
{
  return mesh_->GetSide( nodeids );
}

GEO::CUT::SideHandle * GEO::CutWizard::GetSide( int sid )
{
  return mesh_->GetSide( sid );
}

GEO::CUT::SideHandle * GEO::CutWizard::GetCutSide( int sid, int mi )
{
  return mesh_->GetCutSide( sid, mi );
}

GEO::CUT::ElementHandle * GEO::CutWizard::GetElement( DRT::Element * ele )
{
  return mesh_->GetElement( ele->Id() );
}

GEO::CUT::Node * GEO::CutWizard::GetNode( int nid )
{
  return mesh_->GetNode( nid );
}




/*------------------------------------------------------------------------------------------------*
 * cut routine for parallel framework in XFSI and XFLUIDFLUID                        schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::CutParallel( bool include_inner, std::string VCellgausstype, std::string BCellgausstype )
{
  // for XFSI and XFLUIDFLUID we have communicate node positions and dofset data
  bool communicate = true;

  if(dis_.Comm().NumProc() == 1) communicate = false;

  mesh_->Status();


  // FIRST step: cut the mesh
  mesh_->Cut_Mesh( include_inner );

  // SECOND step: find node positions and create dofset in PARALLEL
  CutParallel_FindPositionDofSets( include_inner, communicate );

  const double t_start = Teuchos::Time::wallTime();
  // THIRD step: perform tessellation or moment fitting on the mesh
  mesh_->Cut_Finalize( include_inner, VCellgausstype, BCellgausstype );

  const double t_end = Teuchos::Time::wallTime()-t_start;
  //if ( backdis_.Comm().MyPID() == 0 )
  //{
    std::cout << "\n XFEM::FluidWizard::Quadrature construction time = " <<t_end<<"\n";
  //}


  mesh_->Status(VCellgausstype);
}


/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::CutParallel_FindPositionDofSets(bool include_inner, bool communicate)
{
  GEO::CUT::Options options;
  mesh_->GetOptions(options);

  if ( options.FindPositions() )
  {

    GEO::CUT::Mesh & m = mesh_->NormalMesh();

    // find inside and outside positions of nodes
    m.FindNodePositions();

    // find undecided nodes
    // * for serial simulations all node positions should be set
    // * for parallel simulations there can be some undecided nodes

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( dis_, m, *mesh_ ) );

    if(communicate) cut_parallel->CommunicateNodePositions();

    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes
    #ifdef DOFSETS_NEW
      mesh_->CreateNodalDofSetNEW( include_inner, dis_);
    #else
      m.FindNodalDOFSets( include_inner );
    #endif

    if(communicate) cut_parallel->CommunicateNodeDofSetNumbers();

  }

}


/*------------------------------------------------------------------------------------------------*
 *cut routine for standard non-parallel framework (only for cuttest)                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::Cut( bool include_inner, std::string VCellgausstype, std::string BCellgausstype )
{
  mesh_->Cut( include_inner, VCellgausstype, BCellgausstype );
}


void GEO::CutWizard::PrintCellStats()
{
  mesh_->PrintCellStats();
}

void GEO::CutWizard::DumpGmshNumDOFSets( bool include_inner)
{
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << filename;

  mesh_->DumpGmshNumDOFSets( str.str(), include_inner, dis_ );
}

void GEO::CutWizard::DumpGmshVolumeCells( bool include_inner )
{
  std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << name
      << ".CUT_volumecells."
      << dis_.Comm().MyPID()
      << ".pos";
  mesh_->DumpGmshVolumeCells( str.str(), include_inner );
}

void GEO::CutWizard::DumpGmshIntegrationCells()
{
  std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << name
      << ".CUT_integrationcells."
      << dis_.Comm().MyPID()
      << ".pos";
  mesh_->DumpGmshIntegrationCells( str.str() );
}

#if 0
Teuchos::RCP<Epetra_CrsGraph> GEO::CutWizard::MatrixGraph( const CutDofSet & dofset, const Epetra_Map & dbcmap )
{
  const Epetra_Map & noderowmap = * dis_.NodeRowMap();
  const Epetra_Map & dofrowmap  = * dofset.DofRowMap();

  int numnode = noderowmap.NumMyElements();

  // build graph

  std::map<int, std::set<int> > graph;

  int row = 0;
  for ( int i=0; i<numnode; ++i )
  {
    DRT::Node * rn = dis_.lRowNode( i );

    int numelements = rn->NumElement();
    DRT::Element** elements = rn->Elements();
    if ( elements==NULL )
      dserror( "no elements at node %d", rn->Id() );

    GEO::CUT::Node * crn = mesh_->GetNode( rn->Id() );
    if ( crn!=NULL )
    {
      std::vector<int> rdofs = dis_.Dof( rn );
      int rplain_numdf = dofset.PlainNumDofPerNode( *rn );
      for ( int k=0; k<numelements; ++k )
      {
        DRT::Element * e = elements[k];
        int numnodes = e->NumNode();
        DRT::Node ** nodes = e->Nodes();

        GEO::CUT::ElementHandle * ce = mesh_->GetElement( e->Id() );
        if ( ce!=NULL )
        {
          std::set<GEO::CUT::VolumeCell*> cells;
          ce->GetVolumeCells( cells );
          for ( std::set<GEO::CUT::VolumeCell*>::iterator i=cells.begin();
                i!=cells.end();
                ++i )
          {
            GEO::CUT::VolumeCell * c = *i;
            const std::vector<int> & nodaldofset = c->NodalDofSet();
            if ( nodaldofset.size()!=numnodes )
            {
              throw std::runtime_error( "number of nodes mismatch" );
            }
            for ( int l=0; l<numnodes; ++l )
            {
              DRT::Node * cn = nodes[l];
              std::vector<int> cdofs;
              dofset.Dof( *cn, nodaldofset[l], cdofs );
            }
          }
        }
        else
        {
          for ( int l=0; l<numnodes; ++l )
          {
            DRT::Node * cn = nodes[l];

            std::vector<int> cdofs = dis_.Dof( cn );
          }
        }
      }
    }
    else
    {
      std::vector<int> rdofs = dis_.Dof( rn );
      for ( unsigned rj=0; rj<rdofs.size(); ++rj, ++row )
      {
        std::set<int> & rowset = graph[rdofs[rj]];
        if ( not dbcmap.MyGID( rdofs[rj] ) )
        {
          // non-Dirichlet row
          for ( int k=0; k<numelements; ++k )
          {
            DRT::Element * e = elements[k];
            int numnodes = e->NumNode();
            DRT::Node ** nodes = e->Nodes();
            for ( int l=0; l<numnodes; ++l )
            {
              DRT::Node * cn = nodes[l];

              std::vector<int> cdofs = dis_.Dof( cn );
              rowset.insert( cdofs.begin(), cdofs.end() );
            }
          }
        }
        else
        {
          // just diagonal entry on Dirichlet rows
          rowset.insert( rdofs[rj] );
        }
      }
    }
  }

  // setup graph with row length and column indices per row

  std::vector<int> sizes;
  sizes.reserve( graph.size() );
  for ( std::map<int, std::set<int> >::iterator i=graph.begin(); i!=graph.end(); ++i )
  {
    std::set<int> & rowset = i->second;
    unsigned s = rowset.size();
    sizes.push_back( s );
  }

  Teuchos::RCP<Epetra_CrsGraph> crsgraph =
    Teuchos::rcp( new Epetra_CrsGraph( Copy, dofrowmap, &sizes[0], true ) );

  for ( std::map<int, std::set<int> >::iterator i=graph.begin(); i!=graph.end(); ++i )
  {
    int gid = i->first;
    std::set<int> & rowset = i->second;
    unsigned s = rowset.size();
    std::vector<int> row;
    row.reserve( s );
    row.assign( rowset.begin(), rowset.end() );

    int err = crsgraph->InsertGlobalIndices( gid, row.size(), &row[0] );
    if ( err )
      dserror( "InsertGlobalIndices failed: %d", err );
  }

  crsgraph->FillComplete();

  return crsgraph;
  return Teuchos::null;
}
#endif

