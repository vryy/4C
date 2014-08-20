/*!-----------------------------------------------------------------------------------------------*
\file geo_intersection.cpp

\brief class that provides the common functionality for a mesh cut based on a level set field or on surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "../drt_fluid_xfluid/xfluid_defines.H"

#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_facet.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_side.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_parallel.H"

#include "geo_utils.H"

#include "geo_intersection.H"


GEO::CutWizard::CutWizard( const DRT::Discretization & dis )
  : dis_( dis ),
    mesh_(Teuchos::null)
{
  myrank_ = dis_.Comm().MyPID();
}

/*
 * Setters
 */

void GEO::CutWizard::SetFindPositions( bool positions )
{
  mesh_->SetFindPositions( positions );
}

/*
 * Getters
 */

GEO::CUT::SideHandle * GEO::CutWizard::GetSide( std::vector<int>& nodeids )
{
  return mesh_->GetSide( nodeids );
}

GEO::CUT::SideHandle * GEO::CutWizard::GetSide( int sid )
{
  return mesh_->GetSide( sid );
}

GEO::CUT::ElementHandle * GEO::CutWizard::GetElement( DRT::Element * ele )
{
  return mesh_->GetElement( ele->Id() );
}

GEO::CUT::Node * GEO::CutWizard::GetNode( int nid )
{
  return mesh_->GetNode( nid );
}

/*
 * Output
 */

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

