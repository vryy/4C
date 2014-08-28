/*!-----------------------------------------------------------------------------------------------*
\file geo_meshintersection.cpp

\brief class that provides to set up a mesh cut based on surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_parallel.H"

#include "geo_meshintersection.H"

/// Constructor
GEO::CutWizardMesh::CutWizardMesh( DRT::Discretization & dis, int numcutmesh )
:   CutWizard( dis )
{
  meshintersection_ = Teuchos::rcp( new GEO::CUT::MeshIntersection( numcutmesh, myrank_ ) );
  mesh_ = meshintersection_;

}

void GEO::CutWizardMesh::AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze )
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );
  meshintersection_->AddCutSide( ele->Id(), nids, xyze, ele->Shape(), mi );
}

void GEO::CutWizardMesh::AddElement( DRT::Element * ele )
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
  meshintersection_->AddElement( ele->Id(), nids, xyze, ele->Shape() );
}

GEO::CUT::SideHandle * GEO::CutWizardMesh::GetCutSide( int sid, int mi )
{
  return meshintersection_->GetCutSide( sid, mi );
}

/*------------------------------------------------------------------------------------------------*
 * cut routine for parallel framework in XFSI and XFLUIDFLUID                        schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::CutParallel( bool include_inner,
                                  INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
                                  INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
                                  bool tetcellsonly,
                                  bool screenoutput)
{
  // for XFSI and XFLUIDFLUID we have communicate node positions and dofset data
  bool communicate = true;

  if(dis_.Comm().NumProc() == 1) communicate = false;

  meshintersection_->Status();

  // just for time measurement
  dis_.Comm().Barrier();

  //----------------------------------------------------------
  // FIRST step (1/3): cut the mesh
  {
    const double t_start = Teuchos::Time::wallTime();

    // cut the mesh
    meshintersection_->Cut_Mesh( include_inner, screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // SECOND step (2/3): find node positions and create dofset in PARALLEL
  {
    const double t_start = Teuchos::Time::wallTime();

    CutParallel_FindPositionDofSets( include_inner, communicate, screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // THIRD step (3/3): perform tessellation or moment fitting on the mesh
  {
    const double t_start = Teuchos::Time::wallTime();

    // perform tessellation or moment fitting on the mesh
    meshintersection_->Cut_Finalize( include_inner, VCellgausstype, BCellgausstype, false, tetcellsonly, screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  meshintersection_->Status(VCellgausstype);
}


/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::CutParallel_FindPositionDofSets(bool include_inner, bool communicate, bool screenoutput)
{


  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 2/3 --- Cut_Positions_Dofsets (parallel)" );


  if(myrank_==0 and screenoutput) IO::cout << "\t * 2/3 Cut_Positions_Dofsets (parallel) ...";

//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  GEO::CUT::Options options;
  mesh_->GetOptions(options);

  if ( options.FindPositions() )
  {

    GEO::CUT::Mesh & m = mesh_->NormalMesh();

    // find inside and outside positions of nodes
    m.FindNodePositions(); // Prevents this function from being in CutWizard.

    // find undecided nodes
    // * for serial simulations all node positions should be set
    // * for parallel simulations there can be some undecided nodes

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( dis_, m, *meshintersection_ ) );

    if(communicate) cut_parallel->CommunicateNodePositions();

    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes
    mesh_->CreateNodalDofSetNEW( include_inner, dis_);

    //DumpGmshNumDOFSets(false);

    if(communicate) cut_parallel->CommunicateNodeDofSetNumbers(include_inner);


  }

//  // just for time measurement
//  dis_.Comm().Barrier();
//
//  //----------------------------------------------------------
//
//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0 )
//  {
//    IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//  }
//
//
//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0  and screenoutput)
//  {
//    IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//  }

}

/*------------------------------------------------------------------------------------------------*
 *cut routine for standard non-parallel framework (only for cuttest)                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardMesh::Cut(
    bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
    INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
    bool screenoutput
)
{
  meshintersection_->Cut( include_inner, VCellgausstype, BCellgausstype,screenoutput );
}
