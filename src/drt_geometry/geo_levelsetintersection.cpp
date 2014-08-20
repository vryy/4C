/*!-----------------------------------------------------------------------------------------------*
\file geo_levelsetintersection.cpp

\brief class that providesthe wizard for the functionality for a mesh cut based on a level set field

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_levelsetintersection.H"
#include "../drt_cut/cut_parallel.H"

#include "geo_levelsetintersection.H"

/// Constructor
GEO::CutWizardLevelSet::CutWizardLevelSet( const DRT::Discretization & dis )
:   CutWizard( dis )
{
  levelsetintersection_ = Teuchos::rcp( new GEO::CUT::LevelSetIntersection( myrank_ ) );
  mesh_ = levelsetintersection_;

}

/*------------------------------------------------------------------------------------------------*
 * check if element is cut and add element to cut libraries if this is the case.     winter 08/14 *
 * myphinp needs to on the node map for the element. This has to be done before the call.         *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardLevelSet::AddElement(DRT::Element * ele, std::vector<double> myphinp)
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
  levelsetintersection_->AddElement( ele->Id(), nids, xyze, &myphinp[0], ele->Shape() );
}

/*------------------------------------------------------------------------------------------------*
 * cut routine for parallel framework of level set cut                               schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardLevelSet::CutParallel( bool include_inner,
                                  INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
                                  INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
                                  bool tetcellsonly,
                                  bool screenoutput)
{
  // for level set cut we have to communicate dofset data
  bool communicate = true;

  if(dis_.Comm().NumProc() == 1) communicate = false;

  levelsetintersection_->Status();

  // just for time measurement
  dis_.Comm().Barrier();

  //----------------------------------------------------------
  // FIRST step (1/3): cut the mesh
  {
    const double t_start = Teuchos::Time::wallTime();

    // cut the mesh
    levelsetintersection_->Cut_Mesh( include_inner,screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // SECOND step (2/3): find node positions and create dofset in PARALLEL
  {
    const double t_start = Teuchos::Time::wallTime();

    CutParallel_FindPositionDofSets( include_inner, communicate ,screenoutput );

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
    levelsetintersection_->Cut_Finalize( include_inner, VCellgausstype, BCellgausstype, true, tetcellsonly, screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  levelsetintersection_->Status(VCellgausstype);
}

/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardLevelSet::CutParallel_FindPositionDofSets(bool include_inner, bool communicate, bool screenoutput)
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
    m.FindLSNodePositions(); // Prevents this function from being in CutWizard.


    //TWO_PHASE_XFEM_FIX_QB: Does this apply for level set cut as well?
    // find undecided nodes
    // * for serial simulations all node positions should be set
    // * for parallel simulations there can be some undecided nodes

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( dis_, m, *levelsetintersection_ ) );

    if(communicate) cut_parallel->CommunicateNodePositions();


    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes
    #ifdef DOFSETS_NEW
      mesh_->CreateNodalDofSetNEW( include_inner, dis_);
    #else
      m.FindNodalDOFSets( include_inner );
    #endif

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

}

/*------------------------------------------------------------------------------------------------*
 *cut routine for standard non-parallel framework (only for cuttest)                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizardLevelSet::Cut(
    bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
    INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
    bool screenoutput
)
{
  levelsetintersection_->Cut( include_inner, screenoutput );
}
