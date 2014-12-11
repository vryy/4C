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


#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_side.H"
#include "../drt_cut/cut_parallel.H"

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


/*------------------------------------------------------------------------------------------------*
 * cut routine for parallel framework in XFSI and XFLUIDFLUID                        schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::Run_Cut(
    bool include_inner,
    INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
    INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
    bool tetcellsonly,
    bool screenoutput)
{

  mesh_->Status();

  // just for time measurement
  dis_.Comm().Barrier();

  //----------------------------------------------------------
  // Selfcut (2/6 Cut_SelfCut)
  {
    const double t_start = Teuchos::Time::wallTime();

    // cut the mesh
    mesh_->Cut_SelfCut(include_inner, screenoutput);

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime() - t_start;
    if (myrank_ == 0)
      IO::cout << "\t\t\t\t... Success (" << t_diff << " secs)" << IO::endl;
  }
  //----------------------------------------------------------
  // Cut Part I: Collision Detection (3/6 Cut_CollisionDetection)
  {
    const double t_start = Teuchos::Time::wallTime();

    // cut the mesh
    mesh_->Cut_CollisionDetection(include_inner, screenoutput);

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime() - t_start;
    if (myrank_ == 0)
      IO::cout << "\t\t... Success (" << t_diff << " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // Cut Part II: Intersection (4/6 Cut_Intersection)
  {
    const double t_start = Teuchos::Time::wallTime();

    // call cut for cutting with a mesh
    mesh_->Cut_MeshIntersection( screenoutput );
    // call cut for cutting with a levelset field
    mesh_->Cut_Mesh( screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << "\t\t\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // Cut Part III & IV: Element Selection and DOF-Set Management (5/6 Cut_Positions_Dofsets)
  {
    const double t_start = Teuchos::Time::wallTime();

    FindPositionDofSets( include_inner, screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << "\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // Cut Part V & VI: Polyhedra Integration and Boundary Tessellation (6/6 Cut_Finalize)
  {
    const double t_start = Teuchos::Time::wallTime();

    // perform tessellation or moment fitting on the mesh
    mesh_->Cut_Finalize( include_inner, VCellgausstype, BCellgausstype, tetcellsonly, screenoutput );

    // just for time measurement
    dis_.Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 ) IO::cout << "\t\t\t\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  mesh_->Status(VCellgausstype);
}


/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::FindPositionDofSets(bool include_inner, bool screenoutput)
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
    FindNodePositions();


    //TWO_PHASE_XFEM_FIX: Does this apply for level set cut as well?
    // find undecided nodes
    // * for serial simulations all node positions should be set
    // * for parallel simulations there can be some undecided nodes

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( dis_, m, *mesh_ ) );

    cut_parallel->CommunicateNodePositions();


    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes
    mesh_->CreateNodalDofSet( include_inner, dis_);

    cut_parallel->CommunicateNodeDofSetNumbers(include_inner);

  }

}

/*------------------------------------------------------------------------------------------------*
 * Print the number of volumecells and boundarycells generated over the whole mesh during the cut *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::PrintCellStats()
{
  mesh_->PrintCellStats();
}


/*------------------------------------------------------------------------------------------------*
 * Write the DOF details of the nodes                                                             *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::DumpGmshNumDOFSets( bool include_inner)
{
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << filename;

  mesh_->DumpGmshNumDOFSets( str.str(), include_inner, dis_ );
}


/*------------------------------------------------------------------------------------------------*
 * Write volumecell output in GMSH format throughout the domain                                   *
 *------------------------------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------------------------------*
 * Write the integrationcells and boundarycells in GMSH format throughout the domain              *
 *------------------------------------------------------------------------------------------------*/
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


