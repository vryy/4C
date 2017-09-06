/*!-----------------------------------------------------------------------------------------------*
\file cut_cutwizard.cpp

\brief class that provides the common functionality for a mesh cut based on a level set field or on surface meshes

\level 3
<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "cut_cutwizard.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_cut/cut_combintersection.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_parallel.H"
#include "../drt_cut/cut_sidehandle.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CutWizard::BackMesh::Init(
    const Teuchos::RCP<const Epetra_Vector> & back_disp_col,
    const Teuchos::RCP<const Epetra_Vector> & back_levelset_col)
{
  back_disp_col_ = back_disp_col;
  back_levelset_col_ = back_levelset_col;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CutWizard::BackMesh::NumMyColElements() const
{
  return back_discret_->NumMyColElements();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Element* GEO::CutWizard::BackMesh::lColElement( int lid ) const
{
  return back_discret_->lColElement( lid );
}

/*-------------------------------------------------------------*
 * constructor
 *-------------------------------------------------------------*/
GEO::CutWizard::CutWizard( const Teuchos::RCP<DRT::Discretization> & backdis )
    : back_mesh_( Teuchos::rcp( new CutWizard::BackMesh( backdis, this ) ) ),
      comm_( backdis->Comm() ),
      myrank_ ( backdis->Comm().MyPID() ),
      intersection_(Teuchos::rcp( new GEO::CUT::CombIntersection( myrank_ ))),
      do_mesh_intersection_(false),
      do_levelset_intersection_(false),
      level_set_sid_(-1),
      VCellgausstype_( INPAR::CUT::VCellGaussPts_Tessellation ),
      BCellgausstype_( INPAR::CUT::BCellGaussPts_Tessellation ),
      gmsh_output_( false ),
      tetcellsonly_( false ),
      screenoutput_( false ),
      lsv_only_plus_domain_(true),
      is_set_options_(false)
{

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CutWizard::CutWizard( const Epetra_Comm & comm )
    : back_mesh_( Teuchos::null ),
      comm_( comm ),
      myrank_ ( comm.MyPID() ),
      intersection_( Teuchos::rcp( new GEO::CUT::CombIntersection( myrank_ ) ) ),
      do_mesh_intersection_(false),
      do_levelset_intersection_(false),
      level_set_sid_(-1),
      VCellgausstype_( INPAR::CUT::VCellGaussPts_Tessellation ),
      BCellgausstype_( INPAR::CUT::BCellGaussPts_Tessellation ),
      gmsh_output_( false ),
      tetcellsonly_( false ),
      screenoutput_( false ),
      lsv_only_plus_domain_( false ),
      is_set_options_( false ),
      is_cut_prepare_performed_(false)
{

}

/*========================================================================*/
//! @name Setters
/*========================================================================*/

/*-------------------------------------------------------------*
 * set options and flags used during the cut
*--------------------------------------------------------------*/
void GEO::CutWizard::SetOptions(
    INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy, //!< strategy for nodal dofset management
    INPAR::CUT::VCellGaussPts VCellgausstype,   //!< Gauss point generation method for Volumecell
    INPAR::CUT::BCellGaussPts BCellgausstype,   //!< Gauss point generation method for Boundarycell
    bool gmsh_output,                           //!< print write gmsh output for cut
    bool positions,                             //!< set inside and outside point, facet and volumecell positions
    bool tetcellsonly,                          //!< generate only tet cells
    bool screenoutput                           //!< print screen output
    )
{

  VCellgausstype_ = VCellgausstype;
  BCellgausstype_ = BCellgausstype;
  gmsh_output_    = gmsh_output;
  tetcellsonly_   = tetcellsonly;
  screenoutput_   = screenoutput;

  // set position option to the intersection class
  intersection_->SetFindPositions( positions );
  intersection_->SetNodalDofSetStrategy( nodal_dofset_strategy );

  is_set_options_ = true;
}


/*-------------------------------------------------------------*
 * set displacement and level-set vectors used during the cut
*--------------------------------------------------------------*/
void GEO::CutWizard::SetBackgroundState(
    Teuchos::RCP<const Epetra_Vector> back_disp_col,      //!< col vector holding background ALE displacements for backdis
    Teuchos::RCP<const Epetra_Vector> back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
    int level_set_sid                                     //!< global id for level-set side
)
{
  // set state vectors used in cut
  back_mesh_->Init( back_disp_col, back_levelset_col );
  level_set_sid_ = level_set_sid;

  do_levelset_intersection_ = back_mesh_->IsLevelSet();

}


/*-------------------------------------------------------------*
 * set displacement and level-set vectors used during the cut
*--------------------------------------------------------------*/
void GEO::CutWizard::AddCutterState(
    const int mc_idx,
    Teuchos::RCP<DRT::Discretization> cutter_dis,
    Teuchos::RCP<const Epetra_Vector> cutter_disp_col
)
{
  AddCutterState(0,cutter_dis, cutter_disp_col, 0);
}

/*-------------------------------------------------------------*
 * set displacement and level-set vectors used during the cut
*--------------------------------------------------------------*/
void GEO::CutWizard::AddCutterState(
    const int mc_idx,
    Teuchos::RCP<DRT::Discretization> cutter_dis,
    Teuchos::RCP<const Epetra_Vector> cutter_disp_col,
    const int                         start_ele_gid
)
{
  std::map<int, Teuchos::RCP<CutterMesh> >::iterator cm = cutter_meshes_.find(mc_idx);

  if(cm!=cutter_meshes_.end()) dserror("cutter mesh with mesh coupling index %i already set", mc_idx);

  cutter_meshes_[mc_idx] = Teuchos::rcp(
      new CutterMesh(
      cutter_dis,
      cutter_disp_col,
      start_ele_gid)
  );

  do_mesh_intersection_ = true;

}

/*-------------------------------------------------------------*
 * Mark surfaces loaded into cut with background surfaces
*--------------------------------------------------------------*/
void GEO::CutWizard::SetMarkedConditionSides(
    //const int mc_idx,                                       //Not needed (for now?)
    Teuchos::RCP<DRT::Discretization> cutter_dis,
    //Teuchos::RCP<const Epetra_Vector> cutter_disp_col,      //Not needed (for now?)
    const int                         start_ele_gid
)
{
  //Set the counter to the gid.
  //  -- Set ids in correspondence to this ID.
  //  -- Loop over the surface elements and find (if it exists) a corresponding side loaded into the cut
  //  ## WARNING: Not sure what happens if it doesn't find a surface?
  for( int lid=0; lid < cutter_dis->NumMyRowElements(); ++lid)
  {
    DRT::Element* cutter_dis_ele = cutter_dis->lRowElement(lid);

    const int numnode = cutter_dis_ele->NumNode();
    const int * nodeids = cutter_dis_ele->NodeIds();
    std::vector<int> node_ids_of_cutterele( nodeids, nodeids+numnode );

    const int eid = cutter_dis_ele->Id();       // id of marked side based on the cutter discretization
    const int marked_sid = eid + start_ele_gid; // id of marked side within the cut library

    //Get sidehandle to corresponding background surface discretization
    // -- if it exists!!!
    GEO::CUT::SideHandle* cut_sidehandle = intersection_->GetMeshHandle().GetSide(node_ids_of_cutterele);

    if(cut_sidehandle != NULL)
    {
      GEO::CUT::plain_side_set cut_sides;
      cut_sidehandle->CollectSides(cut_sides);

      // Set Id's and mark the sides in correspondence with the coupling manager object.
      for (GEO::CUT::plain_side_set::iterator it=cut_sides.begin(); it!=cut_sides.end(); ++it)
      {
        (*it)->SetMarkedSideProperties(marked_sid,GEO::CUT::mark_and_create_boundarycells);
      }
    }
    else
      dserror("If we don't find a marked side it's not sure what happens... You are on your own!");
  }
}

/*========================================================================*/
//! @name main Cut call
/*========================================================================*/
/*-------------------------------------------------------------*
* main Cut call
*--------------------------------------------------------------*/
void GEO::CutWizard::Cut(
    bool include_inner //!< perform cut in the interior of the cutting mesh
)
{

  // safety checks if the cut is initialized correctly
  if(!SafetyChecks(false)) return;

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CutWizard::Cut" );

  if ( myrank_ == 0 and screenoutput_ )
    IO::cout << "\nGEO::CutWizard::Cut:" << IO::endl;

  const double t_start = Teuchos::Time::wallTime();

  /* wirtz 08/14:
   * preprocessing: everything above should only be done once in a simulation;
   *                so it should be moved before the time loop into a preprocessing
   *                step
   * runtime:       everything below should be done in every Newton increment */

  //--------------------------------------
  // perform the actual cut, the intersection
  //--------------------------------------
  Run_Cut( include_inner);


  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( myrank_ == 0  and screenoutput_)
  {
    IO::cout << "\n\t\t\t\t\t\t\t... Success (" << t_end  <<  " secs)\n" << IO::endl;
  }

  //--------------------------------------
  // write statistics and output to screen and files
  //--------------------------------------
  Output(include_inner);

}

/*-------------------------------------------------------------*
* prepare the cut, add background elements and cutting sides
*--------------------------------------------------------------*/
void GEO::CutWizard::Prepare()
{

  // safety checks if the cut is initialized correctly
  if(!SafetyChecks(true)) return;

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/6 --- Cut_Initialize" );

  const double t_start = Teuchos::Time::wallTime();

  if ( myrank_ == 0 and screenoutput_ )
    IO::cout << "\nGEO::CutWizard::Prepare:" << IO::endl;

  if( myrank_ == 0 and screenoutput_ )
    IO::cout << "\n\t * 1/6 Cut_Initialize ...";

  // fill the cutwizard cw with information:
  // build up the mesh_ (normal background mesh) and the cut_mesh_ (cutter mesh) created by the meshhandle:
  // REMARK: DO NOT CHANGE THE ORDER of 1. and 2.
  // 1. Add CutSides (sides of the cutterdiscretization)
  //      -> Update the current position of all cutter-nodes dependent on displacement idispcol
  // 2. Add Elements (elements of the background discretization)

  // Ordering is very important because first we add all cut sides, and create a bounding box which contains
  // all the cut sides
  // Then, when adding elements from background discret, only the elements that intersect this
  // bounding box are added
  // Changing the order would render in problems when all bg-elements on one proc are within the structure,
  // then the bb around the bg-mesh on this proc has no intersection with an an bb around an side element


  // 1. Add CutSides (possible sides of the cutter-discretization and a possible level-set side)
  AddCuttingSides();

  // 2. Add background elements dependent on bounding box created by the CutSides in 1.
  AddBackgroundElements();


  // build the static search tree for the collision detection in the self cut
  intersection_->BuildSelfCutTree();

  // build the static search tree for the collision detection
  intersection_->BuildStaticSearchTree();

  const double t_mid = Teuchos::Time::wallTime()-t_start;
  if ( myrank_ == 0  and screenoutput_)
  {
    IO::cout << "\t\t\t... Success (" << t_mid  <<  " secs)" << IO::endl;
  }

  is_cut_prepare_performed_=true;

}

/*-------------------------------------------------------------*
* add all cutting sides (mesh and level-set sides)
*--------------------------------------------------------------*/
void GEO::CutWizard::AddCuttingSides()
{
  // add all mesh cutting sides
  if(do_mesh_intersection_) AddMeshCuttingSide();

  // add a new level-set side
  if(do_levelset_intersection_) AddLSCuttingSide();
}

/*-------------------------------------------------------------*
* add level-set cutting side
*--------------------------------------------------------------*/
void GEO::CutWizard::AddLSCuttingSide()
{
  // add a new level-set side
  intersection_->AddLevelSetSide(level_set_sid_);
}


/*-------------------------------------------------------------*
* add all mesh-cutting sides of all cutting discretizations
*--------------------------------------------------------------*/
void GEO::CutWizard::AddMeshCuttingSide()
{

  // loop all mesh coupling objects
  for(std::map<int, Teuchos::RCP<CutterMesh> >::iterator it=cutter_meshes_.begin();
      it != cutter_meshes_.end();
      it++)
  {
    Teuchos::RCP<CutterMesh> cutter_mesh= it->second;

    AddMeshCuttingSide(
        cutter_mesh->cutterdis_,
        cutter_mesh->cutter_disp_col_,
        cutter_mesh->start_ele_gid_);
  }
}



/*-------------------------------------------------------------*
* add all cutting sides from the cut-discretization
*--------------------------------------------------------------*/
void GEO::CutWizard::AddMeshCuttingSide(
    Teuchos::RCP<DRT::Discretization> cutterdis,
    Teuchos::RCP<const Epetra_Vector> cutter_disp_col,
    const int start_ele_gid        ///< mesh coupling index
    )
{
  if(cutterdis == Teuchos::null) dserror("cannot add mesh cutting sides for invalid cutter discretiaztion!");

  std::vector<int> lm;
  std::vector<double> mydisp;

  int numcutelements = cutterdis->NumMyColElements();


  for ( int lid = 0; lid < numcutelements; ++lid )
  {
    DRT::Element * element = cutterdis->lColElement(lid);

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      lm.clear();
      mydisp.clear();
      cutterdis->Dof(&node, lm);

      LINALG::Matrix<3, 1> x( node.X() );

      if(cutter_disp_col != Teuchos::null)
      {
        if(lm.size() == 3) // case for BELE3 boundary elements
        {
          DRT::UTILS::ExtractMyValues(*cutter_disp_col,mydisp,lm);
        }
        else if(lm.size() == 4) // case for BELE3_4 boundary elements
        {
          // copy the first three entries for the displacement, the fourth entry should be zero if BELE3_4 is used for cutdis instead of BELE3
          std::vector<int> lm_red; // reduced local map
          lm_red.clear();
          for(int k=0; k< 3; k++) lm_red.push_back(lm[k]);

          DRT::UTILS::ExtractMyValues(*cutter_disp_col,mydisp,lm_red);
        }
        else dserror("wrong number of dofs for node %i", lm.size());

        if (mydisp.size() != 3)
          dserror("we need 3 displacements here");

        LINALG::Matrix<3, 1> disp( &mydisp[0], true );

        //update x-position of cutter node for current time step (update with displacement)
        x.Update( 1, disp, 1 );
      }
      std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
    }

    // add the side of the cutter-discretization
    AddMeshCuttingSide( 0, element, xyze, start_ele_gid );
  }
}

/*-------------------------------------------------------------*
* prepare the cut, add background elements and cutting sides
*--------------------------------------------------------------*/
void GEO::CutWizard::AddMeshCuttingSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze, const int start_ele_gid )
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );

  const int eid = ele->Id();       // id of cutting side based on the cutter discretization
  const int sid = eid + start_ele_gid; // id of cutting side within the cut library

  intersection_->AddMeshCuttingSide( sid, nids, xyze, ele->Shape(), mi );
}

/*-------------------------------------------------------------*
 * add elements from the background discretization
 *-------------------------------------------------------------*/
void GEO::CutWizard::AddBackgroundElements()
{

  // vector with nodal level-set values
  std::vector<double> myphinp;

  // Loop over all Elements to find cut elements and add them to the LevelsetIntersection class
  // Brute force method.
  int numelements = back_mesh_->NumMyColElements();

  for ( int lid = 0; lid < numelements; ++lid )
  {
    const DRT::Element * element = back_mesh_->lColElement( lid );

    LINALG::SerialDenseMatrix xyze;

    GetPhysicalNodalCoordinates( element, xyze );

    if( back_mesh_->IsLevelSet() )
    {
      myphinp.clear();

      DRT::UTILS::ExtractMyNodeBasedValues( element, myphinp, back_mesh_->BackLevelSetCol() );
      AddElement( element, xyze, &myphinp[0], lsv_only_plus_domain_ );
    }
    else
    {
      AddElement( element, xyze, NULL );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CutWizard::GetPhysicalNodalCoordinates(
    const DRT::Element * element, LINALG::SerialDenseMatrix & xyze ) const
{
  std::vector<int> lm;
  std::vector<double> mydisp;

  const int numnode = element->NumNode();
  const DRT::Node * const * nodes = element->Nodes();

  xyze.Shape( 3, numnode );
  for ( int i=0; i < numnode; ++i )
  {
    const DRT::Node & node = *nodes[i];

    LINALG::Matrix<3, 1> x( node.X() );

    if( back_mesh_->IsBackDisp() )
    {
      // castt to DiscretizationXFEM
      Teuchos::RCP<DRT::DiscretizationXFEM>  xbackdis =
          Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(back_mesh_->GetPtr(), true);

      lm.clear();
      mydisp.clear();

      xbackdis->InitialDof(&node, lm); //to get all dofs of background (also not active ones at the moment!)

      if(lm.size() == 3) // case used actually?
      {
        DRT::UTILS::ExtractMyValues( back_mesh_->BackDispCol(), mydisp, lm );
      }
      else if(lm.size() == 4) // case xFluid ... just take the first three
      {
        // copy the first three entries for the displacement, the fourth entry an all others
        std::vector<int> lm_red; // reduced local map
        lm_red.clear();
        for(int k=0; k< 3; k++) lm_red.push_back(lm[k]);

        DRT::UTILS::ExtractMyValues( back_mesh_->BackDispCol(), mydisp, lm_red );
      }
      else
        dserror("wrong number of dofs for node %i", lm.size());

      if (mydisp.size() != 3)
        dserror("we need 3 displacements here");

      LINALG::Matrix<3, 1> disp( &mydisp[0], true );

      //update x-position of cutter node for current time step (update with displacement)
      x.Update( 1, disp, 1 );
    }
    std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
  }
}


/*-------------------------------------------------------------*
* Add this background mesh element to the intersection class
*--------------------------------------------------------------*/
void GEO::CutWizard::AddElement(
    const DRT::Element * ele,
    const Epetra_SerialDenseMatrix & xyze,
    double* myphinp,
    bool lsv_only_plus_domain )
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );

  //If include_inner == false then add elements with negative level set values to discretization.
  intersection_->AddElement( ele->Id(), nids, xyze, ele->Shape(), myphinp, lsv_only_plus_domain );
}

/*------------------------------------------------------------------------------------------------*
 * perform the actual cut, the intersection
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::Run_Cut(
    bool include_inner  //!< perform cut in the interior of the cutting mesh
    )
{

  intersection_->Status();

  // just for time measurement
  comm_.Barrier();

  if(do_mesh_intersection_)
  {
    //----------------------------------------------------------
    // Selfcut (2/6 Cut_SelfCut)
    {
      const double t_start = Teuchos::Time::wallTime();

      // cut the mesh
      intersection_->Cut_SelfCut(include_inner, screenoutput_);

      // just for time measurement
      comm_.Barrier();

      const double t_diff = Teuchos::Time::wallTime() - t_start;
      if (myrank_ == 0 and screenoutput_)
        IO::cout << "\t\t\t\t... Success (" << t_diff << " secs)" << IO::endl;
    }
    //----------------------------------------------------------
    // Cut Part I: Collision Detection (3/6 Cut_CollisionDetection)
    {
      const double t_start = Teuchos::Time::wallTime();

      // cut the mesh
      intersection_->Cut_CollisionDetection(include_inner, screenoutput_);

      // just for time measurement
      comm_.Barrier();

      const double t_diff = Teuchos::Time::wallTime() - t_start;
      if (myrank_ == 0 and screenoutput_)
        IO::cout << "\t\t... Success (" << t_diff << " secs)" << IO::endl;
    }
  }

  //----------------------------------------------------------
  // Cut Part II: Intersection (4/6 Cut_Intersection)
  {
    const double t_start = Teuchos::Time::wallTime();

    intersection_->Cut(screenoutput_);

    // just for time measurement
    comm_.Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 and screenoutput_ ) IO::cout << "\t\t\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // Cut Part III & IV: Element Selection and DOF-Set Management (5/6 Cut_Positions_Dofsets)
  {
    const double t_start = Teuchos::Time::wallTime();

    FindPositionDofSets( include_inner );

    // just for time measurement
    comm_.Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 and screenoutput_ ) IO::cout << "\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // Cut Part V & VI: Polyhedra Integration and Boundary Tessellation (6/6 Cut_Finalize)
  {
    const double t_start = Teuchos::Time::wallTime();

    // perform tessellation or moment fitting on the mesh
    intersection_->Cut_Finalize( include_inner, VCellgausstype_, BCellgausstype_, tetcellsonly_, screenoutput_ );

    // just for time measurement
    comm_.Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 and screenoutput_ ) IO::cout << "\t\t\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  intersection_->Status(VCellgausstype_);

  Post_Run_Cut( include_inner );
}


/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::FindPositionDofSets(bool include_inner)
{

  comm_.Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets (parallel)" );

  if(myrank_==0 and screenoutput_) IO::cout << "\t * 5/6 Cut_Positions_Dofsets (parallel) ...";

//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  GEO::CUT::Options options;
  intersection_->GetOptions( options );

  if ( options.FindPositions() )
  {

    GEO::CUT::Mesh & m = intersection_->NormalMesh();

    bool communicate = (comm_.NumProc() > 1);

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::null;

    if ( communicate )
    {
      cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( back_mesh_->GetPtr(), m, *intersection_ ) );
    }

    // find inside and outside positions of nodes
    // first for mesh cut and distribute data in parallel, after that do the same for the level-set cut

    //--------------------------------------------
    // first, set the position for the mesh cut
    if(do_mesh_intersection_)
    {
      m.FindNodePositions();

      if ( communicate )
        cut_parallel->CommunicateNodePositions();

    }

    //--------------------------------------------
    // second, set the position for the level-set cut (no parallel communication necessary)
    if(do_levelset_intersection_)
    {
      m.FindLSNodePositions();
    }

    if(do_mesh_intersection_)
    {
      m.FindFacetPositions();
    }

    //--------------------------------------------
    comm_.Barrier();

    // find number and connection of dofsets at nodes from cut volumes
    intersection_->CreateNodalDofSet( include_inner, back_mesh_->Get() );

    if (communicate) cut_parallel->CommunicateNodeDofSetNumbers(include_inner);

  }
}


bool GEO::CutWizard::SafetyChecks(bool is_prepare_cut_call)
{

  if( !is_set_options_ )
    dserror( "You have to call SetOptions() before you can use the CutWizard" );

  if(!is_prepare_cut_call and !is_cut_prepare_performed_)
    dserror( "You have to call PrepareCut() before you can call the Cut-routine" );

  if( !do_mesh_intersection_ and !do_levelset_intersection_ )
  {
    if ( myrank_ == 0  and is_prepare_cut_call)
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n";
      std::cout << "WARNING: No mesh intersection and no level-set intersection! \n" <<
                   "         Why do you call the CUT-library?\n";
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n";
    }
    return false;
  }

  return true;

}

/*------------------------------------------------------------------------------------------------*
 * write statistics and output to screen and files
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::Output(bool include_inner)
{
  if(gmsh_output_) DumpGmshNumDOFSets(include_inner);

#ifdef DEBUG
  PrintCellStats();
#endif

  if(gmsh_output_)
  {
    DumpGmshIntegrationCells();
    DumpGmshVolumeCells( include_inner );
  }
}


/*------------------------------------------------------------------------------------------------*
 * Print the number of volumecells and boundarycells generated over the whole mesh during the cut *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::PrintCellStats()
{
  intersection_->PrintCellStats();
}


/*------------------------------------------------------------------------------------------------*
 * Write the DOF details of the nodes                                                             *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::DumpGmshNumDOFSets( bool include_inner)
{
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << filename;

  intersection_->DumpGmshNumDOFSets( str.str(), include_inner, back_mesh_->Get() );
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
      << myrank_
      << ".pos";
  intersection_->DumpGmshVolumeCells( str.str(), include_inner );
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
      << myrank_
      << ".pos";
  intersection_->DumpGmshIntegrationCells( str.str() );
}


/*========================================================================*/
//! @name Getters
/*========================================================================*/

GEO::CUT::SideHandle * GEO::CutWizard::GetSide( std::vector<int>& nodeids )
{
  return intersection_->GetSide( nodeids );
}

GEO::CUT::SideHandle * GEO::CutWizard::GetSide( int sid )
{
  return intersection_->GetSide( sid );
}

GEO::CUT::ElementHandle * GEO::CutWizard::GetElement( const DRT::Element * ele ) const
{
  return intersection_->GetElement( ele->Id() );
}

GEO::CUT::Node * GEO::CutWizard::GetNode( int nid )
{
  return intersection_->GetNode( nid );
}

GEO::CUT::SideHandle * GEO::CutWizard::GetMeshCuttingSide( int sid, int mi )
{
  return intersection_->GetCutSide(sid, mi);
}

bool GEO::CutWizard::HasLSCuttingSide( int sid )
{
  return intersection_->HasLSCuttingSide(sid);
}





