/*!-----------------------------------------------------------------------------------------------*
\file cut_cutwizard.cpp

\brief class that provides the common functionality for a mesh cut based on a level set field or on surface meshes

<pre>
\level 3
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
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


/*-------------------------------------------------------------*
 * constructor
*--------------------------------------------------------------*/
GEO::CutWizard::CutWizard( Teuchos::RCP<DRT::Discretization> dis )
  : backdis_( dis ),
    myrank_ ( backdis_->Comm().MyPID() ),
    intersection_(Teuchos::rcp( new GEO::CUT::CombIntersection( myrank_ ))),
    do_mesh_intersection_(false),
    do_levelset_intersection_(false),
    back_disp_col_(Teuchos::null),
    back_levelset_col_(Teuchos::null),
    level_set_sid_(-1),
    lsv_only_plus_domain_(true),
    is_set_options_(false)
{
  if(backdis_ == Teuchos::null) dserror("null pointer to background dis, invalid!");
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
    const int level_set_sid                               //!< global id for level-set side
)
{
  // set state vectors used in cut
  back_disp_col_     = back_disp_col;
  back_levelset_col_ = back_levelset_col;

  level_set_sid_ = level_set_sid;


  if(back_levelset_col_ != Teuchos::null)
    do_levelset_intersection_ = true;

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
/*
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- BuildDOFCellSets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::Facet::CornerPointsLocal" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets (parallel): distribute DofSetData" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::VCIntegrationRule" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::ListFacets-tmp1" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::ListFacets-tmp2" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::Facet::IsPlanar" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::EqnPlaneOfPolygon" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::ListFacets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::DeleteInlinePts" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::getAreaConvexQuad" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::getAreaTri" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::BoundaryCell::TransformLocalCoords" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::FacetIntegration::DivergenceIntegrationRule" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::VolumeCell::GenerateInternalGaussRule" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::Element::CreateSimpleShapedIntegrationCells" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::VolumeCell::GenerateBoundaryCells" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- creating map" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::LinearElementHandle::VolumeCellSets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- AssignNodalCellSet" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- GEO::CUT::Cmp::Compare(vc,vc)" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FindDOFSetsNEW" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- SortNodalDofSets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- CollectNodalDofSets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- CreateNodalDofSet" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- STEP 2" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FindNodalCellSets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- ConnectNodalDOFSets" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- GEO::CUT::Cmp::operator()" );
  }
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- GEO::CUT::NodalDofSet::NodalDofSet" );
  }
*/


  // safety checks if the cut is initialized correctly
  if(!is_set_options_) dserror("you have call SetOptions() before you can use the CutWizard");

  if(!do_mesh_intersection_ and !do_levelset_intersection_) dserror(" no mesh intersection and no level-set intersection! Why do you call the CUT-library?");


  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CutWizard::Cut" );

  if ( backdis_->Comm().MyPID() == 0 and screenoutput_)
    IO::cout << "\nGEO::CutWizard::Cut:" << IO::endl;

  const double t_start = Teuchos::Time::wallTime();

  //--------------------------------------
  // prepare the cut, add background elements and cutting sides
  //--------------------------------------
  Prepare();

  // wirtz 08/14:
  // preprocessing: everything above should only be done once in a simulation; so it should be moved before the time loop into a preprocessing step
  // runtime:       everything below should be done in every Newton increment

  //--------------------------------------
  // perform the actual cut, the intersection
  //--------------------------------------
  Run_Cut( include_inner);


  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( backdis_->Comm().MyPID() == 0  and screenoutput_)
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

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/6 --- Cut_Initialize" );

  const double t_start = Teuchos::Time::wallTime();

  if(backdis_->Comm().MyPID()==0 and screenoutput_)
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
  if ( backdis_->Comm().MyPID() == 0  and screenoutput_)
  {
    IO::cout << "\t\t\t... Success (" << t_mid  <<  " secs)" << IO::endl;
  }

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
*--------------------------------------------------------------*/
void GEO::CutWizard::AddBackgroundElements()
{

  std::vector<int> lm;
  std::vector<double> mydisp;

  // vector with nodal level-set values
  std::vector<double> myphinp;

  // Loop over all Elements to find cut elements and add them to the LevelsetIntersection class
  // Brute force method.
  int numelements = backdis_->NumMyColElements();



  for ( int lid = 0; lid < numelements; ++lid )
  {
    DRT::Element * element = backdis_->lColElement(lid);

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      LINALG::Matrix<3, 1> x( node.X() );

      if( back_disp_col_ != Teuchos::null)
      {
        // castt to DiscretizationXFEM
        Teuchos::RCP<DRT::DiscretizationXFEM>  xbackdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(backdis_, true);

        lm.clear();
        mydisp.clear();

        xbackdis->InitialDof(&node, lm); //to get all dofs of background (also not active ones at the moment!)

        if(lm.size() == 3) // case used actually?
        {
          DRT::UTILS::ExtractMyValues(*back_disp_col_,mydisp,lm);
        }
        else if(lm.size() == 4) // case xFluid ... just take the first three
        {
          // copy the first three entries for the displacement, the fourth entry an all others
          std::vector<int> lm_red; // reduced local map
          lm_red.clear();
          for(int k=0; k< 3; k++) lm_red.push_back(lm[k]);

          DRT::UTILS::ExtractMyValues(*back_disp_col_,mydisp,lm_red);
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

    if(back_levelset_col_ != Teuchos::null)
    {
      myphinp.clear();

      DRT::UTILS::ExtractMyNodeBasedValues(element, myphinp, *back_levelset_col_);
      AddElement(element, xyze, &myphinp[0], lsv_only_plus_domain_);
    }
    else
    {
      AddElement(element, xyze, NULL);
    }
  }
}

/*-------------------------------------------------------------*
* Add this background mesh element to the intersection class
*--------------------------------------------------------------*/
void GEO::CutWizard::AddElement( DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze, double* myphinp, bool lsv_only_plus_domain )
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
  backdis_->Comm().Barrier();

  if(do_mesh_intersection_)
  {
    //----------------------------------------------------------
    // Selfcut (2/6 Cut_SelfCut)
    {
      const double t_start = Teuchos::Time::wallTime();

      // cut the mesh
      intersection_->Cut_SelfCut(include_inner, screenoutput_);

      // just for time measurement
      backdis_->Comm().Barrier();

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
      backdis_->Comm().Barrier();

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
    backdis_->Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 and screenoutput_ ) IO::cout << "\t\t\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  //----------------------------------------------------------
  // Cut Part III & IV: Element Selection and DOF-Set Management (5/6 Cut_Positions_Dofsets)
  {
    const double t_start = Teuchos::Time::wallTime();

    FindPositionDofSets( include_inner );

    // just for time measurement
    backdis_->Comm().Barrier();

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
    backdis_->Comm().Barrier();

    const double t_diff = Teuchos::Time::wallTime()-t_start;
    if ( myrank_ == 0 and screenoutput_ ) IO::cout << "\t\t\t\t... Success (" << t_diff  <<  " secs)" << IO::endl;
  }

  intersection_->Status(VCellgausstype_);
}


/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::FindPositionDofSets(bool include_inner)
{

  backdis_->Comm().Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets (parallel)" );

  if(myrank_==0 and screenoutput_) IO::cout << "\t * 5/6 Cut_Positions_Dofsets (parallel) ...";

//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  GEO::CUT::Options options;
  intersection_->GetOptions(options);

  if ( options.FindPositions() )
  {

    GEO::CUT::Mesh & m = intersection_->NormalMesh();

    bool communicate = (backdis_->Comm().NumProc() > 1);

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::null;

    if (communicate)
    {
      cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( backdis_, m, *intersection_ ) );
    }

    // find inside and outside positions of nodes
    // first for mesh cut and distribute data in parallel, after that do the same for the level-set cut

    //--------------------------------------------
    // first, set the position for the mesh cut
    if(do_mesh_intersection_)
    {
      m.FindNodePositions();

      if (communicate) cut_parallel->CommunicateNodePositions();

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
    backdis_->Comm().Barrier();

    // find number and connection of dofsets at nodes from cut volumes
    intersection_->CreateNodalDofSet( include_inner, *backdis_);

    if (communicate) cut_parallel->CommunicateNodeDofSetNumbers(include_inner);

  }
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

  intersection_->DumpGmshNumDOFSets( str.str(), include_inner, *backdis_ );
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
      << backdis_->Comm().MyPID()
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
      << backdis_->Comm().MyPID()
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

GEO::CUT::ElementHandle * GEO::CutWizard::GetElement( DRT::Element * ele )
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





