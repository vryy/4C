/*!-----------------------------------------------------------------------------------------------*
\file xfem_fluidwizard.cpp

\brief class that provides the interface that bridges the cut libraries and fluid part of XFEM

<pre>
Maintainer: Benedikt Schott and Magnus Winter
            schott@lnm.mw.tum.de, winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

//#include "../drt_io/io_control.H"
//#include "../drt_lib/drt_utils.H"
//#include "../drt_geometry/integrationcell.H"
//#include "../drt_geometry/geo_intersection.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_geometry/geo_meshintersection.H"
#include "../drt_geometry/geo_levelsetintersection.H"

#include "../drt_io/io_pstream.H"

#include "xfem_fluiddofset.H"

#include "xfem_fluidwizard.H"


/*-------------------------------------------------------------*
* creates a new fluid dofset                                   *
*--------------------------------------------------------------*/
Teuchos::RCP<XFEM::FluidDofSet> XFEM::FluidWizard::DofSet(int maxNumMyReservedDofsperNode)
{
  return Teuchos::rcp( new FluidDofSet( this , maxNumMyReservedDofsperNode, *backdis_ ) );
}

/*-------------------------------------------------------------*
* get the elementhandle created within the cut                 *
*--------------------------------------------------------------*/
GEO::CUT::ElementHandle * XFEM::FluidWizard::GetElement(
    DRT::Element * ele
)
{
  return parentcut_->GetElement( ele );
}

/*-------------------------------------------------------------*
* get the sidehandle created within the cut                 *
*--------------------------------------------------------------*/
GEO::CUT::SideHandle * XFEM::FluidWizard::GetSide( std::vector<int> & nodeids)
{
  return parentcut_->GetSide( nodeids );
}

/*-------------------------------------------------------------*
* get the sidehandle created within the cut                 *
*--------------------------------------------------------------*/
GEO::CUT::SideHandle * XFEM::FluidWizard::GetSideHandle( int sid )
{
  return parentcut_->GetSide( sid );
}

/*-------------------------------------------------------------*
* get the node created within the cut                          *
*--------------------------------------------------------------*/
GEO::CUT::Node * XFEM::FluidWizard::GetNode(
    int nid
)
{
  return parentcut_->GetNode( nid );
}

/*-------------------------------------------------------------*
* Cut routine for the new XFEM framework (XFSI and XFLUIDFLUID)*
*--------------------------------------------------------------*/
void XFEM::FluidWizardMesh::Cut(  bool include_inner,                         //!< perform cut within the structure
                              const Epetra_Vector & idispcol,             //!< col vector holding interface displacements
                              INPAR::CUT::VCellGaussPts VCellgausstype,   //!< Gauss point generation method for Volumecell
                              INPAR::CUT::BCellGaussPts BCellgausstype,   //!< Gauss point generation method for Boundarycell
                              bool parallel,                              //!< use parallel cut algorithms
                              bool gmsh_output,                           //!< print write gmsh output for cut
                              const Teuchos::RCP<Epetra_Vector> & dispnpcol, //!< background displacements
                              bool positions,                             //!< set inside and outside point, facet and volumecell positions
                              bool tetcellsonly,                          //!< generate only tet cells
                              bool screenoutput,                          //!< print screen output
                              bool cutinrefconf                           //!< do not try to update node coordinates
                              )
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizardMesh::Cut" );

  if ( backdis_->Comm().MyPID() == 0 and screenoutput)
    IO::cout << "\nXFEM::FluidWizardMesh::Cut:" << IO::endl;

  const double t_start = Teuchos::Time::wallTime();

  // set a new CutWizardMesh based on the background discretization
  cut_ = Teuchos::rcp( new GEO::CutWizardMesh( *backdis_, 1 ) );
  parentcut_ = cut_;
  cut_->SetFindPositions( positions );
  GEO::CutWizardMesh & cw = *cut_;

  {
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/6 --- Cut_Initialize" );

  if(backdis_->Comm().MyPID()==0 and screenoutput)
    IO::cout << "\n\t * 1/6 Cut_Initialize ...";

  std::vector<int> lm;
  std::vector<double> mydisp;

  bool is_crack = tip_nodes_.size() != 0;


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

  // 1. Add CutSides (sides of the cutterdiscretization)
  int numcutelements = cutterdis_->NumMyColElements();
  for ( int lid = 0; lid < numcutelements; ++lid )
  {
    DRT::Element * element = cutterdis_->lColElement(lid);

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      lm.clear();
      mydisp.clear();
      cutterdis_->Dof(&node, lm);

      LINALG::Matrix<3, 1> x( node.X() );

      if(!cutinrefconf)
      {
        if(lm.size() == 3) // case for BELE3 boundary elements
        {
          DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm);
        }
        else if(lm.size() == 4) // case for BELE3_4 boundary elements
        {
          // copy the first three entries for the displacement, the fourth entry should be zero if BELE3_4 is used for cutdis instead of BELE3
          std::vector<int> lm_red; // reduced local map
          lm_red.clear();
          for(int k=0; k< 3; k++) lm_red.push_back(lm[k]);

          DRT::UTILS::ExtractMyValues(idispcol,mydisp,lm_red);
        }
        else dserror("wrong number of dofs for node %i", lm.size());

        if (mydisp.size() != 3)
          dserror("we need 3 displacements here");

        LINALG::Matrix<3, 1> disp( &mydisp[0], true );



      // ------------------------------------------------------------------------------------------------
      // --- when simulating FSI with crack, nodes that represent crack are treated separately        ---
      // --- These nodes deform until the spring that connects them break. Only after the springs are ---
      // --- broken, we consider it as a real displacement and introduce crack opening there.         ---
      // --- this is advantageous. otherwise we need to introduce pseudo-elements there and should    ---
      // --- deal with them correctly. This becomes problematic when simulating extrinsic cohesive    ---
      // --- zone modeling                                                                            ---
      // ------------------------------------------------------------------------------------------------
      if( is_crack )
      {
        std::map<int, LINALG::Matrix<3,1> >::iterator itt = tip_nodes_.find( node.Id() );
        if( itt != tip_nodes_.end() )
        {
          disp = itt->second;
        }
      }

      //update x-position of cutter node for current time step (update with displacement)
      x.Update( 1, disp, 1 );
    }
    std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
  }

    // add the side of the cutter-discretization to the FluidWizardMesh
    cw.AddCutSide( 0, element, xyze );
  }

  // 2. add background elements dependent on bounding box created by the CutSides in 1.
  int numbackelements = backdis_->NumMyColElements();

  //Cast Discretisation to DiscretisationXFEM just used in case of AleBackground!

  Teuchos::RCP<DRT::DiscretizationXFEM>  xbackdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(backdis_);
   if(!cutinrefconf && dispnpcol != Teuchos::null)
   {
   if (xbackdis == Teuchos::null)
     dserror("XFEM::FluidWizardMesh::Cut: Cast to DiscretizationXFEM failed!");
   }

  for ( int k = 0; k < numbackelements; ++k )
  {
    DRT::Element * element = backdis_->lColElement(k);

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );
//
    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      LINALG::Matrix<3, 1> x( node.X() );

      if(!cutinrefconf && dispnpcol != Teuchos::null)
      {
        lm.clear();
        mydisp.clear();
        xbackdis->InitialDof(&node, lm); //to get all dofs of background (also not active ones at the moment!)

        if(lm.size() == 3) // case used actually?
        {
          DRT::UTILS::ExtractMyValues(*dispnpcol,mydisp,lm);
        }
        else if(lm.size() == 4) // case xFluid ... just take the first three
        {
          // copy the first three entries for the displacement, the fourth entry an all others
          std::vector<int> lm_red; // reduced local map
          lm_red.clear();
          for(int k=0; k< 3; k++) lm_red.push_back(lm[k]);

          DRT::UTILS::ExtractMyValues(*dispnpcol,mydisp,lm_red);
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
    cw.AddElement( element , xyze );
  }

  // build the bounding volume tree for the collision detection in the context of the selfcut
//  cw.BuildBVTree();

  // build the static search tree for the collision detection
  cw.BuildStaticSearchTree();

  const double t_mid = Teuchos::Time::wallTime()-t_start;
  if ( backdis_->Comm().MyPID() == 0  and screenoutput)
  {
    IO::cout << "\t\t\t... Success (" << t_mid  <<  " secs)" << IO::endl;
  }
  }
  // wirtz 08/14:
  // preprocessing: everything above should only be done once in a simulation; so it should be moved before the time loop into a preprocessing step
  // runtime:       everything below should be done in every Newton increment

  // run the (parallel) Cut
  if(parallel)
  {
    cw.CutParallel( include_inner, VCellgausstype, BCellgausstype,tetcellsonly,screenoutput );
  }
  else
  {
    dserror("the non-parallel cutwizard does not support the DofsetNEW framework");
//    cw.Cut( include_inner, VCellgausstype, BCellgausstype );
  }

  // cleanup

  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( backdis_->Comm().MyPID() == 0  and screenoutput)
  {
    IO::cout << "\n\t\t\t\t\t\t\t... Success (" << t_end  <<  " secs)\n" << IO::endl;
  }

  if(gmsh_output) cw.DumpGmshNumDOFSets(include_inner);

#ifdef DEBUG
  cw.PrintCellStats();
#endif

  if(gmsh_output)
  {
    cw.DumpGmshIntegrationCells();
    cw.DumpGmshVolumeCells( include_inner );
  }

}

/*-------------------------------------------------------------*
* get the cut wizard                                           *
*--------------------------------------------------------------*/
GEO::CutWizardMesh & XFEM::FluidWizardMesh::CutWizard()
{
  return *cut_;
}

GEO::CUT::SideHandle * XFEM::FluidWizardMesh::GetCutSide( int sid, int mi )
{
  return cut_->GetCutSide(sid, mi);
}

/*--------------------------------------------------------------------------*
* Cut routine for the new XFEM framework (TPFX)* utilizing level set for cut.
*---------------------------------------------------------------------------*/
void XFEM::FluidWizardLevelSet::Cut(
                              bool include_inner,                         //!< perform cut within the structure
                              const Epetra_Vector & phinpnode,            //!< node based values for the level set function
                              INPAR::CUT::VCellGaussPts VCellgausstype,   //!< Gauss point generation method for Volumecell
                              INPAR::CUT::BCellGaussPts BCellgausstype,   //!< Gauss point generation method for Boundarycell
                              bool parallel,                              //!< use parallel cut algorithms
                              bool gmsh_output,                           //!< print write gmsh output for cut
                              const Teuchos::RCP<Epetra_Vector> & dispnpcol, //!< background displacements
                              bool positions,                             //!< set inside and outside point, facet and volumecell positions
                              bool tetcellsonly,                          //!< generate only tet cells
                              bool screenoutput,                          //!< print screen output
                              bool cutinrefconf                           //!< do not try to perform coordinate update
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizardLevelSet::Cut" );

  //fluiddis_ = backdis_

  if ( backdis_->Comm().MyPID() == 0 and screenoutput)
    IO::cout << "\nXFEM::FluidWizardLevelSet::Cut:" << IO::endl;

  cut_ = Teuchos::rcp( new GEO::CutWizardLevelSet( *backdis_ ) );
  parentcut_ = cut_;
  cut_->SetFindPositions( positions );

  const double t_start = Teuchos::Time::wallTime();

  if (dispnpcol != Teuchos::null)
    dserror("XFEM::FluidWizardLevelSet::Cut: AleDisplacements not implemented here!");

  // Loop over all Elements to find cut elements and add them to the LevelsetIntersection class
  // Brute force method.
  int numelements = backdis_->NumMyColElements();

  std::vector<double> myphinp;

  for ( int lid = 0; lid < numelements; ++lid )
  {
    myphinp.clear();

    DRT::Element * element = backdis_->lColElement(lid);

    DRT::UTILS::ExtractMyNodeBasedValues(element, myphinp, phinpnode);
    cut_->AddElement(element,myphinp,include_inner);
  }

  // run the (parallel) Cut
  if(parallel)
  {
    cut_->CutParallel( include_inner, VCellgausstype, BCellgausstype,tetcellsonly,screenoutput );
  }
  else
  {
//    dserror("the non-parallel cutwizard does not support the DofsetNEW framework");
    cut_->Cut( include_inner, VCellgausstype, BCellgausstype );
  }

  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( backdis_->Comm().MyPID() == 0  and screenoutput)
  {
    IO::cout << "\n\t ... Success (" << t_end  <<  " secs)\n" << IO::endl;
  }

  if(gmsh_output) cut_->DumpGmshNumDOFSets(include_inner);

#ifdef DEBUG
  cut_->PrintCellStats();
#endif

  if(gmsh_output)
  {
    cut_->DumpGmshIntegrationCells();
    cut_->DumpGmshVolumeCells( include_inner );
  }

}

/*-------------------------------------------------------------*
* get the cut wizard                                           *
*--------------------------------------------------------------*/
GEO::CutWizardLevelSet & XFEM::FluidWizardLevelSet::CutWizard()
{
  return *cut_;
}
