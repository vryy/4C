/*!-----------------------------------------------------------------------------------------------*
\file xfem_fluidwizard.cpp

\brief class that provides the interface that bridges the cut libraries and fluid part of XFEM

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/geo_intersection.H"

#include "xfem_fluidwizard.H"
#include "xfem_fluiddofset.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

/*-------------------------------------------------------------*
* Cut routine for the new XFEM framework (XFSI and XFLUIDFLUID)*
*--------------------------------------------------------------*/
void XFEM::FluidWizard::Cut(  bool include_inner,             //!< perform cut within the structure
                              const Epetra_Vector & idispcol, //!< col vector holding interface displacements
                              std::string VCellgausstype,     //!< Gauss point generation method for Volumecell
                              std::string BCellgausstype,     //!< Gauss point generation method for Boundarycell
                              bool parallel,                  //!< use parallel cut algorithms
                              bool gmsh_output,               //!< print write gmsh output for cut
                              bool positions,                 //!< set inside and outside point, facet and volumecell positions
                              bool tetcellsonly,              //!< generate only tet cells
                              bool screenoutput,              //!< print screen output
                              bool cutinrefconf               //!< do not try to update node coordinates
                              )
{
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut" );

  if ( backdis_.Comm().MyPID() == 0 and screenoutput)
    IO::cout << "\nXFEM::FluidWizard::Cut:" << IO::endl;

  const double t_start = Teuchos::Time::wallTime();

  // set a new CutWizard based on the background discretization
  cut_ = Teuchos::rcp( new GEO::CutWizard( backdis_, false, 1 ) );
  cut_->SetFindPositions( positions );
  GEO::CutWizard & cw = *cut_;

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

  // 1. Add CutSides (sides of the cutterdiscretization)
  int numcutelements = cutterdis_.NumMyColElements();
  for ( int lid = 0; lid < numcutelements; ++lid )
  {
    DRT::Element * element = cutterdis_.lColElement(lid);

    const int numnode = element->NumNode();
    DRT::Node ** nodes = element->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      DRT::Node & node = *nodes[i];

      lm.clear();
      mydisp.clear();
      cutterdis_.Dof(&node, lm);

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
    
		// add the side of the cutter-discretization to the FluidWizard
    cw.AddCutSide( 0, element, xyze );
  }

  // 2. add background elements dependent on bounding box created by the CutSides in 1.
  int numbackelements = backdis_.NumMyColElements();
  for ( int k = 0; k < numbackelements; ++k )
  {
    DRT::Element * element = backdis_.lColElement( k );
    cw.AddElement( element );
  }

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
  if ( backdis_.Comm().MyPID() == 0  and screenoutput)
  {
    IO::cout << "\n\t ... Success (" << t_end  <<  " secs)\n" << IO::endl;
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
* creates a new fluid dofset                                   *
*--------------------------------------------------------------*/
Teuchos::RCP<XFEM::FluidDofSet> XFEM::FluidWizard::DofSet(int maxNumMyReservedDofs)
{
  return Teuchos::rcp( new FluidDofSet( this , maxNumMyReservedDofs, backdis_ ) );
}

/*-------------------------------------------------------------*
* get the cut wizard                                           *
*--------------------------------------------------------------*/
GEO::CutWizard & XFEM::FluidWizard::CutWizard()
{
  return *cut_;
}


/*-------------------------------------------------------------*
* get the elementhandle created within the cut                 *
*--------------------------------------------------------------*/
GEO::CUT::ElementHandle * XFEM::FluidWizard::GetElement(
    DRT::Element * ele
)
{
  return cut_->GetElement( ele );
}

/*-------------------------------------------------------------*
* get the sidehandle created within the cut                 *
*--------------------------------------------------------------*/
GEO::CUT::SideHandle * XFEM::FluidWizard::GetSide( std::vector<int> & nodeids)
{
  return cut_->GetSide( nodeids );
}

/*-------------------------------------------------------------*
* get the sidehandle created within the cut                 *
*--------------------------------------------------------------*/
GEO::CUT::SideHandle * XFEM::FluidWizard::GetSideHandle( int sid )
{
  return cut_->GetSide( sid );
}

GEO::CUT::SideHandle * XFEM::FluidWizard::GetCutSide( int sid, int mi )
{
  return cut_->GetCutSide(sid, mi);
}


/*-------------------------------------------------------------*
* get the node created within the cut                          *
*--------------------------------------------------------------*/
GEO::CUT::Node * XFEM::FluidWizard::GetNode(
    int nid
)
{
  return cut_->GetNode( nid );
}

