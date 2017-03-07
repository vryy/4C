/*----------------------------------------------------------------------*/
/*!
\file xfluid_outputservice.cpp

\brief Service class for XFluid-related output.

\level 2

<pre>
\maintainer  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include "xfluid_outputservice.H"
#include "xfluid_state_creator.H"
#include "xfluid_state.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../linalg/linalg_utils.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_discretization_utils.H"
#include "../drt_xfem/xfem_edgestab.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_cutwizard.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_fluid/fluid_utils.H"

FLD::XFluidOutputService::XFluidOutputService(
  const Teuchos::RCP<DRT::DiscretizationXFEM>& discret,
  const Teuchos::RCP<XFEM::ConditionManager>&  cond_manager
) :
  discret_(discret),
  cond_manager_(cond_manager),
  firstoutputofrun_(true),
  restart_count_(0)
{
  // Vector & map extractor for paraview output,
  // mapped to initial fluid dofmap
  dofset_out_ = Teuchos::rcp(new DRT::IndependentDofSet());
  velpressplitter_out_ = Teuchos::rcp(new LINALG::MapExtractor());
  PrepareOutput();
}

void FLD::XFluidOutputService::PrepareOutput()
{
  dofset_out_->Reset();
  dofset_out_->AssignDegreesOfFreedom(*discret_,0,0);
  const int ndim = DRT::Problem::Instance()->NDim();
  // split based on complete fluid field (standard splitter that handles one dofset)
  FLD::UTILS::SetupFluidSplit(*discret_,*dofset_out_,ndim,*velpressplitter_out_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_->DofRowMap(),true);
}

void FLD::XFluidOutputService::Output(
  int step,
  double time,
  bool write_restart_data,
  Teuchos::RCP<const FLD::XFluidState> state,
  Teuchos::RCP<Epetra_Vector>    dispnp,
  Teuchos::RCP<Epetra_Vector>    gridvnp
)
{
  discret_->Writer()->NewStep(step,time);

  // create vector according to the initial row map holding all standard fluid unknowns
  outvec_fluid_->PutScalar(0.0);

  const Epetra_Map* dofrowmap  = dofset_out_->DofRowMap(); // original fluid unknowns
  const Epetra_Map* xdofrowmap = discret_->DofRowMap();   // fluid unknown for current cut

  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    // get row node via local id
    const DRT::Node* xfemnode = discret_->lRowNode(i);

    // the initial dofset contains the original dofs for each row node
    const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));

    // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
    // else copy the right nodes
    const std::vector<int> gdofs_current(discret_->Dof(0,xfemnode));

    if(gdofs_current.size() == 0)
    {
      // std::cout << "no dofs available->hole" << std::endl;
    }
    else if(gdofs_current.size() == gdofs_original.size())
    {
      //std::cout << "same number of dofs available" << std::endl;
    }
    else if(gdofs_current.size() > gdofs_original.size())
    {
      //std::cout << "more dofs available->decide" << std::endl;
    }
    else std::cout << "decide which dofs can be copied and which have to be set to zero" << std::endl;

    if(gdofs_current.size() == 0)
    {
      size_t numdof = gdofs_original.size();
      // no dofs for this node... must be a hole or somethin'
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        //std::cout << dofrowmap->LID(gdofs[idof]) << std::endl;
        (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
      }
    }
    else if(gdofs_current.size() == gdofs_original.size())
    {
      size_t numdof = gdofs_original.size();
      // copy all values
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        //std::cout << dofrowmap->LID(gdofs[idof]) << std::endl;
        (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = (*state->velnp_)[xdofrowmap->LID(gdofs_current[idof])];
      }
    }
    else if(gdofs_current.size() % gdofs_original.size() == 0) //multiple dofsets
    {
      // if there are multiple dofsets we write output for the standard dofset
      GEO::CUT::Node* node = state->Wizard()->GetNode(xfemnode->Id());

      const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > & dofcellsets = node->NodalDofSets();

      int nds = 0;
      bool is_std_set = false;

      // find the standard dofset
      for(std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> >::const_iterator cellsets= dofcellsets.begin();
          cellsets!=dofcellsets.end();
          cellsets++)
      {
        is_std_set = (*cellsets)->Is_Standard_DofSet();

        if(is_std_set) break;

        nds++;
      }

      size_t numdof = gdofs_original.size();
      size_t offset = 0; // no offset in case of no std-dofset means take the first dofset for output

      if(is_std_set)
      {
        offset = numdof*nds;   // offset to start the loop over dofs at the right nds position
      }

      // copy all values
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = (*state->velnp_)[xdofrowmap->LID(gdofs_current[offset+idof])];
      }


      // if there are just ghost values, then we write output for the set with largest physical fluid volume
    }
    else dserror("unknown number of dofs for output");

  };

  // output (hydrodynamic) pressure for visualization

  Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_out_->ExtractCondVector(outvec_fluid_);

  discret_->Writer()->WriteVector("velnp", outvec_fluid_);
  discret_->Writer()->WriteVector("pressure", pressure);

  if (dispnp != Teuchos::null)
  {
    if (gridvnp == Teuchos::null)
      dserror("Missing grid velocities for ALE-xfluid!");

    //write ale displacement for t^{n+1}
    Teuchos::RCP<Epetra_Vector> dispnprm = Teuchos::rcp(new Epetra_Vector(*dispnp));
    dispnprm->ReplaceMap(outvec_fluid_->Map()); //to get dofs starting by 0 ...
    discret_->Writer()->WriteVector("dispnp", dispnprm);

    //write grid velocity for t^{n+1}
    Teuchos::RCP<Epetra_Vector> gridvnprm = Teuchos::rcp(new Epetra_Vector(*gridvnp));
    gridvnprm->ReplaceMap(outvec_fluid_->Map()); //to get dofs starting by 0 ...
    discret_->Writer()->WriteVector("gridv", gridvnprm);

    //write convective velocity for t^{n+1}
    Teuchos::RCP<Epetra_Vector> convvel = Teuchos::rcp(new Epetra_Vector(outvec_fluid_->Map(),true));
    convvel->Update(1.0,*outvec_fluid_,-1.0,*gridvnprm,0.0);
    discret_->Writer()->WriteVector("convel", convvel);
  }

  discret_->Writer()->WriteElementData(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    if(discret_->Comm().MyPID() == 0) IO::cout << "---  write restart... " << IO::endl;

    restart_count_++;

    // velocity/pressure vector
    discret_->Writer()->WriteVector("velnp_res",state->velnp_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    discret_->Writer()->WriteVector("accnp_res",state->accnp_);
    discret_->Writer()->WriteVector("accn_res", state->accn_);
    discret_->Writer()->WriteVector("veln_res", state->veln_);
    discret_->Writer()->WriteVector("velnm_res",state->velnm_);

    if (dispnp != Teuchos::null)
    {
      //write ale displacement for t^{n+1} on full background
      discret_->Writer()->WriteVector("full_dispnp_res", dispnp);

      //write grid velocity for t^{n+1} on full background
      discret_->Writer()->WriteVector("full_gridvnp_res", gridvnp);
    }
  }

  //-----------------------------------------------------------
  // REMARK on "Why to clear the MapCache" for restarts
  //-----------------------------------------------------------
  // every time, when output-vectors are written based on a new(!), still unknown map
  // the map is stored in a mapstack (io.cpp, WriteVector-routine) for efficiency in standard applications.
  // However, in the XFEM for each timestep or FSI-iteration we have to cut and the map is created newly
  // (done by the FillComplete call).
  // This is the reason why the MapStack increases and the storage is overwritten for large problems.
  // Hence, we have clear the MapCache in regular intervals of written restarts.
  // In case of writing paraview-output, here, we use a standard map which does not change over time, that's okay.
  // For the moment, restart_count = 5 is set quite arbitrary, in case that we need for storage, we have to reduce this number

  if(restart_count_ == 5)
  {
    if(discret_->Comm().MyPID() == 0) IO::cout << "\t... Clear MapCache after " << restart_count_ << " written restarts." << IO::endl;

    discret_->Writer()->ClearMapCache(); // clear the output's map-cache
    restart_count_ = 0;
  }

  //-----------------------------------------------------------
  // write paraview output for cutter discretization
  //-----------------------------------------------------------
  cond_manager_->Output(step, time, write_restart_data);
}

FLD::XFluidOutputServiceGmsh::XFluidOutputServiceGmsh(
  Teuchos::ParameterList& params_xfem,
  const Teuchos::RCP<DRT::DiscretizationXFEM>& discret,
  const Teuchos::RCP<XFEM::ConditionManager>&  cond_manager,
  const bool include_inner
) : XFluidOutputService(discret, cond_manager),
  gmsh_sol_out_          ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_SOL_OUT")),
  gmsh_ref_sol_out_      (false),
  gmsh_debug_out_        ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT")),
  gmsh_debug_out_screen_ ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT_SCREEN")),
  gmsh_EOS_out_          ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_EOS_OUT")),
  gmsh_discret_out_      ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DISCRET_OUT")),
  gmsh_step_diff_        (500),
  VolumeCellGaussPointBy_(DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY")),
  include_inner_         (include_inner)
{
  if (! (bool)DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_GMSH"))
    dserror("If GMSH output is globally deactivated, don't create an instance of XFluidOutputServiceGmsh!");
};

void FLD::XFluidOutputServiceGmsh::GmshSolutionOutput(
  const std::string & filename_base,          ///< name for output file
  int step,                                   ///< step number
  const Teuchos::RCP<FLD::XFluidState>& state,///< state
  int count                                   ///< counter (no counter for standard solution output : -1)
)
{
  if (!gmsh_sol_out_ && !gmsh_ref_sol_out_)
    return;

  Teuchos::RCP<const Epetra_Vector> output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_,state->Velnp());

  Teuchos::RCP<const Epetra_Vector> output_col_acc = Teuchos::null;

  if (state->Accnp() != Teuchos::null)
  {
    output_col_acc = DRT::UTILS::GetColVersionOfRowVector(discret_,state->Accnp());
  }

  // no counter for standard solution output : -1
  const std::string prefix("SOL");
  GmshOutput( filename_base, prefix, step, count , state->Wizard(), output_col_vel, output_col_acc );
  cond_manager_->GmshOutput(filename_base, step, gmsh_step_diff_, gmsh_debug_out_screen_);
}

void FLD::XFluidOutputServiceGmsh::GmshSolutionOutputPrevious(
  const std::string & filename_base,          ///< name for output file
  int step,                                   ///< step number
  const Teuchos::RCP<FLD::XFluidState>& state,///< state
  int count
)
{
  if (!gmsh_sol_out_ && !gmsh_ref_sol_out_)
    return;

  Teuchos::RCP<const Epetra_Vector> output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_,state->Veln());

  Teuchos::RCP<const Epetra_Vector> output_col_acc = Teuchos::null;

  if (state->Accn() != Teuchos::null)
  {
    output_col_acc = DRT::UTILS::GetColVersionOfRowVector(discret_,state->Accn());
  }

  const std::string prefix("ref_SOL");
  GmshOutput( filename_base, prefix, step, count , state->Wizard(), output_col_vel, output_col_acc );
}

void FLD::XFluidOutputServiceGmsh::GmshSolutionOutputDebug(
  const std::string & filename_base,          ///< name for output file
  int step,                                   ///< step number
  int count,                                  ///< counter for iterations within a global time step
  const Teuchos::RCP<FLD::XFluidState>& state ///< state
)
{
  if (!gmsh_debug_out_)
    return;

  Teuchos::RCP<const Epetra_Vector> output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_,state->Velnp());
  const std::string prefix("SOL");
  GmshOutput( filename_base, prefix, step, count , state->Wizard(), output_col_vel );
}

void FLD::XFluidOutputServiceGmsh::GmshResidualOutputDebug(
  const std::string & filename_base,          ///< name for output file
  int step,                                   ///< step number
  int count,                                  ///< counter for iterations within a global time step
  const Teuchos::RCP<FLD::XFluidState>& state ///< state
)
{
  if (!gmsh_debug_out_)
    return;

  Teuchos::RCP<const Epetra_Vector> output_col_residual = DRT::UTILS::GetColVersionOfRowVector(discret_,state->Residual());
  const std::string prefix("RES");
  GmshOutput( filename_base, prefix, step, count , state->Wizard(), output_col_residual );
}

void FLD::XFluidOutputServiceGmsh::GmshIncrementOutputDebug(
  const std::string & filename_base,          ///< name for output file
  int step,                                   ///< step number
  int count,                                  ///< counter for iterations within a global time step
  const Teuchos::RCP<FLD::XFluidState>& state ///< state
)
{
  if (!gmsh_debug_out_)
    return;

  Teuchos::RCP<const Epetra_Vector> output_col_incvel = DRT::UTILS::GetColVersionOfRowVector(discret_,state->IncVel());
  const std::string prefix("INC");
  GmshOutput( filename_base, prefix, step, count , state->Wizard(), output_col_incvel );
}

void FLD::XFluidOutputServiceGmsh::GmshOutput(
  const std::string & filename_base,          ///< name for output file
  const std::string & prefix,                 ///< data prefix
  int step,                                   ///< step number
  int count,                                  ///< counter for iterations within a global time step
  const Teuchos::RCP<GEO::CutWizard>& wizard, ///< cut wizard
  Teuchos::RCP<const Epetra_Vector> vel,      ///< vector holding velocity and pressure dofs
  Teuchos::RCP<const Epetra_Vector> acc       ///< vector holding acceleration
)
{
  // Todo: should be private
  int myrank = discret_->Comm().MyPID();

  if(myrank==0) std::cout << "\n\t ... writing Gmsh output...\n" << std::flush;

  bool screen_out = gmsh_debug_out_screen_;

  // output for Element and Node IDs
  std::ostringstream filename_base_vel;
  if(count > -1) filename_base_vel << filename_base << "_" << count << "_vel";
  else        filename_base_vel << filename_base << "_vel";
  const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_vel(filename_vel.c_str());
  gmshfilecontent_vel.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_vel.precision(16);

  std::ostringstream filename_base_press;
  if(count > -1) filename_base_press << filename_base << "_" << count << "_press";
  else        filename_base_press << filename_base << "_press";
  const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_press(filename_press.c_str());
  gmshfilecontent_press.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_press.precision(16);

  std::ostringstream filename_base_acc;
  if(count > -1) filename_base_acc << filename_base << "_" << count << "_acc";
  else        filename_base_acc << filename_base << "_acc";
  const std::string filename_acc = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_acc.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_acc(filename_acc.c_str());
  gmshfilecontent_acc.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_acc.precision(16);

  std::ostringstream filename_base_bound;
  if(count > -1) filename_base_bound << filename_base << "_" << count << "_bound";
  else        filename_base_bound << filename_base << "_bound";
  const std::string filename_bound = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_bound.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_bound(filename_bound.c_str());
  gmshfilecontent_bound.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_bound.precision(16);


  // output for Element and Node IDs
  std::ostringstream filename_base_vel_ghost;
  if(count > -1) filename_base_vel_ghost << filename_base << "_" << count << "_vel_ghost";
  else        filename_base_vel_ghost << filename_base << "_vel_ghost";
  const std::string filename_vel_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel_ghost.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_vel_ghost(filename_vel_ghost.c_str());
  gmshfilecontent_vel_ghost.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_vel_ghost.precision(16);

  std::ostringstream filename_base_press_ghost;
  if(count > -1) filename_base_press_ghost << filename_base << "_" << count << "_press_ghost";
  else        filename_base_press_ghost << filename_base << "_press_ghost";
  const std::string filename_press_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press_ghost.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_press_ghost(filename_press_ghost.c_str());
  gmshfilecontent_press_ghost.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_press_ghost.precision(16);

  std::ostringstream filename_base_acc_ghost;
  if(count > -1) filename_base_acc_ghost << filename_base << "_" << count << "_acc_ghost";
  else        filename_base_acc_ghost << filename_base << "_acc_ghost";
  const std::string filename_acc_ghost = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_acc_ghost.str(), step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_acc_ghost(filename_acc_ghost.c_str());
  gmshfilecontent_acc_ghost.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_acc_ghost.precision(16);


  if(count > -1) // for residual output
  {
     gmshfilecontent_vel         << "View \"" << prefix << "vel "         << count << "\" {\n";
     gmshfilecontent_press       << "View \"" << prefix << "press "       << count << "\" {\n";
     gmshfilecontent_bound       << "View \"" << prefix << "side-normal " << count << "\" {\n";
     gmshfilecontent_vel_ghost   << "View \"" << prefix << "vel_ghost "   << count << "\" {\n";
     gmshfilecontent_press_ghost << "View \"" << prefix << "press_ghost " << count << "\" {\n";
  }
  else
  {
     gmshfilecontent_vel         << "View \"" << prefix << "vel "         << "\" {\n";
     gmshfilecontent_press       << "View \"" << prefix << "press "       << "\" {\n";
     gmshfilecontent_acc         << "View \"" << prefix << "acc  "        << "\" {\n";
     gmshfilecontent_bound       << "View \"" << prefix << "side-normal " << "\" {\n";
     gmshfilecontent_vel_ghost   << "View \"" << prefix << "vel_ghost "   << "\" {\n";
     gmshfilecontent_press_ghost << "View \"" << prefix << "press_ghost " << "\" {\n";
     gmshfilecontent_acc_ghost   << "View \"" << prefix << "acc_ghost  "  << "\" {\n";
  }

  const int numrowele = (discret_)->NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = (discret_)->lRowElement(i);
    GEO::CUT::ElementHandle * e = wizard->GetElement( actele );

    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;

      e->GetVolumeCellsDofSets( cell_sets, nds_sets, include_inner_);


      if(e->IsIntersected())
      {
        int set_counter = 0;

        for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
            s!=cell_sets.end();
            s++)
        {
          GEO::CUT::plain_volumecell_set & cells = *s;

          std::vector<int> & nds = nds_sets[set_counter];


          for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
          {
            GEO::CUT::VolumeCell * vc = *i;

            if ( e->IsCut() )
            {
              GmshOutputVolumeCell( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, e, vc, nds, vel, acc );
              if ( vc->Position()==GEO::CUT::Point::outside )
              {
                if (cond_manager_->HasMeshCoupling())
                  GmshOutputBoundaryCell( *discret_, gmshfilecontent_bound, vc, wizard );
              }
            }
            else
            {
              GmshOutputElement( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
            }
            GmshOutputElement( *discret_, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, nds, vel, acc );
          }
          set_counter += 1;
        }
      }
      else if(cell_sets.size() == 1) // non intersected but one set
      {
        std::vector<int> & nds = nds_sets[0];

        // one standard uncut physical element
        GmshOutputElement( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
      }
      else
      {
        // ghost element
      }
    } // element handle
    else // no element handle
    {
      std::vector<int> nds; // empty vector
      GmshOutputElement( *discret_, gmshfilecontent_vel, gmshfilecontent_press, gmshfilecontent_acc, actele, nds, vel, acc );
      GmshOutputElement( *discret_, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost, gmshfilecontent_acc_ghost, actele, nds, vel, acc );
    }
  } // loop elements

  gmshfilecontent_vel   << "};\n";
  gmshfilecontent_press << "};\n";
  if(count == -1) gmshfilecontent_acc   << "};\n";
  gmshfilecontent_vel_ghost   << "};\n";
  gmshfilecontent_press_ghost << "};\n";
  if(count == -1) gmshfilecontent_acc_ghost   << "};\n";
  gmshfilecontent_bound << "};\n";

  gmshfilecontent_vel.close();
  gmshfilecontent_press.close();
  gmshfilecontent_acc.close();
  gmshfilecontent_bound.close();
  gmshfilecontent_vel_ghost.close();
  gmshfilecontent_press_ghost.close();
  gmshfilecontent_acc_ghost.close();

  if(myrank==0) std::cout << " done\n" << std::flush;
}

/// Gmsh output function for elements without an GEO::CUT::ElementHandle
void FLD::XFluidOutputServiceGmsh::GmshOutputElement(
  DRT::Discretization & discret, ///< background fluid discretization
  std::ofstream & vel_f,         ///< output file stream for velocity
  std::ofstream & press_f,       ///< output file stream for pressure
  std::ofstream & acc_f,         ///< output file stream for acceleration
  DRT::Element * actele,         ///< element
  std::vector<int> & nds,        ///< vector holding the nodal dofsets
  Teuchos::RCP<const Epetra_Vector> vel, ///< vector holding velocity and pressure dofs
  Teuchos::RCP<const Epetra_Vector> acc  ///< vector holding acceleration
)
{
  vel_f.setf(std::ios::scientific,std::ios::floatfield);
  vel_f.precision(16);

  press_f.setf(std::ios::scientific,std::ios::floatfield);
  press_f.precision(16);

  acc_f.setf(std::ios::scientific,std::ios::floatfield);
  acc_f.precision(16);

  //output for accvec ?
  const bool acc_output(acc != Teuchos::null);

  DRT::Element::LocationArray la( 1 );


  if(nds.size()!=0) // for element output of ghost values
  {
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(discret,nds,la,false);
  }
  else
  {
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(discret,la,false);
  }

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  std::vector<double> m_acc(la[0].lm_.size());
  if(acc_output)
  {
    DRT::UTILS::ExtractMyValues(*acc,m_acc,la[0].lm_);
  }

  int numnode = 0;

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    numnode=8;
    vel_f << "VH(";
    press_f << "SH(";
    if(acc_output) acc_f << "VH(";
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    numnode=6;
    vel_f << "VI(";
    press_f << "SI(";
    if(acc_output) acc_f << "VI(";
    break;
  case DRT::Element::pyramid5:
    numnode=5;
    vel_f << "VY(";
    press_f << "SY(";
    if(acc_output) acc_f << "VY(";
    break;
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    numnode=4;
    vel_f << "VS(";
    press_f << "SS(";
    if(acc_output) acc_f << "VS(";
    break;
  default:
  {
    dserror( "unsupported shape" );
    break;
  }
  }

  for ( int i=0; i<numnode; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
      if(acc_output) acc_f << ",";
    }
    const double * x = actele->Nodes()[i]->X();
    vel_f   << x[0] << "," << x[1] << "," << x[2];
    press_f << x[0] << "," << x[1] << "," << x[2];
    if(acc_output) acc_f << x[0] << "," << x[1] << "," << x[2];
  }
  vel_f << "){";
  press_f << "){";
  if(acc_output) acc_f << "){";

  for ( int i=0; i<numnode; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
      if(acc_output) acc_f << ",";
    }
    int j = 4*i;
    vel_f   << m[j] << "," << m[j+1] << "," << m[j+2];
    press_f << m[j+3];
    if(acc_output) acc_f   << m_acc[j] << "," << m_acc[j+1] << "," << m_acc[j+2];
  }

  vel_f << "};\n";
  press_f << "};\n";
  if(acc_output) acc_f << "};\n";
}

/// Gmsh output function for volumecells
void FLD::XFluidOutputServiceGmsh::GmshOutputVolumeCell(
  DRT::Discretization & discret,    ///< background fluid discretization
  std::ofstream & vel_f,            ///< output file stream for velocity
  std::ofstream & press_f,          ///< output file stream for pressure
  std::ofstream & acc_f,            ///< output file stream for acceleration
  DRT::Element * actele,            ///< element
  GEO::CUT::ElementHandle * e,      ///<elementhandle
  GEO::CUT::VolumeCell * vc,        ///< volumecell
  const std::vector<int> & nds,     ///< vector holding the nodal dofsets
  Teuchos::RCP<const Epetra_Vector> velvec,  ///< vector holding velocity and pressure dofs
  Teuchos::RCP<const Epetra_Vector> accvec   ///< vector holding acceleration
)
{
  vel_f.setf(std::ios::scientific,std::ios::floatfield);
  vel_f.precision(16);

  press_f.setf(std::ios::scientific,std::ios::floatfield);
  press_f.precision(16);

  acc_f.setf(std::ios::scientific,std::ios::floatfield);
  acc_f.precision(16);



  //output for accvec ?
  bool acc_output = true;
  if(accvec == Teuchos::null) acc_output=false;

  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,nds,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*velvec,m,la[0].lm_);

  std::vector<double> m_acc(la[0].lm_.size());
  if(acc_output)
  {
    DRT::UTILS::ExtractMyValues(*accvec,m_acc,la[0].lm_);
  }

  Epetra_SerialDenseMatrix vel( 3, actele->NumNode() );
  Epetra_SerialDenseMatrix press( 1, actele->NumNode() );
  Epetra_SerialDenseMatrix acc( 3, actele->NumNode() );

  for ( int i=0; i<actele->NumNode(); ++i )
  {
    vel( 0, i ) = m[4*i+0];
    vel( 1, i ) = m[4*i+1];
    vel( 2, i ) = m[4*i+2];

    press( 0, i ) = m[4*i+3];

    if(acc_output)
    {
     acc( 0, i ) = m_acc[4*i+0];
     acc( 1, i ) = m_acc[4*i+1];
     acc( 2, i ) = m_acc[4*i+2];
    }
  }

  // facet based output for cut volumes
  // integrationcells are not available because tessellation is not used
  if( VolumeCellGaussPointBy_!= INPAR::CUT::VCellGaussPts_Tessellation )
  {
    const GEO::CUT::plain_facet_set & facete = vc->Facets();
    for(GEO::CUT::plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
      // split facet into tri and quad cell
      GEO::CUT::Facet *fe = *i;
      std::vector<std::vector<GEO::CUT::Point*> > split;
      std::vector<GEO::CUT::Point*> corners = fe->CornerPoints();

      if( corners.size()==3 ) // only Tri can be used directly. Quad may be concave
        split.push_back( corners );
      else
      {
        if( !fe->IsFacetSplit() )
          fe->SplitFacet( fe->CornerPoints() );
         split = fe->GetSplitCells();
      }

      for( std::vector<std::vector<GEO::CUT::Point*> >::const_iterator j=split.begin();
                                                                       j!=split.end();j++ )
      {
        std::vector<GEO::CUT::Point*> cell = *j;

        switch ( cell.size() )
        {
        case 3:
          vel_f << "VT(";
          press_f << "ST(";
          if(acc_output) acc_f << "VT(";
          break;
        case 4:
          vel_f << "VQ(";
          press_f << "SQ(";
          if(acc_output) acc_f << "VQ(";
          break;
        default:
        {
          dserror( "splitting facets failed" );
          break;
        }
        }

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
            if(acc_output) acc_f << ",";
          }
          const double * x = cell[k]->X();
          vel_f   << x[0] << "," << x[1] << "," << x[2];
          press_f << x[0] << "," << x[1] << "," << x[2];
          if(acc_output) acc_f   << x[0] << "," << x[1] << "," << x[2];
        }
        vel_f << "){";
        press_f << "){";
        if(acc_output) acc_f << "){";

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          LINALG::Matrix<3,1> v( true );
          LINALG::Matrix<1,1> p( true );
          LINALG::Matrix<3,1> a( true );

          GEO::CUT::Point * point = cell[k];
          const LINALG::Matrix<3,1> & rst = e->LocalCoordinates( point );

          switch ( actele->Shape() )
          {
          case DRT::Element::hex8:
          {
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );
            LINALG::Matrix<3,numnodes> acceleration( acc, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
            break;
          }
          case DRT::Element::hex20:
          {
            // TODO: check the output for hex20
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex20 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );
            LINALG::Matrix<3,numnodes> acceleration( acc, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
            break;
          }
          case DRT::Element::hex27:
          {
            // TODO: check the output for hex27
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex27 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );
            LINALG::Matrix<3,numnodes> acceleration( acc, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
            break;
          }
          default:
          {
            dserror( "unsupported shape" );
            break;
          }
          }

          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
            if(acc_output) acc_f << ",";
          }
          vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
          press_f << p( 0 );
          if(acc_output) acc_f   << a( 0 ) << "," << a( 1 ) << "," << a( 2 );
        }

        vel_f << "};\n";
        press_f << "};\n";
        if(acc_output) acc_f << "};\n";
      }
    }
  }

  // integrationcells based output for tessellation
  else
  {
    const GEO::CUT::plain_integrationcell_set & intcells = vc->IntegrationCells();
    for ( GEO::CUT::plain_integrationcell_set::const_iterator i=intcells.begin();
          i!=intcells.end();
          ++i )
    {
      GEO::CUT::IntegrationCell * ic = *i;

      const std::vector<GEO::CUT::Point*> & points = ic->Points();
  //    Epetra_SerialDenseMatrix values( 4, points.size() );

      switch ( ic->Shape() )
      {
      case DRT::Element::hex8:
        vel_f << "VH(";
        press_f << "SH(";
        if(acc_output) acc_f << "VH(";
        break;
      case DRT::Element::tet4:
        vel_f << "VS(";
        press_f << "SS(";
        if(acc_output) acc_f << "VS(";
        break;
      default:
      {
        dserror( "unsupported shape" );
        break;
      }
      }

      for ( unsigned i=0; i<points.size(); ++i )
      {
        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
          if(acc_output) acc_f << ",";
        }
        const double * x = points[i]->X();
        vel_f   << x[0] << "," << x[1] << "," << x[2];
        press_f << x[0] << "," << x[1] << "," << x[2];
        if(acc_output) acc_f   << x[0] << "," << x[1] << "," << x[2];
      }
      vel_f << "){";
      press_f << "){";
      if(acc_output) acc_f << "){";

      for ( unsigned i=0; i<points.size(); ++i )
      {
        LINALG::Matrix<3,1> v( true );
        LINALG::Matrix<1,1> p( true );
        LINALG::Matrix<3,1> a( true );

        GEO::CUT::Point * point = points[i];
        const LINALG::Matrix<3,1> & rst = e->LocalCoordinates( point );

        switch ( actele->Shape() )
        {
        case DRT::Element::hex8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::hex20:
        {
          // TODO: check the output for hex20
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex20 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::hex27:
        {
          // TODO: check the output for hex27
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex27 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::wedge6:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge6 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::wedge15:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge15>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge15 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::tet4:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::tet4 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        case DRT::Element::tet10:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::tet10 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );
          LINALG::Matrix<3,numnodes> acceleration( acc, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          if(acc_output) a.Multiply( 1, acceleration, funct, 1 );
          break;
        }
        default:
        {
          dserror( "unsupported shape" );
          break;
        }
        }


        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
          if(acc_output) acc_f << ",";
        }
        vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
        press_f << p( 0 );
        if(acc_output) acc_f   << a( 0 ) << "," << a( 1 ) << "," << a( 2 );
      }

      vel_f << "};\n";
      press_f << "};\n";
      if(acc_output) acc_f << "};\n";
    }
  }
}

/// Gmsh output function for boundarycells
void FLD::XFluidOutputServiceGmsh::GmshOutputBoundaryCell(
  DRT::Discretization & discret,    ///< background fluid discretization
  std::ofstream & bound_f,          ///< output file stream for boundary mesh
  GEO::CUT::VolumeCell * vc,        ///< volumecell
  const Teuchos::RCP<GEO::CutWizard>& wizard ///< cut wizard
)
{
  bound_f.setf(std::ios::scientific,std::ios::floatfield);
  bound_f.precision(16);

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<2,2> metrictensor;
  double drs;

  std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
  vc->GetBoundaryCells( bcells );
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::iterator i=bcells.begin();
        i!=bcells.end();
        ++i )
  {
    int sid = i->first;
    std::vector<GEO::CUT::BoundaryCell*> & bcs = i->second;

    if(!cond_manager_->IsMeshCoupling(sid)) continue;

    DRT::Element * side = cond_manager_->GetSide( sid );

    GEO::CUT::SideHandle * s = wizard->GetMeshCuttingSide( sid, 0 );

    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    for ( std::vector<GEO::CUT::BoundaryCell*>::iterator i=bcs.begin();
          i!=bcs.end();
          ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      //Issue with boundary cell outputs for marked background sides
      if(bc->GetFacet()->OnMarkedBackgroundSide()) continue;

      switch ( bc->Shape() )
      {
      case DRT::Element::quad4:
        bound_f << "VQ(";
        break;
      case DRT::Element::tri3:
        bound_f << "VT(";
        break;
      default:
//        dserror( "unsupported shape" );
        break;
      }

      const std::vector<GEO::CUT::Point*> & points = bc->Points();
      for ( std::vector<GEO::CUT::Point*>::const_iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        GEO::CUT::Point * p = *i;

        if ( i!=points.begin() )
          bound_f << ",";

        const double * x = p->X();
        bound_f << x[0] << "," << x[1] << "," << x[2];
      }

      bound_f << "){";

      for ( std::vector<GEO::CUT::Point*>::const_iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        GEO::CUT::Point * p = *i;

        // the bc corner points will always lie on the respective side
        const LINALG::Matrix<2,1> & eta = s->LocalCoordinates( p );

        switch ( side->Shape() )
        {
        case DRT::Element::quad4:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad4>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        case DRT::Element::tri3:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::tri3 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::tri3>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        case DRT::Element::quad8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad8>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad8 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad8>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        case DRT::Element::quad9:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad9>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          LINALG::Matrix<2,numnodes> deriv;
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad9 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad9>( xyze, deriv, metrictensor, drs, &normal );
          break;
        }
        default:
        {
          dserror( "unsupported side shape %d", side->Shape() );
          break;
        }
        }

        if ( i!=points.begin() )
          bound_f << ",";

        // side's outward point normal vector (not the bc's normal vector)
        bound_f << normal( 0 ) << "," << normal( 1 ) << "," << normal( 2 );
      }
      bound_f << "};\n";
    }
  }
}

void FLD::XFluidOutputServiceGmsh::GmshOutputDiscretization(
  bool print_faces,
  int step
  )
{
  if (!gmsh_discret_out_)
    return;

  if (discret_->Comm().MyPID() == 0) std::cout << "Discretization output " << discret_->Name() << std::endl;

  // cast to DiscretizationFaces
  Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

  // output for Element and Node IDs
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("DISCRET", step, gmsh_step_diff_, gmsh_debug_out_screen_, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent.precision(16);

  XFEM::UTILS::PrintDiscretizationToStream(discret_,
     discret_->Name(), true, false, true, false, print_faces, false, gmshfilecontent);

  // append other discretizations involved (cutter surface discretization, coupling discretization, etc.)
  cond_manager_->GmshOutputDiscretization(gmshfilecontent);

  gmshfilecontent.close();
}

void FLD::XFluidOutputServiceGmsh::GmshOutputEOS(
  int step,
  Teuchos::RCP<XFEM::XFEM_EdgeStab> edge_stab
)
{
  if(!gmsh_EOS_out_ || edge_stab == Teuchos::null)
    return;

  // cast to DiscretizationXFEM
  Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

  // output for Element and Node IDs
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("EOS", step, gmsh_step_diff_, gmsh_debug_out_screen_, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent.precision(16);

  if( xdiscret->FilledExtension() == true) // stabilization output
  {

    std::map<int,int> & ghost_penalty_map = edge_stab->GetGhostPenaltyMap();

    if (!edge_stab->GetGhostPenaltyMap().empty())
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" " << "ghost penalty stabilized \" {\n";
      for (int i=0; i<xdiscret->NumMyRowFaces(); ++i)
      {
        const DRT::Element* actele = xdiscret->lRowFace(i);


        std::map<int,int>::iterator it = ghost_penalty_map.find(actele->Id());
        if(it != ghost_penalty_map.end())
        {
          int ghost_penalty = it->second;

          if(ghost_penalty) IO::GMSH::elementAtInitialPositionToStream(double(ghost_penalty),actele, gmshfilecontent);
        }
        else dserror("face %d in map not found", actele->Id());
      }
      gmshfilecontent << "};\n";
    }

    std::map<int,int> & edge_based_map = edge_stab->GetEdgeBasedMap();


    if (!edge_stab->GetEdgeBasedMap().empty())
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" " << "edgebased stabilized \" {\n";

      for (int i=0; i<xdiscret->NumMyRowFaces(); ++i)
      {
        const DRT::Element* actele = xdiscret->lRowFace(i);
        std::map<int,int>::iterator it = edge_based_map.find(actele->Id());

        if(it != edge_based_map.end())
        {
          int edge_stab =it->second;

          if(edge_stab) IO::GMSH::elementAtInitialPositionToStream(double(edge_stab),actele, gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    }
  } // end stabilization output

  gmshfilecontent.close();
}
