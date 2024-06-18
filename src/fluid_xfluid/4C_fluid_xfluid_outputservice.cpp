/*----------------------------------------------------------------------*/
/*! \file

\brief Service class for XFluid-related output.

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_xfluid_outputservice.hpp"

#include "4C_cut_cutwizard.hpp"
#include "4C_cut_elementhandle.hpp"
#include "4C_cut_integrationcell.hpp"
#include "4C_cut_sidehandle.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_dofset_transparent_independent.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fluid_xfluid_state.hpp"
#include "4C_fluid_xfluid_state_creator.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_utils_parameter_list.hpp"
#include "4C_xfem_condition_manager.hpp"
#include "4C_xfem_discretization.hpp"
#include "4C_xfem_discretization_utils.hpp"
#include "4C_xfem_edgestab.hpp"

FOUR_C_NAMESPACE_OPEN

FLD::XFluidOutputService::XFluidOutputService(const Teuchos::RCP<XFEM::DiscretizationXFEM>& discret,
    const Teuchos::RCP<XFEM::ConditionManager>& cond_manager)
    : discret_(discret), cond_manager_(cond_manager), firstoutputofrun_(true), restart_count_(0)
{
  // Vector & map extractor for paraview output,
  // mapped to initial fluid dofmap
  dofset_out_ = Teuchos::rcp(new Core::DOFSets::IndependentDofSet());
  velpressplitter_out_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  prepare_output();
}

void FLD::XFluidOutputService::prepare_output()
{
  dofset_out_->Reset();
  dofset_out_->assign_degrees_of_freedom(*discret_, 0, 0);
  const int ndim = Global::Problem::Instance()->NDim();
  // split based on complete fluid field (standard splitter that handles one dofset)
  Core::LinAlg::CreateMapExtractorFromDiscretization(
      *discret_, *dofset_out_, ndim, *velpressplitter_out_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = Core::LinAlg::CreateVector(*dofset_out_->dof_row_map(), true);
}

void FLD::XFluidOutputService::output(int step, double time, bool write_restart_data,
    Teuchos::RCP<const FLD::XFluidState> state, Teuchos::RCP<Epetra_Vector> dispnp,
    Teuchos::RCP<Epetra_Vector> gridvnp)
{
  discret_->Writer()->new_step(step, time);

  // create vector according to the initial row map holding all standard fluid unknowns
  outvec_fluid_->PutScalar(0.0);

  const Epetra_Map* dofrowmap = dofset_out_->dof_row_map();  // original fluid unknowns
  const Epetra_Map* xdofrowmap = discret_->dof_row_map();    // fluid unknown for current cut

  for (int i = 0; i < discret_->NumMyRowNodes(); ++i)
  {
    // get row node via local id
    const Core::Nodes::Node* xfemnode = discret_->lRowNode(i);

    // the initial dofset contains the original dofs for each row node
    const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));

    // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
    // else copy the right nodes
    const std::vector<int> gdofs_current(discret_->Dof(0, xfemnode));

    if (gdofs_current.size() == 0)
    {
      // std::cout << "no dofs available->hole" << std::endl;
    }
    else if (gdofs_current.size() == gdofs_original.size())
    {
      // std::cout << "same number of dofs available" << std::endl;
    }
    else if (gdofs_current.size() > gdofs_original.size())
    {
      // std::cout << "more dofs available->decide" << std::endl;
    }
    else
      std::cout << "decide which dofs can be copied and which have to be set to zero" << std::endl;

    if (gdofs_current.size() == 0)
    {
      size_t numdof = gdofs_original.size();
      // no dofs for this node... must be a hole or somethin'
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        // std::cout << dofrowmap->LID(gdofs[idof]) << std::endl;
        (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
      }
    }
    else if (gdofs_current.size() == gdofs_original.size())
    {
      size_t numdof = gdofs_original.size();
      // copy all values
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        // std::cout << dofrowmap->LID(gdofs[idof]) << std::endl;
        (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] =
            (*state->velnp_)[xdofrowmap->LID(gdofs_current[idof])];
      }
    }
    else if (gdofs_current.size() % gdofs_original.size() == 0)  // multiple dofsets
    {
      // if there are multiple dofsets we write output for the standard dofset
      Core::Geo::Cut::Node* node = state->Wizard()->GetNode(xfemnode->Id());

      const std::vector<Teuchos::RCP<Core::Geo::Cut::NodalDofSet>>& dofcellsets =
          node->NodalDofSets();

      int nds = 0;
      bool is_std_set = false;

      // find the standard dofset
      for (std::vector<Teuchos::RCP<Core::Geo::Cut::NodalDofSet>>::const_iterator cellsets =
               dofcellsets.begin();
           cellsets != dofcellsets.end(); cellsets++)
      {
        is_std_set = (*cellsets)->Is_Standard_DofSet();

        if (is_std_set) break;

        nds++;
      }

      size_t numdof = gdofs_original.size();
      size_t offset =
          0;  // no offset in case of no std-dofset means take the first dofset for output

      if (is_std_set)
      {
        offset = numdof * nds;  // offset to start the loop over dofs at the right nds position
      }

      // copy all values
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] =
            (*state->velnp_)[xdofrowmap->LID(gdofs_current[offset + idof])];
      }


      // if there are just ghost values, then we write output for the set with largest physical
      // fluid volume
    }
    else
      FOUR_C_THROW("unknown number of dofs for output");
  };

  // output (hydrodynamic) pressure for visualization

  Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_out_->ExtractCondVector(outvec_fluid_);

  discret_->Writer()->write_vector("velnp", outvec_fluid_);
  discret_->Writer()->write_vector("pressure", pressure);

  if (dispnp != Teuchos::null)
  {
    if (gridvnp == Teuchos::null) FOUR_C_THROW("Missing grid velocities for ALE-xfluid!");

    // write ale displacement for t^{n+1}
    Teuchos::RCP<Epetra_Vector> dispnprm = Teuchos::rcp(new Epetra_Vector(*dispnp));
    dispnprm->ReplaceMap(outvec_fluid_->Map());  // to get dofs starting by 0 ...
    discret_->Writer()->write_vector("dispnp", dispnprm);

    // write grid velocity for t^{n+1}
    Teuchos::RCP<Epetra_Vector> gridvnprm = Teuchos::rcp(new Epetra_Vector(*gridvnp));
    gridvnprm->ReplaceMap(outvec_fluid_->Map());  // to get dofs starting by 0 ...
    discret_->Writer()->write_vector("gridv", gridvnprm);

    // write convective velocity for t^{n+1}
    Teuchos::RCP<Epetra_Vector> convvel =
        Teuchos::rcp(new Epetra_Vector(outvec_fluid_->Map(), true));
    convvel->Update(1.0, *outvec_fluid_, -1.0, *gridvnprm, 0.0);
    discret_->Writer()->write_vector("convel", convvel);
  }

  discret_->Writer()->write_element_data(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    if (discret_->Comm().MyPID() == 0) Core::IO::cout << "---  write restart... " << Core::IO::endl;

    restart_count_++;

    // velocity/pressure vector
    discret_->Writer()->write_vector("velnp_res", state->velnp_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    discret_->Writer()->write_vector("accnp_res", state->accnp_);
    discret_->Writer()->write_vector("accn_res", state->accn_);
    discret_->Writer()->write_vector("veln_res", state->veln_);
    discret_->Writer()->write_vector("velnm_res", state->velnm_);

    if (dispnp != Teuchos::null)
    {
      // write ale displacement for t^{n+1} on full background
      discret_->Writer()->write_vector("full_dispnp_res", dispnp);

      // write grid velocity for t^{n+1} on full background
      discret_->Writer()->write_vector("full_gridvnp_res", gridvnp);
    }
  }

  //-----------------------------------------------------------
  // REMARK on "Why to clear the MapCache" for restarts
  //-----------------------------------------------------------
  // every time, when output-vectors are written based on a new(!), still unknown map
  // the map is stored in a mapstack (io.cpp, write_vector-routine) for efficiency in standard
  // applications. However, in the XFEM for each timestep or FSI-iteration we have to cut and the
  // map is created newly (done by the fill_complete call). This is the reason why the MapStack
  // increases and the storage is overwritten for large problems. Hence, we have clear the MapCache
  // in regular intervals of written restarts. In case of writing paraview-output, here, we use a
  // standard map which does not change over time, that's okay. For the moment, restart_count = 5 is
  // set quite arbitrary, in case that we need for storage, we have to reduce this number

  if (restart_count_ == 5)
  {
    if (discret_->Comm().MyPID() == 0)
      Core::IO::cout << "\t... Clear MapCache after " << restart_count_ << " written restarts."
                     << Core::IO::endl;

    discret_->Writer()->clear_map_cache();  // clear the output's map-cache
    restart_count_ = 0;
  }

  //-----------------------------------------------------------
  // write paraview output for cutter discretization
  //-----------------------------------------------------------
  cond_manager_->output(step, time, write_restart_data);
}

FLD::XFluidOutputServiceGmsh::XFluidOutputServiceGmsh(Teuchos::ParameterList& params_xfem,
    const Teuchos::RCP<XFEM::DiscretizationXFEM>& discret,
    const Teuchos::RCP<XFEM::ConditionManager>& cond_manager, const bool include_inner)
    : XFluidOutputService(discret, cond_manager),
      gmsh_sol_out_((bool)Core::UTILS::IntegralValue<int>(params_xfem, "GMSH_SOL_OUT")),
      gmsh_ref_sol_out_((bool)Core::UTILS::IntegralValue<int>(params_xfem, "GMSH_TIMINT_OUT")),
      gmsh_debug_out_((bool)Core::UTILS::IntegralValue<int>(params_xfem, "GMSH_DEBUG_OUT")),
      gmsh_debug_out_screen_(
          (bool)Core::UTILS::IntegralValue<int>(params_xfem, "GMSH_DEBUG_OUT_SCREEN")),
      gmsh_eos_out_((bool)Core::UTILS::IntegralValue<int>(params_xfem, "GMSH_EOS_OUT")),
      gmsh_discret_out_((bool)Core::UTILS::IntegralValue<int>(params_xfem, "GMSH_DISCRET_OUT")),
      gmsh_step_diff_(500),
      volume_cell_gauss_point_by_(Core::UTILS::IntegralValue<Inpar::Cut::VCellGaussPts>(
          params_xfem, "VOLUME_GAUSS_POINTS_BY")),
      include_inner_(include_inner)
{
  if (!(bool)Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->IOParams(), "OUTPUT_GMSH"))
    FOUR_C_THROW(
        "If GMSH output is globally deactivated, don't create an instance of "
        "XFluidOutputServiceGmsh!");
};

void FLD::XFluidOutputServiceGmsh::GmshSolutionOutput(
    const std::string& filename_base,             ///< name for output file
    int step,                                     ///< step number
    const Teuchos::RCP<FLD::XFluidState>& state,  ///< state
    int count  ///< counter (no counter for standard solution output : -1)
)
{
  if (!gmsh_sol_out_) return;

  Teuchos::RCP<const Epetra_Vector> output_col_vel =
      Core::Rebalance::GetColVersionOfRowVector(discret_, state->Velnp());

  Teuchos::RCP<const Epetra_Vector> output_col_acc = Teuchos::null;

  if (state->Accnp() != Teuchos::null)
  {
    output_col_acc = Core::Rebalance::GetColVersionOfRowVector(discret_, state->Accnp());
  }

  Teuchos::RCP<const Epetra_Vector> dispnp_col = Teuchos::null;

  if (state->dispnp_ != Teuchos::null)
    dispnp_col = Core::Rebalance::GetColVersionOfRowVector(discret_, state->dispnp_);


  // no counter for standard solution output : -1
  const std::string prefix("SOL");
  GmshOutput(filename_base, prefix, step, count, state->Wizard(), output_col_vel, output_col_acc,
      dispnp_col);
  cond_manager_->GmshOutput(filename_base, step, gmsh_step_diff_, gmsh_debug_out_screen_);
}

void FLD::XFluidOutputServiceGmsh::gmsh_solution_output_previous(
    const std::string& filename_base,             ///< name for output file
    int step,                                     ///< step number
    const Teuchos::RCP<FLD::XFluidState>& state,  ///< state
    int count)
{
  if (!gmsh_ref_sol_out_) return;

  Teuchos::RCP<const Epetra_Vector> output_col_vel =
      Core::Rebalance::GetColVersionOfRowVector(discret_, state->Veln());

  Teuchos::RCP<const Epetra_Vector> output_col_acc = Teuchos::null;

  if (state->Accn() != Teuchos::null)
  {
    output_col_acc = Core::Rebalance::GetColVersionOfRowVector(discret_, state->Accn());
  }

  Teuchos::RCP<const Epetra_Vector> dispnp_col = Teuchos::null;

  if (state->dispnp_ != Teuchos::null)
    dispnp_col = Core::Rebalance::GetColVersionOfRowVector(discret_, state->dispnp_);


  const std::string prefix("ref_SOL");
  GmshOutput(filename_base, prefix, step, count, state->Wizard(), output_col_vel, output_col_acc,
      dispnp_col);
}

void FLD::XFluidOutputServiceGmsh::gmsh_solution_output_debug(
    const std::string& filename_base,  ///< name for output file
    int step,                          ///< step number
    int count,                         ///< counter for iterations within a global time step
    const Teuchos::RCP<FLD::XFluidState>& state  ///< state
)
{
  if (!gmsh_debug_out_) return;

  Teuchos::RCP<const Epetra_Vector> dispnp_col = Teuchos::null;

  if (state->dispnp_ != Teuchos::null)
    dispnp_col = Core::Rebalance::GetColVersionOfRowVector(discret_, state->dispnp_);

  Teuchos::RCP<const Epetra_Vector> output_col_vel =
      Core::Rebalance::GetColVersionOfRowVector(discret_, state->Velnp());
  const std::string prefix("SOL");
  GmshOutput(filename_base, prefix, step, count, state->Wizard(), output_col_vel, Teuchos::null,
      dispnp_col);
}

void FLD::XFluidOutputServiceGmsh::gmsh_residual_output_debug(
    const std::string& filename_base,  ///< name for output file
    int step,                          ///< step number
    int count,                         ///< counter for iterations within a global time step
    const Teuchos::RCP<FLD::XFluidState>& state  ///< state
)
{
  if (!gmsh_debug_out_) return;

  Teuchos::RCP<const Epetra_Vector> dispnp_col = Teuchos::null;

  if (state->dispnp_ != Teuchos::null)
    dispnp_col = Core::Rebalance::GetColVersionOfRowVector(discret_, state->dispnp_);


  Teuchos::RCP<const Epetra_Vector> output_col_residual =
      Core::Rebalance::GetColVersionOfRowVector(discret_, state->Residual());
  const std::string prefix("RES");
  GmshOutput(filename_base, prefix, step, count, state->Wizard(), output_col_residual,
      Teuchos::null, dispnp_col);
}

void FLD::XFluidOutputServiceGmsh::gmsh_increment_output_debug(
    const std::string& filename_base,  ///< name for output file
    int step,                          ///< step number
    int count,                         ///< counter for iterations within a global time step
    const Teuchos::RCP<FLD::XFluidState>& state  ///< state
)
{
  if (!gmsh_debug_out_) return;

  Teuchos::RCP<const Epetra_Vector> dispnp_col = Teuchos::null;

  if (state->dispnp_ != Teuchos::null)
    dispnp_col = Core::Rebalance::GetColVersionOfRowVector(discret_, state->dispnp_);

  Teuchos::RCP<const Epetra_Vector> output_col_incvel =
      Core::Rebalance::GetColVersionOfRowVector(discret_, state->IncVel());
  const std::string prefix("INC");
  GmshOutput(filename_base, prefix, step, count, state->Wizard(), output_col_incvel, Teuchos::null,
      dispnp_col);
}

void FLD::XFluidOutputServiceGmsh::GmshOutput(
    const std::string& filename_base,  ///< name for output file
    const std::string& prefix,         ///< data prefix
    int step,                          ///< step number
    int count,                         ///< counter for iterations within a global time step
    const Teuchos::RCP<Core::Geo::CutWizard>& wizard,  ///< cut wizard
    Teuchos::RCP<const Epetra_Vector> vel,    ///< vector holding velocity and pressure dofs
    Teuchos::RCP<const Epetra_Vector> acc,    ///< vector holding acceleration
    Teuchos::RCP<const Epetra_Vector> dispnp  ///< vector holding acceleration
)
{
  // Todo: should be private
  int myrank = discret_->Comm().MyPID();

  if (myrank == 0) std::cout << "\n\t ... writing Gmsh output...\n" << std::flush;

  bool screen_out = gmsh_debug_out_screen_;

  auto file_name = discret_->Writer()->output()->file_name();

  // output for Element and Node IDs
  std::ostringstream filename_base_vel;
  if (count > -1)
    filename_base_vel << filename_base << "_" << count << "_vel";
  else
    filename_base_vel << filename_base << "_vel";
  const std::string filename_vel = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_vel.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_vel(filename_vel.c_str());
  gmshfilecontent_vel.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_vel.precision(16);

  std::ostringstream filename_base_press;
  if (count > -1)
    filename_base_press << filename_base << "_" << count << "_press";
  else
    filename_base_press << filename_base << "_press";
  const std::string filename_press = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_press.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_press(filename_press.c_str());
  gmshfilecontent_press.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_press.precision(16);

  std::ostringstream filename_base_acc;
  if (count > -1)
    filename_base_acc << filename_base << "_" << count << "_acc";
  else
    filename_base_acc << filename_base << "_acc";
  const std::string filename_acc = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_acc.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_acc(filename_acc.c_str());
  gmshfilecontent_acc.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_acc.precision(16);

  std::ostringstream filename_base_bound;
  if (count > -1)
    filename_base_bound << filename_base << "_" << count << "_bound";
  else
    filename_base_bound << filename_base << "_bound";
  const std::string filename_bound = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_bound.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_bound(filename_bound.c_str());
  gmshfilecontent_bound.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_bound.precision(16);


  // output for Element and Node IDs
  std::ostringstream filename_base_vel_ghost;
  if (count > -1)
    filename_base_vel_ghost << filename_base << "_" << count << "_vel_ghost";
  else
    filename_base_vel_ghost << filename_base << "_vel_ghost";
  const std::string filename_vel_ghost = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_vel_ghost.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_vel_ghost(filename_vel_ghost.c_str());
  gmshfilecontent_vel_ghost.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_vel_ghost.precision(16);

  std::ostringstream filename_base_press_ghost;
  if (count > -1)
    filename_base_press_ghost << filename_base << "_" << count << "_press_ghost";
  else
    filename_base_press_ghost << filename_base << "_press_ghost";
  const std::string filename_press_ghost = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_press_ghost.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_press_ghost(filename_press_ghost.c_str());
  gmshfilecontent_press_ghost.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_press_ghost.precision(16);

  std::ostringstream filename_base_acc_ghost;
  if (count > -1)
    filename_base_acc_ghost << filename_base << "_" << count << "_acc_ghost";
  else
    filename_base_acc_ghost << filename_base << "_acc_ghost";
  const std::string filename_acc_ghost = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(
      filename_base_acc_ghost.str(), file_name, step, gmsh_step_diff_, screen_out, myrank);
  if (gmsh_debug_out_screen_) std::cout << std::endl;
  std::ofstream gmshfilecontent_acc_ghost(filename_acc_ghost.c_str());
  gmshfilecontent_acc_ghost.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_acc_ghost.precision(16);


  if (count > -1)  // for residual output
  {
    gmshfilecontent_vel << "View \"" << prefix << "vel " << count << "\" {\n";
    gmshfilecontent_press << "View \"" << prefix << "press " << count << "\" {\n";
    gmshfilecontent_bound << "View \"" << prefix << "side-normal " << count << "\" {\n";
    gmshfilecontent_vel_ghost << "View \"" << prefix << "vel_ghost " << count << "\" {\n";
    gmshfilecontent_press_ghost << "View \"" << prefix << "press_ghost " << count << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel << "View \"" << prefix << "vel "
                        << "\" {\n";
    gmshfilecontent_press << "View \"" << prefix << "press "
                          << "\" {\n";
    gmshfilecontent_acc << "View \"" << prefix << "acc  "
                        << "\" {\n";
    gmshfilecontent_bound << "View \"" << prefix << "side-normal "
                          << "\" {\n";
    gmshfilecontent_vel_ghost << "View \"" << prefix << "vel_ghost "
                              << "\" {\n";
    gmshfilecontent_press_ghost << "View \"" << prefix << "press_ghost "
                                << "\" {\n";
    gmshfilecontent_acc_ghost << "View \"" << prefix << "acc_ghost  "
                              << "\" {\n";
  }

  const int numrowele = (discret_)->NumMyRowElements();
  for (int i = 0; i < numrowele; ++i)
  {
    Core::Elements::Element* actele = (discret_)->lRowElement(i);
    Core::Geo::Cut::ElementHandle* e = wizard->GetElement(actele);

    if (e != nullptr)
    {
      std::vector<Core::Geo::Cut::plain_volumecell_set> cell_sets;
      std::vector<std::vector<int>> nds_sets;

      e->get_volume_cells_dof_sets(cell_sets, nds_sets, include_inner_);


      if (e->IsIntersected())
      {
        int set_counter = 0;

        for (std::vector<Core::Geo::Cut::plain_volumecell_set>::iterator s = cell_sets.begin();
             s != cell_sets.end(); s++)
        {
          Core::Geo::Cut::plain_volumecell_set& cells = *s;

          std::vector<int>& nds = nds_sets[set_counter];


          for (Core::Geo::Cut::plain_volumecell_set::iterator i = cells.begin(); i != cells.end();
               ++i)
          {
            Core::Geo::Cut::VolumeCell* vc = *i;

            if (e->IsCut())
            {
              gmsh_output_volume_cell(*discret_, gmshfilecontent_vel, gmshfilecontent_press,
                  gmshfilecontent_acc, actele, e, vc, nds, vel, acc);
              if (vc->Position() == Core::Geo::Cut::Point::outside)
              {
                if (cond_manager_->HasMeshCoupling())
                  gmsh_output_boundary_cell(*discret_, gmshfilecontent_bound, vc, wizard);
              }
            }
            else
            {
              gmsh_output_element(*discret_, gmshfilecontent_vel, gmshfilecontent_press,
                  gmshfilecontent_acc, actele, nds, vel, acc, dispnp);
            }
            gmsh_output_element(*discret_, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost,
                gmshfilecontent_acc_ghost, actele, nds, vel, acc, dispnp);
          }
          set_counter += 1;
        }
      }
      else if (cell_sets.size() == 1)  // non intersected but one set
      {
        std::vector<int>& nds = nds_sets[0];

        // one standard uncut physical element
        gmsh_output_element(*discret_, gmshfilecontent_vel, gmshfilecontent_press,
            gmshfilecontent_acc, actele, nds, vel, acc, dispnp);
      }
      else
      {
        // ghost element
      }
    }     // element handle
    else  // no element handle
    {
      std::vector<int> nds;  // empty vector
      gmsh_output_element(*discret_, gmshfilecontent_vel, gmshfilecontent_press,
          gmshfilecontent_acc, actele, nds, vel, acc, dispnp);
      gmsh_output_element(*discret_, gmshfilecontent_vel_ghost, gmshfilecontent_press_ghost,
          gmshfilecontent_acc_ghost, actele, nds, vel, acc, dispnp);
    }
  }  // loop elements

  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";
  if (count == -1) gmshfilecontent_acc << "};\n";
  gmshfilecontent_vel_ghost << "};\n";
  gmshfilecontent_press_ghost << "};\n";
  if (count == -1) gmshfilecontent_acc_ghost << "};\n";
  gmshfilecontent_bound << "};\n";

  gmshfilecontent_vel.close();
  gmshfilecontent_press.close();
  gmshfilecontent_acc.close();
  gmshfilecontent_bound.close();
  gmshfilecontent_vel_ghost.close();
  gmshfilecontent_press_ghost.close();
  gmshfilecontent_acc_ghost.close();

  if (myrank == 0) std::cout << " done\n" << std::flush;
}

/// Gmsh output function for elements without an Core::Geo::Cut::ElementHandle
void FLD::XFluidOutputServiceGmsh::gmsh_output_element(
    Core::FE::Discretization& discret,        ///< background fluid discretization
    std::ofstream& vel_f,                     ///< output file stream for velocity
    std::ofstream& press_f,                   ///< output file stream for pressure
    std::ofstream& acc_f,                     ///< output file stream for acceleration
    Core::Elements::Element* actele,          ///< element
    std::vector<int>& nds,                    ///< vector holding the nodal dofsets
    Teuchos::RCP<const Epetra_Vector> vel,    ///< vector holding velocity and pressure dofs
    Teuchos::RCP<const Epetra_Vector> acc,    ///< vector holding acceleration
    Teuchos::RCP<const Epetra_Vector> dispnp  ///< vector holding ale displacements
)
{
  vel_f.setf(std::ios::scientific, std::ios::floatfield);
  vel_f.precision(16);

  press_f.setf(std::ios::scientific, std::ios::floatfield);
  press_f.precision(16);

  acc_f.setf(std::ios::scientific, std::ios::floatfield);
  acc_f.precision(16);

  // output for accvec ?
  const bool acc_output(acc != Teuchos::null);

  Core::Elements::Element::LocationArray la(1);


  if (nds.size() != 0)  // for element output of ghost values
  {
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(discret, nds, la, false);
  }
  else
  {
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(discret, la, false);
  }

  std::vector<double> m(la[0].lm_.size());
  Core::FE::ExtractMyValues(*vel, m, la[0].lm_);

  std::vector<double> m_acc(la[0].lm_.size());
  if (acc_output)
  {
    Core::FE::ExtractMyValues(*acc, m_acc, la[0].lm_);
  }

  const bool ale_output(dispnp != Teuchos::null);

  std::vector<double> m_disp(la[0].lm_.size());
  if (ale_output)
  {
    Core::FE::ExtractMyValues(*dispnp, m_disp, la[0].lm_);
  }

  int numnode = 0;

  switch (actele->Shape())
  {
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
      numnode = 8;
      vel_f << "VH(";
      press_f << "SH(";
      if (acc_output) acc_f << "VH(";
      break;
    case Core::FE::CellType::wedge6:
    case Core::FE::CellType::wedge15:
      numnode = 6;
      vel_f << "VI(";
      press_f << "SI(";
      if (acc_output) acc_f << "VI(";
      break;
    case Core::FE::CellType::pyramid5:
      numnode = 5;
      vel_f << "VY(";
      press_f << "SY(";
      if (acc_output) acc_f << "VY(";
      break;
    case Core::FE::CellType::tet4:
    case Core::FE::CellType::tet10:
      numnode = 4;
      vel_f << "VS(";
      press_f << "SS(";
      if (acc_output) acc_f << "VS(";
      break;
    default:
    {
      FOUR_C_THROW("unsupported shape");
      break;
    }
  }

  for (int i = 0; i < numnode; ++i)
  {
    if (i > 0)
    {
      vel_f << ",";
      press_f << ",";
      if (acc_output) acc_f << ",";
    }
    int j = 4 * i;
    const auto& x = actele->Nodes()[i]->X();
    vel_f << x[0] + m_disp[j + 0] << "," << x[1] + m_disp[j + 1] << "," << x[2] + m_disp[j + 2];
    press_f << x[0] + m_disp[j + 0] << "," << x[1] + m_disp[j + 1] << "," << x[2] + m_disp[j + 2];
    if (acc_output)
      acc_f << x[0] + m_disp[j + 0] << "," << x[1] + m_disp[j + 1] << "," << x[2] + m_disp[j + 2];
  }
  vel_f << "){";
  press_f << "){";
  if (acc_output) acc_f << "){";

  for (int i = 0; i < numnode; ++i)
  {
    if (i > 0)
    {
      vel_f << ",";
      press_f << ",";
      if (acc_output) acc_f << ",";
    }
    int j = 4 * i;
    vel_f << m[j] << "," << m[j + 1] << "," << m[j + 2];
    press_f << m[j + 3];
    if (acc_output) acc_f << m_acc[j] << "," << m_acc[j + 1] << "," << m_acc[j + 2];
  }

  vel_f << "};\n";
  press_f << "};\n";
  if (acc_output) acc_f << "};\n";
}

/// Gmsh output function for volumecells
void FLD::XFluidOutputServiceGmsh::gmsh_output_volume_cell(
    Core::FE::Discretization& discret,         ///< background fluid discretization
    std::ofstream& vel_f,                      ///< output file stream for velocity
    std::ofstream& press_f,                    ///< output file stream for pressure
    std::ofstream& acc_f,                      ///< output file stream for acceleration
    Core::Elements::Element* actele,           ///< element
    Core::Geo::Cut::ElementHandle* e,          ///< elementhandle
    Core::Geo::Cut::VolumeCell* vc,            ///< volumecell
    const std::vector<int>& nds,               ///< vector holding the nodal dofsets
    Teuchos::RCP<const Epetra_Vector> velvec,  ///< vector holding velocity and pressure dofs
    Teuchos::RCP<const Epetra_Vector> accvec   ///< vector holding acceleration
)
{
  vel_f.setf(std::ios::scientific, std::ios::floatfield);
  vel_f.precision(16);

  press_f.setf(std::ios::scientific, std::ios::floatfield);
  press_f.precision(16);

  acc_f.setf(std::ios::scientific, std::ios::floatfield);
  acc_f.precision(16);



  // output for accvec ?
  bool acc_output = true;
  if (accvec == Teuchos::null) acc_output = false;

  Core::Elements::Element::LocationArray la(1);

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret, nds, la, false);

  std::vector<double> m(la[0].lm_.size());
  Core::FE::ExtractMyValues(*velvec, m, la[0].lm_);

  std::vector<double> m_acc(la[0].lm_.size());
  if (acc_output)
  {
    Core::FE::ExtractMyValues(*accvec, m_acc, la[0].lm_);
  }

  Core::LinAlg::SerialDenseMatrix vel(3, actele->num_node());
  Core::LinAlg::SerialDenseMatrix press(1, actele->num_node());
  Core::LinAlg::SerialDenseMatrix acc(3, actele->num_node());

  for (int i = 0; i < actele->num_node(); ++i)
  {
    vel(0, i) = m[4 * i + 0];
    vel(1, i) = m[4 * i + 1];
    vel(2, i) = m[4 * i + 2];

    press(0, i) = m[4 * i + 3];

    if (acc_output)
    {
      acc(0, i) = m_acc[4 * i + 0];
      acc(1, i) = m_acc[4 * i + 1];
      acc(2, i) = m_acc[4 * i + 2];
    }
  }

  // facet based output for cut volumes
  // integrationcells are not available because tessellation is not used
  if (volume_cell_gauss_point_by_ != Inpar::Cut::VCellGaussPts_Tessellation)
  {
    const Core::Geo::Cut::plain_facet_set& facete = vc->Facets();
    for (Core::Geo::Cut::plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
    {
      // split facet into tri and quad cell
      Core::Geo::Cut::Facet* fe = *i;
      std::vector<std::vector<Core::Geo::Cut::Point*>> split;
      std::vector<Core::Geo::Cut::Point*> corners = fe->CornerPoints();

      if (corners.size() == 3)  // only Tri can be used directly. Quad may be concave
        split.push_back(corners);
      else
      {
        if (!fe->IsFacetSplit()) fe->SplitFacet(fe->CornerPoints());
        split = fe->GetSplitCells();
      }

      for (std::vector<std::vector<Core::Geo::Cut::Point*>>::const_iterator j = split.begin();
           j != split.end(); j++)
      {
        std::vector<Core::Geo::Cut::Point*> cell = *j;

        switch (cell.size())
        {
          case 3:
            vel_f << "VT(";
            press_f << "ST(";
            if (acc_output) acc_f << "VT(";
            break;
          case 4:
            vel_f << "VQ(";
            press_f << "SQ(";
            if (acc_output) acc_f << "VQ(";
            break;
          default:
          {
            FOUR_C_THROW("splitting facets failed");
            break;
          }
        }

        for (unsigned k = 0; k < cell.size(); ++k)
        {
          if (k > 0)
          {
            vel_f << ",";
            press_f << ",";
            if (acc_output) acc_f << ",";
          }
          const double* x = cell[k]->X();
          vel_f << x[0] << "," << x[1] << "," << x[2];
          press_f << x[0] << "," << x[1] << "," << x[2];
          if (acc_output) acc_f << x[0] << "," << x[1] << "," << x[2];
        }
        vel_f << "){";
        press_f << "){";
        if (acc_output) acc_f << "){";

        for (unsigned k = 0; k < cell.size(); ++k)
        {
          Core::LinAlg::Matrix<3, 1> v(true);
          Core::LinAlg::Matrix<1, 1> p(true);
          Core::LinAlg::Matrix<3, 1> a(true);

          Core::Geo::Cut::Point* point = cell[k];
          const Core::LinAlg::Matrix<3, 1>& rst = e->local_coordinates(point);

          switch (actele->Shape())
          {
            case Core::FE::CellType::hex8:
            {
              const int numnodes = Core::FE::num_nodes<Core::FE::CellType::hex8>;
              Core::LinAlg::Matrix<numnodes, 1> funct;
              Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::hex8);
              Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
              Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
              Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

              v.Multiply(1, velocity, funct, 1);
              p.Multiply(1, pressure, funct, 1);
              if (acc_output) a.Multiply(1, acceleration, funct, 1);
              break;
            }
            case Core::FE::CellType::hex20:
            {
              // TODO: check the output for hex20
              const int numnodes = Core::FE::num_nodes<Core::FE::CellType::hex20>;
              Core::LinAlg::Matrix<numnodes, 1> funct;
              Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::hex20);
              Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
              Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
              Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

              v.Multiply(1, velocity, funct, 1);
              p.Multiply(1, pressure, funct, 1);
              if (acc_output) a.Multiply(1, acceleration, funct, 1);
              break;
            }
            case Core::FE::CellType::hex27:
            {
              // TODO: check the output for hex27
              const int numnodes = Core::FE::num_nodes<Core::FE::CellType::hex27>;
              Core::LinAlg::Matrix<numnodes, 1> funct;
              Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::hex27);
              Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
              Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
              Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

              v.Multiply(1, velocity, funct, 1);
              p.Multiply(1, pressure, funct, 1);
              if (acc_output) a.Multiply(1, acceleration, funct, 1);
              break;
            }
            default:
            {
              FOUR_C_THROW("unsupported shape");
              break;
            }
          }

          if (k > 0)
          {
            vel_f << ",";
            press_f << ",";
            if (acc_output) acc_f << ",";
          }
          vel_f << v(0) << "," << v(1) << "," << v(2);
          press_f << p(0);
          if (acc_output) acc_f << a(0) << "," << a(1) << "," << a(2);
        }

        vel_f << "};\n";
        press_f << "};\n";
        if (acc_output) acc_f << "};\n";
      }
    }
  }

  // integrationcells based output for tessellation
  else
  {
    const Core::Geo::Cut::plain_integrationcell_set& intcells = vc->IntegrationCells();
    for (Core::Geo::Cut::plain_integrationcell_set::const_iterator i = intcells.begin();
         i != intcells.end(); ++i)
    {
      Core::Geo::Cut::IntegrationCell* ic = *i;

      const std::vector<Core::Geo::Cut::Point*>& points = ic->Points();

      switch (ic->Shape())
      {
        case Core::FE::CellType::hex8:
          vel_f << "VH(";
          press_f << "SH(";
          if (acc_output) acc_f << "VH(";
          break;
        case Core::FE::CellType::tet4:
          vel_f << "VS(";
          press_f << "SS(";
          if (acc_output) acc_f << "VS(";
          break;
        default:
        {
          FOUR_C_THROW("unsupported shape");
          break;
        }
      }

      for (unsigned i = 0; i < points.size(); ++i)
      {
        if (i > 0)
        {
          vel_f << ",";
          press_f << ",";
          if (acc_output) acc_f << ",";
        }
        const double* x = points[i]->X();
        vel_f << x[0] << "," << x[1] << "," << x[2];
        press_f << x[0] << "," << x[1] << "," << x[2];
        if (acc_output) acc_f << x[0] << "," << x[1] << "," << x[2];
      }
      vel_f << "){";
      press_f << "){";
      if (acc_output) acc_f << "){";

      for (unsigned i = 0; i < points.size(); ++i)
      {
        Core::LinAlg::Matrix<3, 1> v(true);
        Core::LinAlg::Matrix<1, 1> p(true);
        Core::LinAlg::Matrix<3, 1> a(true);

        Core::Geo::Cut::Point* point = points[i];
        const Core::LinAlg::Matrix<3, 1>& rst = e->local_coordinates(point);

        switch (actele->Shape())
        {
          case Core::FE::CellType::hex8:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::hex8>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::hex8);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            // TODO: check the output for hex20
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::hex20>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::hex20);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            // TODO: check the output for hex27
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::hex27>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::hex27);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          case Core::FE::CellType::wedge6:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::wedge6>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::wedge6);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          case Core::FE::CellType::wedge15:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::wedge15>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::wedge15);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::tet4>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::tet4);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::tet10>;
            Core::LinAlg::Matrix<numnodes, 1> funct;
            Core::FE::shape_function_3D(funct, rst(0), rst(1), rst(2), Core::FE::CellType::tet10);
            Core::LinAlg::Matrix<3, numnodes> velocity(vel, true);
            Core::LinAlg::Matrix<1, numnodes> pressure(press, true);
            Core::LinAlg::Matrix<3, numnodes> acceleration(acc, true);

            v.Multiply(1, velocity, funct, 1);
            p.Multiply(1, pressure, funct, 1);
            if (acc_output) a.Multiply(1, acceleration, funct, 1);
            break;
          }
          default:
          {
            FOUR_C_THROW("unsupported shape");
            break;
          }
        }


        if (i > 0)
        {
          vel_f << ",";
          press_f << ",";
          if (acc_output) acc_f << ",";
        }
        vel_f << v(0) << "," << v(1) << "," << v(2);
        press_f << p(0);
        if (acc_output) acc_f << a(0) << "," << a(1) << "," << a(2);
      }

      vel_f << "};\n";
      press_f << "};\n";
      if (acc_output) acc_f << "};\n";
    }
  }
}

/// Gmsh output function for boundarycells
void FLD::XFluidOutputServiceGmsh::gmsh_output_boundary_cell(
    Core::FE::Discretization& discret,                ///< background fluid discretization
    std::ofstream& bound_f,                           ///< output file stream for boundary mesh
    Core::Geo::Cut::VolumeCell* vc,                   ///< volumecell
    const Teuchos::RCP<Core::Geo::CutWizard>& wizard  ///< cut wizard
)
{
  bound_f.setf(std::ios::scientific, std::ios::floatfield);
  bound_f.precision(16);

  Core::LinAlg::Matrix<3, 1> normal;
  Core::LinAlg::Matrix<2, 2> metrictensor;
  double drs;

  std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>> bcells;
  vc->GetBoundaryCells(bcells);
  for (std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>::iterator i = bcells.begin();
       i != bcells.end(); ++i)
  {
    int sid = i->first;
    std::vector<Core::Geo::Cut::BoundaryCell*>& bcs = i->second;

    if (!cond_manager_->IsMeshCoupling(sid)) continue;

    Core::Elements::Element* side = cond_manager_->get_side(sid);

    Core::Geo::Cut::SideHandle* s = wizard->GetMeshCuttingSide(sid, 0);

    const int numnodes = side->num_node();
    Core::Nodes::Node** nodes = side->Nodes();
    Core::LinAlg::SerialDenseMatrix side_xyze(3, numnodes);
    for (int i = 0; i < numnodes; ++i)
    {
      const double* x = nodes[i]->X().data();
      std::copy(x, x + 3, &side_xyze(0, i));
    }

    for (std::vector<Core::Geo::Cut::BoundaryCell*>::iterator i = bcs.begin(); i != bcs.end(); ++i)
    {
      Core::Geo::Cut::BoundaryCell* bc = *i;

      // Issue with boundary cell outputs for marked background sides
      if (bc->GetFacet()->on_marked_background_side()) continue;

      switch (bc->Shape())
      {
        case Core::FE::CellType::quad4:
          bound_f << "VQ(";
          break;
        case Core::FE::CellType::tri3:
          bound_f << "VT(";
          break;
        default:
          //        FOUR_C_THROW( "unsupported shape" );
          break;
      }

      const std::vector<Core::Geo::Cut::Point*>& points = bc->Points();
      for (std::vector<Core::Geo::Cut::Point*>::const_iterator i = points.begin();
           i != points.end(); ++i)
      {
        Core::Geo::Cut::Point* p = *i;

        if (i != points.begin()) bound_f << ",";

        const double* x = p->X();
        bound_f << x[0] << "," << x[1] << "," << x[2];
      }

      bound_f << "){";

      for (std::vector<Core::Geo::Cut::Point*>::const_iterator i = points.begin();
           i != points.end(); ++i)
      {
        Core::Geo::Cut::Point* p = *i;

        // the bc corner points will always lie on the respective side
        const Core::LinAlg::Matrix<2, 1>& eta = s->local_coordinates(p);

        switch (side->Shape())
        {
          case Core::FE::CellType::quad4:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::quad4>;
            Core::LinAlg::Matrix<3, numnodes> xyze(side_xyze, true);
            Core::LinAlg::Matrix<2, numnodes> deriv;
            Core::FE::shape_function_2D_deriv1(deriv, eta(0), eta(1), Core::FE::CellType::quad4);
            Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::quad4>(
                xyze, deriv, metrictensor, drs, &normal);
            break;
          }
          case Core::FE::CellType::tri3:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::tri3>;
            Core::LinAlg::Matrix<3, numnodes> xyze(side_xyze, true);
            Core::LinAlg::Matrix<2, numnodes> deriv;
            Core::FE::shape_function_2D_deriv1(deriv, eta(0), eta(1), Core::FE::CellType::tri3);
            Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::tri3>(
                xyze, deriv, metrictensor, drs, &normal);
            break;
          }
          case Core::FE::CellType::quad8:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::quad8>;
            Core::LinAlg::Matrix<3, numnodes> xyze(side_xyze, true);
            Core::LinAlg::Matrix<2, numnodes> deriv;
            Core::FE::shape_function_2D_deriv1(deriv, eta(0), eta(1), Core::FE::CellType::quad8);
            Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::quad8>(
                xyze, deriv, metrictensor, drs, &normal);
            break;
          }
          case Core::FE::CellType::quad9:
          {
            const int numnodes = Core::FE::num_nodes<Core::FE::CellType::quad9>;
            Core::LinAlg::Matrix<3, numnodes> xyze(side_xyze, true);
            Core::LinAlg::Matrix<2, numnodes> deriv;
            Core::FE::shape_function_2D_deriv1(deriv, eta(0), eta(1), Core::FE::CellType::quad9);
            Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::quad9>(
                xyze, deriv, metrictensor, drs, &normal);
            break;
          }
          default:
          {
            FOUR_C_THROW("unsupported side shape %d", side->Shape());
            break;
          }
        }

        if (i != points.begin()) bound_f << ",";

        // side's outward point normal vector (not the bc's normal vector)
        bound_f << normal(0) << "," << normal(1) << "," << normal(2);
      }
      bound_f << "};\n";
    }
  }
}

void FLD::XFluidOutputServiceGmsh::gmsh_output_discretization(
    bool print_faces, int step, std::map<int, Core::LinAlg::Matrix<3, 1>>* curr_pos)
{
  if (!gmsh_discret_out_) return;

  if (discret_->Comm().MyPID() == 0)
    std::cout << "discretization output " << discret_->Name() << std::endl;

  // cast to DiscretizationFaces
  Teuchos::RCP<Core::FE::DiscretizationFaces> xdiscret =
      Teuchos::rcp_dynamic_cast<Core::FE::DiscretizationFaces>(discret_, true);
  if (xdiscret == Teuchos::null)
    FOUR_C_THROW(
        "Failed to cast Core::FE::Discretization to "
        "Core::FE::DiscretizationFaces.");

  // output for Element and Node IDs
  const std::string filename = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles("DISCRET",
      discret_->Writer()->output()->file_name(), step, gmsh_step_diff_, gmsh_debug_out_screen_,
      discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  gmshfilecontent.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent.precision(16);

  XFEM::UTILS::PrintDiscretizationToStream(discret_, discret_->Name(), true, false, true, false,
      print_faces, false, gmshfilecontent, curr_pos);

  // append other discretizations involved (cutter surface discretization, coupling discretization,
  // etc.)
  cond_manager_->gmsh_output_discretization(gmshfilecontent);

  gmshfilecontent.close();
}

void FLD::XFluidOutputServiceGmsh::GmshOutputEOS(
    int step, Teuchos::RCP<XFEM::XfemEdgeStab> edge_stab)
{
  if (!gmsh_eos_out_ || edge_stab == Teuchos::null) return;

  // cast to DiscretizationXFEM
  Teuchos::RCP<Core::FE::DiscretizationFaces> xdiscret =
      Teuchos::rcp_dynamic_cast<Core::FE::DiscretizationFaces>(discret_, true);
  if (xdiscret == Teuchos::null)
    FOUR_C_THROW(
        "Failed to cast Core::FE::Discretization to "
        "Core::FE::DiscretizationFaces.");

  // output for Element and Node IDs
  const std::string filename = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles("EOS",
      discret_->Writer()->output()->file_name(), step, gmsh_step_diff_, gmsh_debug_out_screen_,
      discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  gmshfilecontent.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent.precision(16);

  if (xdiscret->FilledExtension() == true)  // stabilization output
  {
    std::map<int, int>& ghost_penalty_map = edge_stab->GetGhostPenaltyMap();

    if (!edge_stab->GetGhostPenaltyMap().empty())
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" "
                      << "ghost penalty stabilized \" {\n";
      for (int i = 0; i < xdiscret->NumMyRowFaces(); ++i)
      {
        const Core::Elements::Element* actele = xdiscret->lRowFace(i);


        std::map<int, int>::iterator it = ghost_penalty_map.find(actele->Id());
        if (it != ghost_penalty_map.end())
        {
          int ghost_penalty = it->second;

          if (ghost_penalty)
            Core::IO::Gmsh::elementAtInitialPositionToStream(
                double(ghost_penalty), actele, gmshfilecontent);
        }
        else
          FOUR_C_THROW("face %d in map not found", actele->Id());
      }
      gmshfilecontent << "};\n";
    }

    std::map<int, int>& edge_based_map = edge_stab->GetEdgeBasedMap();


    if (!edge_stab->GetEdgeBasedMap().empty())
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" "
                      << "edgebased stabilized \" {\n";

      for (int i = 0; i < xdiscret->NumMyRowFaces(); ++i)
      {
        const Core::Elements::Element* actele = xdiscret->lRowFace(i);
        std::map<int, int>::iterator it = edge_based_map.find(actele->Id());

        if (it != edge_based_map.end())
        {
          int edge_stab = it->second;

          if (edge_stab)
            Core::IO::Gmsh::elementAtInitialPositionToStream(
                double(edge_stab), actele, gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    }
  }  // end stabilization output

  gmshfilecontent.close();
}

FOUR_C_NAMESPACE_CLOSE
