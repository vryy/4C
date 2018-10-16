/*----------------------------------------------------------------------*/
/*!
\file xfem_coupling_fpi_mesh.cpp

\brief manages mesh based coupling of fluid and porous media
xfluid class and the cut-library

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "xfem_coupling_fpi_mesh.H"

#include "xfem_utils.H"
#include "xfem_interface_utils.H"
#include "xfem_discretization_utils.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/fluidporo.H"

//! constructor
XFEM::MeshCouplingFPI::MeshCouplingFPI(
    Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<DRT::Discretization>&
        cond_dis,           ///< discretization from which cutter discretization can be derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step,         ///< time step
    MeshCouplingFPI::coupled_field field  ///< which field is coupled to the fluid
    )
    : MeshVolCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step,
          (field == MeshCouplingFPI::ps_ps
                  ? "_ps_ps"
                  : (field == MeshCouplingFPI::ps_pf
                            ? "_ps_pf"
                            : (field == MeshCouplingFPI::pf_ps ? "_pf_ps" : "_pf_pf")))),
      coupled_field_(field)
{
  // TODO: init here, but set in SetConditionSpecificParameters
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::SetCouplingName()
{
  // Set couling name dependent on the field type ... to be able to differentiate between the two
  // fpi coupling objects!
  std::stringstream str;
  str << cond_name_
      << (coupled_field_ == MeshCouplingFPI::ps_ps
                 ? "_ps_ps"
                 : (coupled_field_ == MeshCouplingFPI::ps_pf
                           ? "_ps_pf"
                           : (coupled_field_ == MeshCouplingFPI::pf_ps ? "_pf_ps" : "_pf_pf")));

  // replace the set condname by its specification given by the coupling field
  coupl_name_ = str.str();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::InitStateVectors()
{
  XFEM::MeshCoupling::InitStateVectors();

  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();
  const Epetra_Map* cutterdofcolmap = cutter_dis_->DofColMap();

  itrueresidual_ = LINALG::CreateVector(*cutterdofrowmap, true);
  iforcecol_ = LINALG::CreateVector(*cutterdofcolmap, true);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::DoConditionSpecificSetup()
{
  // required for parallel computing - deactived due to linker problem for post processor & cut_test
  // ager 12/16
#if (0)
  // XFEM::MeshCoupling::DoConditionSpecificSetup();

  // We do ghosting just for the first fpi coupling object, which is PS_PS, then this is done an we
  // just reconnect to the parent pointers for every cutter_dis_ Actually this reconnecting step
  // should be done in an setup routine, to guarantee, that it is done late enought (no ghosting
  // afterwards...)
  if (coupled_field_ == MeshCouplingFPI::ps_ps)
    POROELAST::UTILS::CreateVolumeGhosting(*cutter_dis_);
  else
    POROELAST::UTILS::ReconnectParentPointers(*cutter_dis_, *cond_dis_);

#endif
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::CompleteStateVectors()
{
  //-------------------------------------------------------------------------------
  // finalize itrueresidual vector

  // need to export the interface forces
  Epetra_Vector iforce_tmp(itrueresidual_->Map(), true);
  Epetra_Export exporter_iforce(iforcecol_->Map(), iforce_tmp.Map());
  int err1 = iforce_tmp.Export(*iforcecol_, exporter_iforce, Add);
  if (err1) dserror("Export using exporter returned err=%d", err1);

  // scale the interface trueresidual with -1.0 to get the forces acting on structural side (no
  // residual-scaling!)
  itrueresidual_->Update(-1.0, iforce_tmp, 0.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::SetupConfigurationMap()
{
  switch (coupled_field_)
  {
    case MeshCouplingFPI::ps_ps:
    {
      // Configuration of Consistency Terms

      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool, double>(true, 1.0);
      }
      else
      {
        configuration_map_[INPAR::XFEM::F_Con_n_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Con_n_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool, double>(true, 1.0);
      }
      // configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0-porosity);
      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);

      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] =
            std::pair<bool, double>(true, 1.0);  // Here we need alpha BJ finally!!! Todo
      }

      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool, double>(true, 1.0);

      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Con_t_Row] =
            std::pair<bool, double>(true, 1.0);  //+sign for penalty!
        configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool, double>(true, 1.0);
        configuration_map_[INPAR::XFEM::X_Con_t_Row] =
            std::pair<bool, double>(true, 1.0);  //+sign for penalty!
      }

#ifdef INFLOW_STAB
      configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF1] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF2] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_Row_linF3] = std::pair<bool, double>(true, 1.0);
#endif
      break;
    }
    case MeshCouplingFPI::ps_pf:
    {
      // Configuration of Adjount Consistency Terms
      configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool, double>(true, 1.0);

      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool, double>(true, 1.0);
        if (full_BJ_)
          configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool, double>(true, 1.0);
      }

      //    //Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool, double>(true, 1.0);
      if (full_BJ_)
        configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool, double>(true, 1.0);
      break;
    }
    case MeshCouplingFPI::pf_ps:
    {
      // Configuration of Consistency Terms
      configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Con_n_Col] = std::pair<bool, double>(true, 1.0);
      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      break;
    }
    case MeshCouplingFPI::pf_pf:
    {
      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool, double>(true, 1.0);
      configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool, double>(true, 1.0);
      break;
    }
  }
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
    double& visc_m,          //< master sided dynamic viscosity
    double& visc_s,          //< slave sided dynamic viscosity
    double& density_m,       //< master sided density
    double& visc_stab_tang,  //< viscous tangential NIT Penalty scaling
    double& full_stab, const LINALG::Matrix<3, 1>& x, const DRT::Condition* cond,
    DRT::Element* ele,   //< Element
    DRT::Element* bele,  //< Boundary Element
    double* funct,       //< local shape function for Gauss Point (from fluid element)
    double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
    LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
    LINALG::Matrix<3, 1>& normal,     //< normal at gp
    LINALG::Matrix<3, 1>& vel_m       //< master velocity at gp
)
{
  double J = 0;
  double porosity = CalcPorosity(bele, rst_slave, J);

  double trperm = CalctrPermeability(bele, porosity, J);


  static double sliplength = trperm / (BJ_coeff_);

  static double dynvisc = (kappa_m * visc_m + (1.0 - kappa_m) * visc_s);

  static double stabnit = 0.0;
  static double stabadj = 0.0;

  XFEM::UTILS::GetNavierSlipStabilizationParameters(
      visc_stab_tang, dynvisc, sliplength, stabnit, stabadj);

  // Overall there are 9 coupling blocks to evaluate for fpi:
  // 1 - ps_ps --> ff,fps,psf,psps
  // 2 - ps_pf --> fpf,pspf
  // 3 - pf_ps --> pff, pfps
  // 4 - pf_pf --> pfpf
  switch (coupled_field_)
  {
    case MeshCouplingFPI::ps_ps:
    {
      configuration_map_[INPAR::XFEM::X_Adj_n_Col].second = 1.0 - porosity;
      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = stabadj;
        configuration_map_[INPAR::XFEM::X_Adj_t_Col].second = 1.0 - full_BJ_ * porosity;
        configuration_map_[INPAR::XFEM::FStr_Adj_t_Col].second = sliplength;
      }
      configuration_map_[INPAR::XFEM::F_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Col].second = 1.0 - porosity;
      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = stabnit;
        configuration_map_[INPAR::XFEM::X_Pen_t_Row].second = stabnit;
        configuration_map_[INPAR::XFEM::F_Con_t_Row].second = -stabnit;  //+sign for penalty!
        configuration_map_[INPAR::XFEM::F_Con_t_Col].second = sliplength / dynvisc;
        configuration_map_[INPAR::XFEM::X_Con_t_Row].second = -stabnit;  //+sign for penalty!
      }
      else
      {
        configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = 1. / sliplength;
        configuration_map_[INPAR::XFEM::X_Pen_t_Row].second = 1. / sliplength;

        // does nothing but should just be done in case we don't use the adjoint
        configuration_map_[INPAR::XFEM::F_Adj_t_Col].second = 1.0;
        configuration_map_[INPAR::XFEM::X_Adj_t_Col].second = 1.0 - full_BJ_ * porosity;
      }

      configuration_map_[INPAR::XFEM::X_Pen_t_Col].second = 1.0 - full_BJ_ * porosity;
      break;
    }
    case MeshCouplingFPI::ps_pf:
    {
      configuration_map_[INPAR::XFEM::X_Adj_n_Col].second = porosity;
      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = stabadj;
        configuration_map_[INPAR::XFEM::X_Adj_t_Col].second = full_BJ_ * porosity;
      }
      configuration_map_[INPAR::XFEM::F_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Col].second = porosity;
      if (!Sub_tang_)
      {
        configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = stabnit;
        configuration_map_[INPAR::XFEM::X_Pen_t_Row].second = stabnit;
      }
      else
      {
        configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = 1. / sliplength;
        configuration_map_[INPAR::XFEM::X_Pen_t_Row].second = 1. / sliplength;

        // does nothing but should just be done in case we don't use the adjoint
        configuration_map_[INPAR::XFEM::X_Adj_t_Col].second = full_BJ_ * porosity;
      }
      configuration_map_[INPAR::XFEM::X_Pen_t_Col].second = full_BJ_ * porosity;
      break;
    }
    case MeshCouplingFPI::pf_ps:
    {
      configuration_map_[INPAR::XFEM::X_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Col].second = 1.0 - porosity;

      // does nothing but should just be done in case we don't use the adjoint
      configuration_map_[INPAR::XFEM::F_Adj_n_Col].second =
          configuration_map_[INPAR::XFEM::F_Pen_n_Col].second;
      configuration_map_[INPAR::XFEM::X_Adj_n_Col].second =
          configuration_map_[INPAR::XFEM::X_Pen_n_Col].second;
      break;
    }
    case MeshCouplingFPI::pf_pf:
    {
      // Configuration of Penalty Terms
      configuration_map_[INPAR::XFEM::X_Pen_n_Row].second = full_stab;
      configuration_map_[INPAR::XFEM::X_Pen_n_Col].second = porosity;

      // does nothing but should just be done in case we don't use the adjoint
      configuration_map_[INPAR::XFEM::X_Adj_n_Col].second =
          configuration_map_[INPAR::XFEM::X_Pen_n_Col].second;
      break;
    }
  }

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::ZeroStateVectors_FPI()
{
  itrueresidual_->PutScalar(0.0);
  iforcecol_->PutScalar(0.0);
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFPI::ReadRestart(const int step)
{
  if (myrank_) IO::cout << "ReadRestart for boundary discretization " << IO::endl;

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");
  //  const int    step = boundaryreader.ReadInt("step");

  if (myrank_ == 0)
  {
    IO::cout << "time: " << time << IO::endl;
    IO::cout << "step: " << step << IO::endl;
  }

  boundaryreader.ReadVector(iveln_, "iveln_res");
  boundaryreader.ReadVector(idispn_, "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.ReadVector(ivelnp_, "ivelnp_res");
  boundaryreader.ReadVector(idispnp_, "idispnp_res");
  boundaryreader.ReadVector(idispnpi_, "idispnpi_res");

  if (not(cutter_dis_->DofRowMap())->SameAs(ivelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(iveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispn_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not(cutter_dis_->DofRowMap())->SameAs(idispnpi_->Map()))
    dserror("Global dof numbering in maps does not match");
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::GmshOutput(const std::string& filename_base, const int step,
    const int gmsh_step_diff, const bool gmsh_debug_out_screen)
{
  std::ostringstream filename_base_fsi;
  filename_base_fsi << filename_base << "_force";

  // compute the current boundary position
  std::map<int, LINALG::Matrix<3, 1>> currinterfacepositions;
  XFEM::UTILS::ExtractNodeVectors(cutter_dis_, currinterfacepositions, idispnp_);


  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(
      filename_base_fsi.str(), step, gmsh_step_diff, gmsh_debug_out_screen, myrank_);

  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "iforce \" {" << std::endl;
    // draw vector field 'force' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, itrueresidual_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "idispnp \" {" << std::endl;
    // draw vector field 'idispnp' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, idispnp_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "ivelnp \" {" << std::endl;
    // draw vector field 'ivelnp' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(
        cutter_dis_, ivelnp_, currinterfacepositions, gmshfilecontent, 3, 3);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::GmshOutputDiscretization(std::ostream& gmshfilecontent)
{
  // print surface discretization
  XFEM::MeshCoupling::GmshOutputDiscretization(gmshfilecontent);

  // compute the current solid and boundary position
  std::map<int, LINALG::Matrix<3, 1>> currsolidpositions;

  // write dis with zero solid displacements here!
  Teuchos::RCP<Epetra_Vector> solid_dispnp = LINALG::CreateVector(*cond_dis_->DofRowMap(), true);

  XFEM::UTILS::ExtractNodeVectors(cond_dis_, currsolidpositions, solid_dispnp);

  XFEM::UTILS::PrintDiscretizationToStream(cond_dis_, cond_dis_->Name(), true, false, true, false,
      false, false, gmshfilecontent, &currsolidpositions);
}

void XFEM::MeshCouplingFPI::Output(const int step, const double time, const bool write_restart_data)
{
  // output for interface
  cutter_output_->NewStep(step, time);

  cutter_output_->WriteVector("ivelnp", ivelnp_);
  cutter_output_->WriteVector("idispnp", idispnp_);
  cutter_output_->WriteVector("itrueresnp", itrueresidual_);

  cutter_output_->WriteElementData(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    cutter_output_->WriteVector("iveln_res", iveln_);
    cutter_output_->WriteVector("idispn_res", idispn_);
    cutter_output_->WriteVector("ivelnp_res", ivelnp_);
    cutter_output_->WriteVector("idispnp_res", idispnp_);
    cutter_output_->WriteVector("idispnpi_res", idispnpi_);
  }
}

void XFEM::MeshCouplingFPI::SetConditionSpecificParameters()
{
  std::vector<DRT::Condition*> conditions_XFPI;
  cutter_dis_->GetCondition(cond_name_, conditions_XFPI);

  // Create maps for easy extraction at gausspoint level
  for (std::vector<DRT::Condition*>::iterator i = conditions_XFPI.begin();
       i != conditions_XFPI.end(); ++i)
  {
    const bool full_BJ = (*(*i)->Get<std::string>("Variant") == "BJ");
    const bool Sub_tang = (*(*i)->Get<std::string>("Method") == "SUB");
    if (i != conditions_XFPI.begin())
    {
      if (fabs(BJ_coeff_ - (*i)->GetDouble("bj_coeff")) > 1e-16)
        dserror(
            "XFEM::MeshCouplingFPI::SetConditionSpecificParameters: You defined two FPI "
            "conditions, with different BJ_coeff!");

      if (full_BJ_ != full_BJ)
        dserror(
            "XFEM::MeshCouplingFPI::SetConditionSpecificParameters: You defined two FPI "
            "conditions, with different BJ Variant!");

      if (Sub_tang_ != Sub_tang)
        dserror(
            "XFEM::MeshCouplingFPI::SetConditionSpecificParameters: You defined two FPI "
            "conditions, with different BJ Method!");
    }

    BJ_coeff_ = (*i)->GetDouble("bj_coeff");
    full_BJ_ = full_BJ;
    Sub_tang_ = Sub_tang;
  }

  if (!full_BJ_)
  {
    if (coupled_field_ == XFEM::MeshCouplingFPI::ps_ps)
      std::cout << "==| XFEM::MeshCouplingFPI: Actual FPI Formulation is Beavers Joseph Saffmann"
                << std::flush;
  }
  else
  {
    if (coupled_field_ == XFEM::MeshCouplingFPI::ps_ps)
      std::cout << "==| XFEM::MeshCouplingFPI: Actual FPI Formulation is Beavers Joseph"
                << std::flush;
  }
  if (!Sub_tang_)
  {
    if (coupled_field_ == XFEM::MeshCouplingFPI::ps_ps)
      std::cout << " -- by tangential Nitsche formulation! |==" << std::endl;
  }
  else
  {
    if (coupled_field_ == XFEM::MeshCouplingFPI::ps_ps)
      std::cout << " -- by tangential Substitution formulation! |==" << std::endl;
  }
}

//----------------------------------------------------------------------
// LiftDrag                                                  chfoe 11/07
//----------------------------------------------------------------------
// calculate lift&drag forces
//
// Lift and drag forces are based upon the right hand side true-residual entities
// of the corresponding nodes. The contribution of the end node of a line is entirely
// added to a present L&D force.
/*----------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::LiftDrag(const int step, const double time) const
{
  // get forces on all procs
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<const Epetra_Vector> iforcecol =
      DRT::UTILS::GetColVersionOfRowVector(cutter_dis_, itrueresidual_);

  if (myrank_ == 0)
  {
    // compute force components
    const int nsd = 3;
    const Epetra_Map* dofcolmap = cutter_dis_->DofColMap();
    LINALG::Matrix<3, 1> c(true);
    for (int inode = 0; inode < cutter_dis_->NumMyColNodes(); ++inode)
    {
      const DRT::Node* node = cutter_dis_->lColNode(inode);
      const std::vector<int> dof = cutter_dis_->Dof(node);
      for (int isd = 0; isd < nsd; ++isd)
      {
        // [// minus to get correct sign of lift and drag (force acting on the body) ]
        c(isd) += (*iforcecol)[dofcolmap->LID(dof[isd])];
      }
    }

    // print to file
    std::ostringstream s;
    std::ostringstream header;

    header << std::left << std::setw(10) << "Time" << std::right << std::setw(16) << "F_x"
           << std::right << std::setw(16) << "F_y" << std::right << std::setw(16) << "F_z";
    s << std::left << std::setw(10) << std::scientific << time << std::right << std::setw(16)
      << std::scientific << c(0) << std::right << std::setw(16) << std::scientific << c(1)
      << std::right << std::setw(16) << std::scientific << c(2);

    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName() +
                              ".liftdrag." + cond_name_ + ".txt";
    if (step <= 1)
    {
      f.open(fname.c_str(), std::fstream::trunc);
      f << header.str() << std::endl;
    }
    else
    {
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
    }
    f << s.str() << "\n";
    f.close();

    std::cout << header.str() << std::endl << s.str() << std::endl;
  }
}

// ------------------------------------------------------------------------
// Caluculate the normalized trace of permeability matrix
//        for J,porosity pair on this FaceElement               ager 12/17
// ------------------------------------------------------------------------
double XFEM::MeshCouplingFPI::CalctrPermeability(DRT::Element* ele, double& porosity, double& J)
{
  // Calculate normalized trace of permeability matrix
  DRT::FaceElement* fele = dynamic_cast<DRT::FaceElement*>(ele);
  if (!fele) dserror("Cast to Faceele failed!");
  DRT::Element* coupl_ele = fele->ParentElement();
  if (coupl_ele == NULL) dserror("No coupl_ele!");
  Teuchos::RCP<MAT::FluidPoro> poromat;
  // access second material in structure element
  if (coupl_ele->NumMaterial() > 1)
    poromat = Teuchos::rcp_dynamic_cast<MAT::FluidPoro>(coupl_ele->Material(1));
  else
    dserror("no second material defined for element %i", ele->Id());

  static LINALG::Matrix<3, 3> reactiontensor(true);
  poromat->ComputeReactionTensor(reactiontensor, J, porosity);

  return sqrt((1. / reactiontensor(0, 0) + 1. / reactiontensor(1, 1) + 1. / reactiontensor(2, 2)) /
              (poromat->Viscosity() * 3));
}

// --------------------------------------------------------------------
// Caluculate the Porosity for this FaceElement Gausspoint   ager 12/16
// --------------------------------------------------------------------
double XFEM::MeshCouplingFPI::CalcPorosity(
    DRT::Element* ele, LINALG::Matrix<3, 1>& rst_slave, double& J)
{
  DRT::FaceElement* fele = dynamic_cast<DRT::FaceElement*>(ele);
  if (!fele) dserror("Cast to Faceele failed!");

  DRT::Element* coupl_ele = fele->ParentElement();
  if (coupl_ele == NULL) dserror("No coupl_ele!");

  double pres = 0.0;
  J = ComputeJacobianandPressure(ele, rst_slave, pres);

  Teuchos::RCP<MAT::StructPoro> poromat;
  // access second material in structure element
  if (coupl_ele->NumMaterial() > 1)
  {
    poromat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(coupl_ele->Material(0));
    if (poromat->MaterialType() != INPAR::MAT::m_structporo and
        poromat->MaterialType() != INPAR::MAT::m_structpororeaction and
        poromat->MaterialType() != INPAR::MAT::m_structpororeactionECM)
      dserror("invalid structure material for poroelasticity");
  }
  else
    dserror("no second material defined for element %i", ele->Id());

  Teuchos::ParameterList params;  // empty parameter list;
  double porosity;
  poromat->ComputePorosity(params, pres, J,
      1,  // not used
      porosity, false);
  return porosity;
}

// ---------------------------------------------------------------------------------------
// Compute Jacobian and extract PoroFluidPressure this FaceElement Gausspoint   ager 12/17
// ------------------------------------------------------------------------------------------
double XFEM::MeshCouplingFPI::ComputeJacobianandPressure(
    DRT::Element* ele, LINALG::Matrix<3, 1>& rst_slave, double& pres)
{
  DRT::FaceElement* fele = dynamic_cast<DRT::FaceElement*>(ele);
  if (!fele) dserror("Cast to Faceele failed!");

  DRT::Element* coupl_ele = fele->ParentElement();

  if (fele->Shape() == DRT::Element::quad4)
  {
    pres = 0.0;

    const unsigned int SLAVE_NUMDOF = 3;

    DRT::UTILS::CollectedGaussPoints intpoints =
        DRT::UTILS::CollectedGaussPoints(1);  // reserve just for 1 entry ...
    intpoints.Append(rst_slave(0, 0), rst_slave(1, 0), 0.0, 1.0);

    // get coordinates of gauss point w.r.t. local parent coordinate system
    LINALG::SerialDenseMatrix pqxg(1, SLAVE_NUMDOF);
    LINALG::Matrix<SLAVE_NUMDOF, SLAVE_NUMDOF> derivtrafo(true);

    DRT::UTILS::BoundaryGPToParentGP<SLAVE_NUMDOF>(
        pqxg, derivtrafo, intpoints, coupl_ele->Shape(), fele->Shape(), fele->FaceParentNumber());

    LINALG::Matrix<SLAVE_NUMDOF, 1> pxsi(true);

    // coordinates of the current integration point in parent coordinate system
    for (uint idim = 0; idim < SLAVE_NUMDOF; idim++)
    {
      pxsi(idim) = pqxg(0, idim);
    }
    if (coupl_ele->Shape() == DRT::Element::hex8)
    {
      const size_t PARENT_NEN =
          DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::Matrix<PARENT_NEN, 1> pfunc_loc(
          true);  // derivatives of parent element shape functions in parent element coordinate
                  // system
      LINALG::Matrix<SLAVE_NUMDOF, PARENT_NEN> pderiv_loc(
          true);  // derivatives of parent element shape functions in parent element coordinate
                  // system

      // evaluate derivatives of parent element shape functions at current integration point in
      // parent coordinate system
      DRT::UTILS::shape_function<DRT::Element::hex8>(pxsi, pfunc_loc);
      DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(pxsi, pderiv_loc);
      //
      // get Jacobian matrix and determinant w.r.t. spatial configuration
      //
      // |J| = det(xjm) * det(Jmat^-1) = det(xjm) * 1/det(Jmat)
      //
      //    _                     _
      //   |  x_1,1  x_2,1  x_3,1  |           d x_i
      //   |  x_1,2  x_2,2  x_3,2  | = xjm  = --------
      //   |_ x_1,3  x_2,3  x_3,3 _|           d s_j
      //    _                     _
      //   |  X_1,1  X_2,1  X_3,1  |           d X_i
      //   |  X_1,2  X_2,2  X_3,2  | = Jmat = --------
      //   |_ X_1,3  X_2,3  X_3,3 _|           d s_j
      //
      LINALG::Matrix<SLAVE_NUMDOF, SLAVE_NUMDOF> xjm;
      LINALG::Matrix<SLAVE_NUMDOF, SLAVE_NUMDOF> Jmat;

      LINALG::Matrix<SLAVE_NUMDOF, PARENT_NEN> xrefe(true);  // material coord. of parent element
      LINALG::Matrix<SLAVE_NUMDOF, PARENT_NEN> xcurr(true);  // current  coord. of parent element

      // update element geometry of parent element
      {
        DRT::Node** nodes = coupl_ele->Nodes();
        for (uint inode = 0; inode < PARENT_NEN; ++inode)
        {
          for (unsigned int idof = 0; idof < SLAVE_NUMDOF; ++idof)
          {
            int lid = fulldispnp_->Map().LID(GetCondDis()->Dof(0, coupl_ele->Nodes()[inode], idof));

            const double* x = nodes[inode]->X();
            xrefe(idof, inode) = x[idof];

            if (lid != -1)
              xcurr(idof, inode) = xrefe(idof, inode) + fulldispnp_->operator[](lid);
            else
              dserror("Local ID for dispnp not found (lid = -1)!");
          }
          int lidp = fullpres_->Map().LID(
              lm_struct_x_lm_pres_.operator[](GetCondDis()->Dof(0, coupl_ele->Nodes()[inode], 2)));

          if (lidp != -1)
            pres += fullpres_->operator[](lidp) * pfunc_loc(inode);
          else
            dserror("Local ID for pressure not found (lid = -1)!");
        }
      }

      xjm.MultiplyNT(pderiv_loc, xcurr);
      Jmat.MultiplyNT(pderiv_loc, xrefe);
      double det = xjm.Determinant();
      double detJ = Jmat.Determinant();
      double J = det / detJ;
      return J;
    }
    else
      dserror(
          "TDetDeformationGradient for type %s not yet implemented, just add your element type!",
          (DRT::DistypeToString(coupl_ele->Shape())).c_str());
    return -1.0;
  }
  else
    dserror("TDetDeformationGradient for type %s not yet implemented, just add your element type!",
        (DRT::DistypeToString(fele->Shape())).c_str());
  return -1.0;
}
