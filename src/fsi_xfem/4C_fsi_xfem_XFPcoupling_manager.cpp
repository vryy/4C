/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Manager for eXtended Fluid Poro Coupling

\level 3


*----------------------------------------------------------------------*/
#include "4C_fsi_xfem_XFPcoupling_manager.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_xfluid.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_poroelast_base.hpp"
#include "4C_xfem_condition_manager.hpp"

FOUR_C_NAMESPACE_OPEN

XFEM::XfpCouplingManager::XfpCouplingManager(Teuchos::RCP<XFEM::ConditionManager> condmanager,
    Teuchos::RCP<PoroElast::PoroBase> poro, Teuchos::RCP<FLD::XFluid> xfluid, std::vector<int> idx)
    : CouplingCommManager(poro->structure_field()->discretization(),
          poro->fluid_field()->discretization(), "XFEMSurfFPIMono", 0, 3),
      poro_(poro),
      xfluid_(xfluid),
      cond_name_ps_ps_("XFEMSurfFPIMono_ps_ps"),
      cond_name_ps_pf_("XFEMSurfFPIMono_ps_pf"),
      cond_name_pf_ps_("XFEMSurfFPIMono_pf_ps"),
      cond_name_pf_pf_("XFEMSurfFPIMono_pf_pf"),
      idx_(idx)
{
  // Coupling_Comm_Manager create all Coupling Objects now with Structure has idx = 0, Fluid has idx
  // = 1!

  mcfpi_ps_ps_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(
      condmanager->GetMeshCoupling(cond_name_ps_ps_));
  if (mcfpi_ps_ps_ == Teuchos::null)
    FOUR_C_THROW(" Failed to get MeshCouplingFPI for Porostructure!");
  mcfpi_ps_ps_->initialize_struc_pres_map(poro_->fluid_structure_coupling().SlaveDofMap(),
      poro_->fluid_structure_coupling().PermMasterDofMap());

  mcfpi_ps_pf_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(
      condmanager->GetMeshCoupling(cond_name_ps_pf_));
  if (mcfpi_ps_pf_ == Teuchos::null) FOUR_C_THROW(" Failed to get MeshCouplingFPI for Porofluid!");
  mcfpi_ps_pf_->initialize_struc_pres_map(poro_->fluid_structure_coupling().SlaveDofMap(),
      poro_->fluid_structure_coupling().PermMasterDofMap());

  mcfpi_pf_ps_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(
      condmanager->GetMeshCoupling(cond_name_pf_ps_));
  if (mcfpi_pf_ps_ == Teuchos::null) FOUR_C_THROW(" Failed to get MeshCouplingFPI for Porofluid!");
  mcfpi_pf_ps_->initialize_struc_pres_map(poro_->fluid_structure_coupling().SlaveDofMap(),
      poro_->fluid_structure_coupling().PermMasterDofMap());

  mcfpi_pf_pf_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(
      condmanager->GetMeshCoupling(cond_name_pf_pf_));
  if (mcfpi_pf_pf_ == Teuchos::null) FOUR_C_THROW(" Failed to get MeshCouplingFPI for Porofluid!");
  mcfpi_pf_pf_->initialize_struc_pres_map(poro_->fluid_structure_coupling().SlaveDofMap(),
      poro_->fluid_structure_coupling().PermMasterDofMap());

  // safety check
  if (!mcfpi_ps_ps_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    FOUR_C_THROW("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (psps)!");
  if (!mcfpi_ps_pf_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    FOUR_C_THROW("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (pspf)!");
  if (!mcfpi_pf_ps_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    FOUR_C_THROW("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (pfps)!");
  if (!mcfpi_pf_pf_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    FOUR_C_THROW("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (pfpf)!");

  // storage of the resulting Robin-type structural forces from the old timestep
  // Recovering of Lagrange multiplier happens on fluid field
  lambda_ps_ = Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->Map(1), true));
  lambda_pf_ = Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->Map(1), true));
}

void XFEM::XfpCouplingManager::InitCouplingStates()
{
  mcfpi_ps_ps_->reconnect_parent_pointers();
  mcfpi_pf_ps_->reconnect_parent_pointers();
  mcfpi_ps_pf_->reconnect_parent_pointers();
  mcfpi_pf_pf_->reconnect_parent_pointers();
}

void XFEM::XfpCouplingManager::SetCouplingStates()
{
  // 1 Set Displacement on both mesh couplings ... we get them from the structure field!
  InsertVector(0, poro_->structure_field()->Dispnp(), 0, mcfpi_ps_ps_->IDispnp(),
      CouplingCommManager::full_to_partial);
  InsertVector(0, poro_->structure_field()->Dispnp(), 0, mcfpi_ps_pf_->IDispnp(),
      CouplingCommManager::full_to_partial);
  InsertVector(0, poro_->structure_field()->Dispnp(), 0, mcfpi_pf_ps_->IDispnp(),
      CouplingCommManager::full_to_partial);
  InsertVector(0, poro_->structure_field()->Dispnp(), 0, mcfpi_pf_pf_->IDispnp(),
      CouplingCommManager::full_to_partial);

  // As interfaces embedded into the background mesh are fully ghosted, we don't know which
  Teuchos::RCP<Epetra_Map> sfulldofmap =
      Core::LinAlg::AllreduceEMap(*poro_->structure_field()->discretization()->dof_row_map());
  Teuchos::RCP<Epetra_Vector> dispnp_col = Teuchos::rcp(new Epetra_Vector(*sfulldofmap, true));
  Core::LinAlg::Export(*poro_->structure_field()->Dispnp(), *dispnp_col);
  Teuchos::RCP<Epetra_Map> ffulldofmap =
      Core::LinAlg::AllreduceEMap(*poro_->fluid_field()->discretization()->dof_row_map());
  Teuchos::RCP<Epetra_Vector> velnp_col = Teuchos::rcp(new Epetra_Vector(*ffulldofmap, true));
  Core::LinAlg::Export(*poro_->fluid_field()->Velnp(), *velnp_col);

  mcfpi_ps_ps_->SetFullState(dispnp_col, velnp_col);
  mcfpi_ps_pf_->SetFullState(dispnp_col, velnp_col);
  mcfpi_pf_ps_->SetFullState(dispnp_col, velnp_col);
  mcfpi_pf_pf_->SetFullState(dispnp_col, velnp_col);


  // 2 Set Structural Velocity onto ps mesh coupling
  InsertVector(0, poro_->structure_field()->Velnp(), 0, mcfpi_ps_ps_->IVelnp(),
      CouplingCommManager::full_to_partial);
  InsertVector(0, poro_->structure_field()->Velnp(), 0, mcfpi_pf_ps_->IVelnp(),
      CouplingCommManager::full_to_partial);
  //  poro_->structure_field()->Velnp()->Print(std::cout);

  //  InsertVector(1,poro_->fluid_field()->GridVel(),0,mcfpi_ps_ps_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
  //  InsertVector(1,poro_->fluid_field()->GridVel(),0,mcfpi_pf_ps_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
  // 3 Set Fluid Velocity onto pf mesh coupling
  InsertVector(1, poro_->fluid_field()->Velnp(), 0, mcfpi_ps_pf_->IVelnp(),
      CouplingCommManager::full_to_partial);
  InsertVector(1, poro_->fluid_field()->Velnp(), 0, mcfpi_pf_pf_->IVelnp(),
      CouplingCommManager::full_to_partial);
}

void XFEM::XfpCouplingManager::AddCouplingMatrix(
    Core::LinAlg::BlockSparseMatrixBase& systemmatrix, double scaling)
{
  const double scaling_disp_vel =
      1 / ((1 - poro_->structure_field()->TimIntParam()) * poro_->structure_field()->Dt());
  const double dt = poro_->fluid_field()->Dt();
  if (idx_.size() == 2)  // assum that the poro field is not split and we just have a blockmatrix
                         // P/F
  {
    Core::LinAlg::SparseMatrix& C_ss_block = (systemmatrix)(idx_[0], idx_[0]);
    Core::LinAlg::SparseMatrix& C_fs_block = (systemmatrix)(idx_[1], idx_[0]);
    Core::LinAlg::SparseMatrix& C_sf_block = (systemmatrix)(idx_[0], idx_[1]);

    // 1// Add Blocks f-ps(2), ps-f(3), ps-ps(4)
    C_ss_block.Add(*xfluid_->C_ss_Matrix(cond_name_ps_ps_), false, scaling * scaling_disp_vel, 1.0);
    C_sf_block.Add(*xfluid_->C_sx_Matrix(cond_name_ps_ps_), false, scaling, 1.0);
    C_fs_block.Add(*xfluid_->C_xs_Matrix(cond_name_ps_ps_), false, scaling * scaling_disp_vel, 1.0);

    // 2// Add Blocks f-pf(5), ps-pf(6)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> C_ps_pf = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        xfluid_->C_ss_Matrix(cond_name_ps_pf_)->RowMap(), 81, false));
    InsertMatrix(-1, 0, *xfluid_->C_ss_Matrix(cond_name_ps_pf_), 1, *C_ps_pf,
        CouplingCommManager::col, 1, true, false);
    C_ps_pf->Complete(*GetMapExtractor(1)->Map(1), *GetMapExtractor(0)->Map(1));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> C_f_pf = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        xfluid_->C_xs_Matrix(cond_name_ps_pf_)->RowMap(), 81, false));
    InsertMatrix(-1, 0, *xfluid_->C_xs_Matrix(cond_name_ps_pf_), 1, *C_f_pf,
        CouplingCommManager::col, 1, true, false);
    C_f_pf->Complete(*GetMapExtractor(1)->Map(1), C_fs_block.RangeMap());
    C_fs_block.Add(*C_f_pf, false, scaling, 1.0);
    C_ss_block.Add(*C_ps_pf, false, scaling, 1.0);

    // 3// Add Blocks pf-f(7), pf-ps(8)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> C_pf_ps =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*GetMapExtractor(1)->Map(1), 81, false));
    Teuchos::RCP<Core::LinAlg::SparseMatrix> C_pf_f =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*GetMapExtractor(1)->Map(1), 81, false));
    InsertMatrix(-1, 0, *xfluid_->C_ss_Matrix(cond_name_pf_ps_), 1, *C_pf_ps,
        CouplingCommManager::row, 1, true, false);
    C_pf_ps->Complete(*GetMapExtractor(0)->Map(1), *GetMapExtractor(1)->Map(1));
    InsertMatrix(-1, 0, *xfluid_->C_sx_Matrix(cond_name_pf_ps_), 1, *C_pf_f,
        CouplingCommManager::row, 1, true, false);
    C_pf_f->Complete(*xfluid_->dof_row_map(), *GetMapExtractor(1)->Map(1));
    C_ss_block.Add(*C_pf_ps, false, scaling * scaling_disp_vel * dt, 1.0);
    C_sf_block.Add(*C_pf_f, false, scaling * dt, 1.0);

    // 4// Add Block pf-pf(9)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> C_pf_pf =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*GetMapExtractor(1)->Map(1), 81, false));
    InsertMatrix(-1, 0, *xfluid_->C_ss_Matrix(cond_name_pf_pf_), 1, *C_pf_pf,
        CouplingCommManager::row_and_col);
    C_pf_pf->Complete(*GetMapExtractor(1)->Map(1), *GetMapExtractor(1)->Map(1));
    C_ss_block.Add(*C_pf_pf, false, scaling * dt, 1.0);
  }
  else if (idx_.size() == 3)
  {
    Core::LinAlg::SparseMatrix& C_ss_block = (systemmatrix)(idx_[0], idx_[0]);
    Core::LinAlg::SparseMatrix& C_fs_block = (systemmatrix)(idx_[1], idx_[0]);
    Core::LinAlg::SparseMatrix& C_sf_block = (systemmatrix)(idx_[0], idx_[1]);

    Core::LinAlg::SparseMatrix& C_pfpfblock = (systemmatrix)(idx_[2], idx_[2]);
    Core::LinAlg::SparseMatrix& C_fpf_block = (systemmatrix)(idx_[1], idx_[2]);
    Core::LinAlg::SparseMatrix& C_pff_block = (systemmatrix)(idx_[2], idx_[1]);

    Core::LinAlg::SparseMatrix& C_pfs_block = (systemmatrix)(idx_[2], idx_[0]);
    Core::LinAlg::SparseMatrix& C_spf_block = (systemmatrix)(idx_[0], idx_[2]);

    // 1// Add Blocks f-ps(2), ps-f(3), ps-ps(4)
    C_ss_block.Add(*xfluid_->C_ss_Matrix(cond_name_ps_ps_), false, scaling * scaling_disp_vel, 1.0);
    C_sf_block.Add(*xfluid_->C_sx_Matrix(cond_name_ps_ps_), false, scaling, 1.0);
    C_fs_block.Add(*xfluid_->C_xs_Matrix(cond_name_ps_ps_), false, scaling * scaling_disp_vel, 1.0);

    // 2// Add Blocks f-pf(5), ps-pf(6)
    InsertMatrix(-1, 0, *xfluid_->C_ss_Matrix(cond_name_ps_pf_), 1, C_spf_block,
        CouplingCommManager::col, scaling, true, true);
    InsertMatrix(-1, 0, *xfluid_->C_xs_Matrix(cond_name_ps_pf_), 1, C_fpf_block,
        CouplingCommManager::col, scaling, true, true);

    // 3// Add Blocks pf-f(7), pf-ps(8)
    InsertMatrix(-1, 0, *xfluid_->C_ss_Matrix(cond_name_pf_ps_), 1, C_pfs_block,
        CouplingCommManager::row, scaling * scaling_disp_vel * dt, true, true);
    InsertMatrix(-1, 0, *xfluid_->C_sx_Matrix(cond_name_pf_ps_), 1, C_pff_block,
        CouplingCommManager::row, scaling * dt, true, true);

    // 4// Add Block pf-pf(9)
    InsertMatrix(-1, 0, *xfluid_->C_ss_Matrix(cond_name_pf_pf_), 1, C_pfpfblock,
        CouplingCommManager::row_and_col, scaling * dt, true, true);
  }
  else
    FOUR_C_THROW(
        "XFPCoupling_Manager::AddCouplingMatrix: Not implemented for number of blocks = %d",
        idx_.size());
}

///*-----------------------------------------------------------------------------------------*
//| Add the coupling rhs                                                        ager 06/2016 |
//*-----------------------------------------------------------------------------------------*/
void XFEM::XfpCouplingManager::AddCouplingRHS(
    Teuchos::RCP<Epetra_Vector> rhs, const Core::LinAlg::MultiMapExtractor& me, double scaling)
{
  const double dt = poro_->fluid_field()->Dt();
  if (idx_.size() == 2)  // assum that the poro field is not split and we just have a blockmatrix
                         // P/F
  {
    Teuchos::RCP<const Epetra_Vector> rhs_C_ps_ps = xfluid_->RHS_s_Vec(cond_name_ps_ps_);
    Teuchos::RCP<const Epetra_Vector> rhs_C_ps_pf = xfluid_->RHS_s_Vec(cond_name_ps_pf_);
    Teuchos::RCP<const Epetra_Vector> rhs_C_pf_ps = xfluid_->RHS_s_Vec(cond_name_pf_ps_);
    Teuchos::RCP<const Epetra_Vector> rhs_C_pf_pf = xfluid_->RHS_s_Vec(cond_name_pf_pf_);

    Teuchos::RCP<Epetra_Vector> prhs = Teuchos::rcp(new Epetra_Vector(*me.Map(idx_[0]), true));

    InsertVector(0, rhs_C_ps_ps, 0, prhs, CouplingCommManager::partial_to_global, true, scaling);
    InsertVector(0, rhs_C_ps_pf, 0, prhs, CouplingCommManager::partial_to_global, true, scaling);

    InsertVector(
        0, rhs_C_pf_ps, 1, prhs, CouplingCommManager::partial_to_global, true, scaling * dt);
    InsertVector(
        0, rhs_C_pf_pf, 1, prhs, CouplingCommManager::partial_to_global, true, scaling * dt);

    // Add lambda contribution
    if (lambda_ps_ != Teuchos::null && lambda_pf_ != Teuchos::null)
    {
      /*----------------------------------------------------------------------*/
      // get time integration parameters of structure and fluid time integrators
      // to enable consistent time integration among the fields
      /*----------------------------------------------------------------------*/

      /*----------------------------------------------------------------------*/
      // this is the interpolation weight for quantities from last time step
      // alpha_f for genalpha and (1-theta) for OST (weighting of the old time step n for
      // displacements) TimeIntegration for poro needs to be consistent!
      const double stiparam =
          poro_->structure_field()->TimIntParam();  // (1-theta) for OST and alpha_f for Genalpha

      // scale factor for the structure system matrix w.r.t the new time step
      const double scaling_S = 1.0 / (1.0 - stiparam);  // 1/(1-alpha_F) = 1/weight^S_np

      InsertVector(0, lambda_ps_, 0, prhs, CouplingCommManager::partial_to_global, true,
          stiparam * scaling_S);
      InsertVector(0, lambda_pf_, 1, prhs, CouplingCommManager::partial_to_global, true,
          stiparam * scaling_S);
    }

    me.AddVector(prhs, idx_[0], rhs);
  }
  else if (idx_.size() == 3)
  {
    Teuchos::RCP<const Epetra_Vector> rhs_C_ps_ps = xfluid_->RHS_s_Vec(cond_name_ps_ps_);
    Teuchos::RCP<const Epetra_Vector> rhs_C_ps_pf = xfluid_->RHS_s_Vec(cond_name_ps_pf_);
    Teuchos::RCP<const Epetra_Vector> rhs_C_pf_ps = xfluid_->RHS_s_Vec(cond_name_pf_ps_);
    Teuchos::RCP<const Epetra_Vector> rhs_C_pf_pf = xfluid_->RHS_s_Vec(cond_name_pf_pf_);

    Teuchos::RCP<Epetra_Vector> srhs = Teuchos::rcp(new Epetra_Vector(*me.Map(idx_[0]), true));
    Teuchos::RCP<Epetra_Vector> pfrhs = Teuchos::rcp(new Epetra_Vector(*me.Map(idx_[2]), true));

    InsertVector(0, rhs_C_ps_ps, 0, srhs, CouplingCommManager::partial_to_full, true, scaling);
    InsertVector(0, rhs_C_ps_pf, 0, srhs, CouplingCommManager::partial_to_full, true, scaling);

    InsertVector(
        0, rhs_C_pf_ps, 1, pfrhs, CouplingCommManager::partial_to_full, true, scaling * dt);
    InsertVector(
        0, rhs_C_pf_pf, 1, pfrhs, CouplingCommManager::partial_to_full, true, scaling * dt);

    // Add lambda contribution
    if (lambda_ps_ != Teuchos::null && lambda_pf_ != Teuchos::null)
    {
      /*----------------------------------------------------------------------*/
      // get time integration parameters of structure and fluid time integrators
      // to enable consistent time integration among the fields
      /*----------------------------------------------------------------------*/

      /*----------------------------------------------------------------------*/
      // this is the interpolation weight for quantities from last time step
      // alpha_f for genalpha and (1-theta) for OST (weighting of the old time step n for
      // displacements) TimeIntegration for poro needs to be consistent!
      const double stiparam =
          poro_->structure_field()->TimIntParam();  // (1-theta) for OST and alpha_f for Genalpha

      // scale factor for the structure system matrix w.r.t the new time step
      const double scaling_S = 1.0 / (1.0 - stiparam);  // 1/(1-alpha_F) = 1/weight^S_np

      InsertVector(
          0, lambda_ps_, 0, srhs, CouplingCommManager::partial_to_full, true, stiparam * scaling_S);
      InsertVector(0, lambda_pf_, 1, pfrhs, CouplingCommManager::partial_to_full, true,
          stiparam * scaling_S);
    }

    me.AddVector(srhs, idx_[0], rhs);
    me.AddVector(pfrhs, idx_[2], rhs);
  }
  else
    FOUR_C_THROW("XFPCoupling_Manager::AddCouplingRHS: Not implemented for number of blocks = %d",
        idx_.size());
}

/*----------------------------------------------------------------------*/
/* Store the Coupling RHS of the Old Timestep in lambda     ager 06/2016 |
 *----------------------------------------------------------------------*/
void XFEM::XfpCouplingManager::Update(double scaling)
{
  /*----------------------------------------------------------------------*/
  // we directly store the fluid-unscaled rhs_C_s residual contribution from the fluid solver which
  // corresponds to the actual acting forces

  // scaling for the structural residual is done when it is added to the global residual vector
  // get the coupling rhs from the xfluid, this vector is based on the boundary dis which is part of
  // the structure dis
  lambda_ps_->Update(scaling, *xfluid_->RHS_s_Vec(cond_name_ps_ps_), 0.0);
  lambda_ps_->Update(scaling, *xfluid_->RHS_s_Vec(cond_name_ps_pf_), 1.0);

  const double dt = poro_->fluid_field()->Dt();
  lambda_pf_->Update(scaling * dt, *xfluid_->RHS_s_Vec(cond_name_pf_ps_), 0.0);
  lambda_pf_->Update(scaling * dt, *xfluid_->RHS_s_Vec(cond_name_pf_pf_), 1.0);
  return;
}

/*----------------------------------------------------------------------*/
/* Write Output                                             ager 06/2016 |
 *-----------------------------------------------------------------------*/
void XFEM::XfpCouplingManager::Output(Core::IO::DiscretizationWriter& writer)
{
  //--------------------------------
  // output for Lagrange multiplier field (ie forces onto the structure, Robin-type forces
  // consisting of fluid forces and the Nitsche penalty term contribution)
  //--------------------------------
  Teuchos::RCP<Epetra_Vector> lambdafull =
      Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(), true));
  InsertVector(0, lambda_ps_, 0, lambdafull, CouplingCommManager::partial_to_full);
  writer.write_vector("fpilambda_ps", lambdafull);

  lambdafull = Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(), true));
  InsertVector(0, lambda_pf_, 0, lambdafull, CouplingCommManager::partial_to_full);
  writer.write_vector("fpilambda_pf", lambdafull);
  return;
}
/*----------------------------------------------------------------------*/
/* Read Restart on the interface                            ager 06/2016 |
 *-----------------------------------------------------------------------*/
void XFEM::XfpCouplingManager::read_restart(Core::IO::DiscretizationReader& reader)
{
  Teuchos::RCP<Epetra_Vector> lambdafull =
      Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(), true));
  reader.read_vector(lambdafull, "fpilambda_ps");
  InsertVector(0, lambdafull, 0, lambda_ps_, CouplingCommManager::full_to_partial);

  lambdafull = Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(), true));
  reader.read_vector(lambdafull, "fpilambda_pf");
  InsertVector(0, lambdafull, 0, lambda_pf_, CouplingCommManager::full_to_partial);
  return;
}

/*-----------------------------------------------------------------------------------------*
| Get Timeface on the interface (for OST this is 1/(theta dt))                ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
double XFEM::XfpCouplingManager::get_interface_timefac()
{
  FOUR_C_THROW("Check if you really want this!");
  /*
   * Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double dt = xfluid_->Dt();
  if (interface_second_order_)
    return 2. / dt;
  else
    return 1. / dt;
}

FOUR_C_NAMESPACE_CLOSE
