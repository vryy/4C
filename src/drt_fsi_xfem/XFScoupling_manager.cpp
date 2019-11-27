/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Manager for eXtended Fluid Structural Coupling

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

*----------------------------------------------------------------------*/
#include "XFScoupling_manager.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"

#include "../drt_fluid_xfluid/xfluidfluid.H"  //Todo: remove me finally

/*-----------------------------------------------------------------------------------------*
| Constructor                                                                 ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
XFEM::XFSCoupling_Manager::XFSCoupling_Manager(Teuchos::RCP<ConditionManager> condmanager,
    Teuchos::RCP<ADAPTER::Structure> structure, Teuchos::RCP<FLD::XFluid> xfluid,
    std::vector<int> idx)
    : Coupling_Comm_Manager(structure->Discretization(), "XFEMSurfFSIMono", 0, 3),
      struct_(structure),
      xfluid_(xfluid),
      cond_name_("XFEMSurfFSIMono"),
      idx_(idx),
      interface_second_order_(false)
{
  if (idx_.size() != 2) dserror("XFSCoupling_Manager required two block ( 2 != %d)", idx_.size());

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  interface_second_order_ = DRT::INPUT::IntegralValue<int>(fsidyn, "SECONDORDER");

  // Coupling_Comm_Manager create all Coupling Objects now with Structure has idx = 0, Fluid has idx
  // = 1!
  mcfsi_ =
      Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(condmanager->GetMeshCoupling(cond_name_));
  if (mcfsi_ == Teuchos::null) dserror(" Failed to get MeshCouplingFSI for Structure!");

  mcfsi_->SetTimeFac(1. / GetInterfaceTimefac());

  // safety check
  if (!mcfsi_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    dserror("XFSCoupling_Manager: Maps of Condition and Mesh Coupling do not fit!");

  // storage of the resulting Robin-type structural forces from the old timestep
  // Recovering of Lagrange multiplier happens on fluid field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*mcfsi_->GetCouplingDis()->DofRowMap(), true));
}

/*-----------------------------------------------------------------------------------------*
| Set required displacement & velocity states in the coupling object          ager 04/2017 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::InitCouplingStates()
{
  // 1 Set Displacement on both mesh couplings ... we get them from the structure field!
  InsertVector(0, struct_->Dispn(), 0, mcfsi_->IDispn(), Coupling_Comm_Manager::full_to_partial);
  InsertVector(0, struct_->Dispn(), 0, mcfsi_->IDispnp(), Coupling_Comm_Manager::full_to_partial);

  // 2 Set Displacement on both mesh couplings ... we get them from the structure field!
  InsertVector(0, struct_->Veln(), 0, mcfsi_->IVeln(), Coupling_Comm_Manager::full_to_partial);
  InsertVector(0, struct_->Veln(), 0, mcfsi_->IVelnp(), Coupling_Comm_Manager::full_to_partial);
}

/*-----------------------------------------------------------------------------------------*
| Set required displacement & velocity states in the coupling object          ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::SetCouplingStates()
{
  // 1 update last increment, before we set new idispnp
  mcfsi_->UpdateDisplacementIterationVectors();

  // 2 Set Displacement on both mesh couplings ... we get them from the structure field!
  InsertVector(0, struct_->Dispnp(), 0, mcfsi_->IDispnp(), Coupling_Comm_Manager::full_to_partial);

  // get interface velocity at t(n)
  Teuchos::RCP<Epetra_Vector> velnp =
      Teuchos::rcp(new Epetra_Vector(mcfsi_->IVelnp()->Map(), true));
  velnp->Update(1.0, *mcfsi_->IDispnp(), -1.0, *mcfsi_->IDispn(), 0.0);

  // inverse of FSI (1st order, 2nd order) scaling
  const double scaling_FSI = GetInterfaceTimefac();  // 1/(theta_FSI * dt) =  1/weight^FSI_np
  const double dt = xfluid_->Dt();

  // v^{n+1} = -(1-theta)/theta * v^{n} - 1/(theta*dt)*(d^{n+1}-d^{n0})
  velnp->Update(-(dt - 1 / scaling_FSI) * scaling_FSI, *mcfsi_->IVeln(), scaling_FSI);

  // 3 Set Structural Velocity onto ps mesh coupling
  InsertVector(0, velnp, 0, mcfsi_->IVelnp(), Coupling_Comm_Manager::partial_to_partial);

  // 4 Set Structural Velocity onto the structural discretization
  if (mcfsi_->GetAveragingStrategy() != INPAR::XFEM::Xfluid_Sided)
  {
    // Set Dispnp (used to calc local coord of gausspoints)
    struct_->Discretization()->SetState("dispnp", struct_->Dispnp());
    // Set Velnp (used for interface integration)
    Teuchos::RCP<Epetra_Vector> fullvelnp =
        Teuchos::rcp(new Epetra_Vector(struct_->Velnp()->Map(), true));
    fullvelnp->Update(1.0, *struct_->Dispnp(), -1.0, *struct_->Dispn(), 0.0);
    fullvelnp->Update(-(dt - 1 / scaling_FSI) * scaling_FSI, *struct_->Veln(), scaling_FSI);
    struct_->Discretization()->SetState("velaf", fullvelnp);
  }
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling matrixes to the global systemmatrix                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::AddCouplingMatrix(
    LINALG::BlockSparseMatrixBase& systemmatrix, double scaling)
{
  /*----------------------------------------------------------------------*/
  // Coupling blocks C_sf, C_fs and C_ss
  /*----------------------------------------------------------------------*/
  LINALG::SparseMatrix& C_ss_block = (systemmatrix)(idx_[0], idx_[0]);
  /*----------------------------------------------------------------------*/
  // scaling factor for displacement <-> velocity conversion (FSI)
  // inverse of FSI (1st order, 2nd order) scaling
  const double scaling_FSI = GetInterfaceTimefac();  // 1/(theta_FSI * dt) =  1/weight^FSI_np

  // * all the coupling matrices are scaled with the weighting of the fluid w.r.t new time step np
  //    -> Unscale the blocks with (1/(theta_f*dt) = 1/weight(t^f_np))
  // * additionally the C_*s blocks (C_ss and C_fs) have to include the conversion from structural
  // displacements to structural velocities
  //    -> Scale these blocks with (1/(theta_FSI*dt) = 1/weight(t^FSI_np))
  //
  // REMARK that Scale() scales the original coupling matrix in xfluid

  // C_ss_block scaled with 1/(theta_f*dt) * 1/(theta_FSI*dt) = 1/weight(t^f_np) *
  // 1/weight(t^FSI_np) add the coupling block C_ss on the already existing diagonal block
  C_ss_block.Add(*xfluid_->C_ss_Matrix(cond_name_), false, scaling * scaling_FSI, 1.0);


  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  // Todo: Need to eighter split fluid matrixes in the fsi algo or change the maps of the coupling
  // matrixes(merged)
  bool is_xfluidfluid =
      Teuchos::rcp_dynamic_cast<FLD::XFluidFluid>(xfluid_, false) != Teuchos::null;

  if (probtype == prb_fsi_xfem && !is_xfluidfluid)  // use assign for off diagonal blocks
  {
    // scale the off diagonal coupling blocks
    xfluid_->C_sx_Matrix(cond_name_)
        ->Scale(scaling);  //<   1/(theta_f*dt)                    = 1/weight(t^f_np)
    xfluid_->C_xs_Matrix(cond_name_)
        ->Scale(scaling * scaling_FSI);  //<   1/(theta_f*dt) * 1/(theta_FSI*dt) = 1/weight(t^f_np)
                                         //* 1/weight(t^FSI_np)

    systemmatrix.Assign(idx_[0], idx_[1], LINALG::View, *xfluid_->C_sx_Matrix(cond_name_));
    systemmatrix.Assign(idx_[1], idx_[0], LINALG::View, *xfluid_->C_xs_Matrix(cond_name_));
  }
  else if (probtype == prb_fpsi_xfem || is_xfluidfluid)
  {
    LINALG::SparseMatrix& C_fs_block = (systemmatrix)(idx_[1], idx_[0]);
    LINALG::SparseMatrix& C_sf_block = (systemmatrix)(idx_[0], idx_[1]);

    C_sf_block.Add(*xfluid_->C_sx_Matrix(cond_name_), false, scaling, 1.0);
    C_fs_block.Add(*xfluid_->C_xs_Matrix(cond_name_), false, scaling * scaling_FSI, 1.0);
  }
  else
  {
    dserror("XFSCoupling_Manager: Want to use me for other problemtype --> check and add me!");
  }
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling rhs                                                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::AddCouplingRHS(
    Teuchos::RCP<Epetra_Vector> rhs, const LINALG::MultiMapExtractor& me, double scaling)
{
  Teuchos::RCP<Epetra_Vector> coup_rhs_sum = Teuchos::rcp(new Epetra_Vector(*xfluid_->RHS_s_Vec(
      cond_name_)));  // REMARK: Copy this vector to store the correct lambda_ in update!
  /// Lagrange multiplier \lambda_\Gamma^n at the interface (ie forces onto the structure,
  /// Robin-type forces consisting of fluid forces and the Nitsche penalty term contribution)
  if (lambda_ != Teuchos::null)
  {
    /*----------------------------------------------------------------------*/
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    /*----------------------------------------------------------------------*/

    /*----------------------------------------------------------------------*/
    // this is the interpolation weight for quantities from last time step
    // alpha_f for genalpha and (1-theta) for OST (weighting of the old time step n for
    // displacements)
    const double stiparam = struct_->TimIntParam();  // (1-theta) for OST and alpha_f for Genalpha

    // scale factor for the structure system matrix w.r.t the new time step
    const double scaling_S = 1.0 / (1.0 - stiparam);  // 1/(1-alpha_F) = 1/weight^S_np
    // add Lagrange multiplier (structural forces from t^n)
    int err = coup_rhs_sum->Update(stiparam * scaling_S, *lambda_, scaling);
    if (err) dserror("Update of Nit_Struct_FSI RHS failed with errcode = %d!", err);
  }
  else
  {
    coup_rhs_sum->Scale(scaling);
  }

  Teuchos::RCP<Epetra_Vector> coup_rhs = Teuchos::rcp(new Epetra_Vector(*me.Map(idx_[0]), true));
  LINALG::Export(*coup_rhs_sum, *coup_rhs);  // use this command as long as poro ist not split into
                                             // two bocks in the monolithic algorithm!
  // InsertVector(0,coup_rhs_sum,0,coup_rhs,Coupling_Comm_Manager::partial_to_full);
  me.AddVector(coup_rhs, idx_[0], rhs);
}

/*----------------------------------------------------------------------*/
/* Store the Coupling RHS of the Old Timestep in lambda     ager 06/2016 |
 *----------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::Update(double scaling)
{
  /*----------------------------------------------------------------------*/
  // we directly store the fluid-unscaled rhs_C_s residual contribution from the fluid solver which
  // corresponds to the actual acting forces

  // scaling for the structural residual is done when it is added to the global residual vector
  // get the coupling rhs from the xfluid, this vector is based on the boundary dis which is part of
  // the structure dis
  lambda_->Update(scaling, *xfluid_->RHS_s_Vec(cond_name_), 0.0);
  return;
}

/*----------------------------------------------------------------------*/
/* Write Output                                             ager 06/2016 |
 *-----------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::Output(IO::DiscretizationWriter& writer)
{
  //--------------------------------
  // output for Lagrange multiplier field (ie forces onto the structure, Robin-type forces
  // consisting of fluid forces and the Nitsche penalty term contribution)
  //--------------------------------
  Teuchos::RCP<Epetra_Vector> lambdafull =
      Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(), true));
  InsertVector(0, lambda_, 0, lambdafull, Coupling_Comm_Manager::partial_to_full);
  writer.WriteVector("fsilambda", lambdafull);
  return;
}
/*----------------------------------------------------------------------*/
/* Read Restart on the interface                            ager 06/2016 |
 *-----------------------------------------------------------------------*/
void XFEM::XFSCoupling_Manager::ReadRestart(IO::DiscretizationReader& reader)
{
  Teuchos::RCP<Epetra_Vector> lambdafull =
      Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(), true));
  reader.ReadVector(lambdafull, "fsilambda");
  InsertVector(0, lambdafull, 0, lambda_, Coupling_Comm_Manager::full_to_partial);
  return;
}

/*-----------------------------------------------------------------------------------------*
| Get Timeface on the interface (for OST this is 1/(theta dt))                ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
double XFEM::XFSCoupling_Manager::GetInterfaceTimefac()
{
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
