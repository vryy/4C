/*----------------------------------------------------------------------*/
/*! \file
\brief Coupling Manager for eXtended Fluid Fluid Coupling

\level 3

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

*----------------------------------------------------------------------*/
#include "XFFcoupling_manager.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"


#include "../drt_lib/drt_discret_xfem.H"
//#include "../drt_io/io.H"

/*-----------------------------------------------------------------------------------------*
| Constructor                                                                 ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
XFEM::XFFCoupling_Manager::XFFCoupling_Manager(Teuchos::RCP<ConditionManager> condmanager,
    Teuchos::RCP<FLD::XFluid> xfluid, Teuchos::RCP<FLD::XFluid> fluid, std::vector<int> idx)
    : Coupling_Comm_Manager(fluid->Discretization(), "XFEMSurfFluidFluid", 0, 3),
      fluid_(fluid),
      xfluid_(xfluid),
      cond_name_("XFEMSurfFluidFluid"),
      idx_(idx)
{
  if (idx_.size() != 2) dserror("XFFCoupling_Manager required two block ( 2 != %d)", idx_.size());

  // Coupling_Comm_Manager create all Coupling Objects now with Fluid has idx = 0, Fluid has idx =
  // 1!
  mcffi_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFluidFluid>(
      condmanager->GetMeshCoupling(cond_name_));
  if (mcffi_ == Teuchos::null) dserror(" Failed to get MeshCouplingFFI for embedded fluid!");
}

/*-----------------------------------------------------------------------------------------*
| Set required displacement & velocity states in the coupling object          ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFFCoupling_Manager::InitCouplingStates() {}


/*-----------------------------------------------------------------------------------------*
| Set required displacement & velocity states in the coupling object          ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFFCoupling_Manager::SetCouplingStates()
{
  std::cout << "SetCouplingStates in XFFCoupling_Manager" << std::endl;

  /// free the fluid-fluid interface
  mcffi_->SetInterfaceFree();

  mcffi_->UpdateDisplacementIterationVectors();  // update last iteration interface displacements
  LINALG::Export(*fluid_->Dispnp(), *mcffi_->IDispnp());
  LINALG::Export(*fluid_->Velnp(), *mcffi_->IVelnp());
  LINALG::Export(*fluid_->Veln(), *mcffi_->IVeln());

  Teuchos::RCP<Epetra_Vector> tmp_diff =
      Teuchos::rcp(new Epetra_Vector((*mcffi_->IDispnp()).Map()));
  tmp_diff->Update(1.0, *mcffi_->IDispnp(), -1.0, *mcffi_->IDispnpi(), 0.0);

  double norm = 0.0;
  tmp_diff->NormInf(&norm);

  if (norm < 1e-12)
    std::cout << "No change in XFF interface position!!!" << std::endl;
  else
  {
    std::cout << "Change in XFF interface position??? with infnorm " << norm << std::endl;
  }


  //  std::cout << "mcffi-IDispnp()" << *mcffi_->IDispnp()  << std::endl;
  //  std::cout << "mcffi-IVelpnp()" << *mcffi_->IVelnp()   << std::endl;

  //  //1 update last increment, before we set new idispnp
  //  mcffi_->UpdateDisplacementIterationVectors();
  //
  //  //2 Set Displacement on both mesh couplings ... we get them from the embedded fluid field!
  //  InsertVector(0,fluid_->Dispnp(),0,mcffi_->IDispnp(),Coupling_Comm_Manager::full_to_partial);
  //
  //
  //  InsertVector(0,fluid_->Velnp(),0,mcffi_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
  //

  return;
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling matrixes to the global systemmatrix                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFFCoupling_Manager::AddCouplingMatrix(
    LINALG::BlockSparseMatrixBase& systemmatrix, double scaling)
{
  /*----------------------------------------------------------------------*/
  // Coupling blocks C_fxf, C_xff and C_ff
  /*----------------------------------------------------------------------*/
  LINALG::SparseMatrix& C_ff_block = (systemmatrix)(idx_[0], idx_[0]);
  /*----------------------------------------------------------------------*/

  // add the coupling block C_ss on the already existing diagonal block
  C_ff_block.Add(*xfluid_->C_ss_Matrix(cond_name_), false, scaling, 1.0);

  LINALG::SparseMatrix& C_xff_block = (systemmatrix)(idx_[1], idx_[0]);
  LINALG::SparseMatrix& C_fxf_block = (systemmatrix)(idx_[0], idx_[1]);

  C_fxf_block.Add(*xfluid_->C_sx_Matrix(cond_name_), false, scaling, 1.0);
  C_xff_block.Add(*xfluid_->C_xs_Matrix(cond_name_), false, scaling, 1.0);
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling rhs                                                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFFCoupling_Manager::AddCouplingRHS(
    Teuchos::RCP<Epetra_Vector> rhs, const LINALG::MultiMapExtractor& me, double scaling)
{
  // REMARK: Copy this vector to store the correct lambda_ in update!
  Teuchos::RCP<Epetra_Vector> coup_rhs_sum =
      Teuchos::rcp(new Epetra_Vector(*xfluid_->RHS_s_Vec(cond_name_)));

  coup_rhs_sum->Scale(scaling);

  Teuchos::RCP<Epetra_Vector> coup_rhs = Teuchos::rcp(new Epetra_Vector(*me.Map(idx_[0]), true));
  LINALG::Export(*coup_rhs_sum, *coup_rhs);
  me.AddVector(coup_rhs, idx_[0], rhs);

  return;
}
