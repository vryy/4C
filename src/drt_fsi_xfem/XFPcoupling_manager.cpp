/*!----------------------------------------------------------------------
\file XFPcoupling_manager.cpp
\brief  Coupling Manager for eXtended Fluid Poro Coupling

\level 3

<pre>
\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>

*----------------------------------------------------------------------*/
#include "XFPcoupling_manager.H"

#include "../drt_xfem/xfem_condition_manager.H"

#include "../drt_poroelast/poro_base.H"

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"

#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_io/io.H"

//remove me
#include "../drt_adapter/adapter_coupling.H"

XFEM::XFPCoupling_Manager::XFPCoupling_Manager(Teuchos::RCP<XFEM::ConditionManager> condmanager,
    Teuchos::RCP<POROELAST::PoroBase> poro, Teuchos::RCP<FLD::XFluid> xfluid, std::vector<int> idx)
:Coupling_Comm_Manager(poro->StructureField()->Discretization(),poro->FluidField()->Discretization(),"XFEMSurfFPIMono",0,3),
 poro_(poro),
 xfluid_(xfluid),
 cond_name_ps_ps_("XFEMSurfFPIMono_ps_ps"),
 cond_name_ps_pf_("XFEMSurfFPIMono_ps_pf"),
 cond_name_pf_ps_("XFEMSurfFPIMono_pf_ps"),
 cond_name_pf_pf_("XFEMSurfFPIMono_pf_pf"),
 idx_(idx)
{
  //Coupling_Comm_Manager create all Coupling Objects now with Structure has idx = 0, Fluid has idx = 1!

  mcfpi_ps_ps_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(condmanager->GetMeshCoupling(cond_name_ps_ps_));
  if (mcfpi_ps_ps_ == Teuchos::null) dserror(" Failed to get MeshCouplingFPI for Porostructure!");
  mcfpi_ps_ps_->InitializeStrucPresMap(poro_->FluidStructureCoupling().SlaveDofMap(),poro_->FluidStructureCoupling().PermMasterDofMap());

  mcfpi_ps_pf_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(condmanager->GetMeshCoupling(cond_name_ps_pf_));
  if (mcfpi_ps_pf_ == Teuchos::null) dserror(" Failed to get MeshCouplingFPI for Porofluid!");
  mcfpi_ps_pf_->InitializeStrucPresMap(poro_->FluidStructureCoupling().SlaveDofMap(),poro_->FluidStructureCoupling().PermMasterDofMap());

  mcfpi_pf_ps_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(condmanager->GetMeshCoupling(cond_name_pf_ps_));
  if (mcfpi_pf_ps_ == Teuchos::null) dserror(" Failed to get MeshCouplingFPI for Porofluid!");
  mcfpi_pf_ps_->InitializeStrucPresMap(poro_->FluidStructureCoupling().SlaveDofMap(),poro_->FluidStructureCoupling().PermMasterDofMap());

  mcfpi_pf_pf_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(condmanager->GetMeshCoupling(cond_name_pf_pf_));
  if (mcfpi_pf_pf_ == Teuchos::null) dserror(" Failed to get MeshCouplingFPI for Porofluid!");
  mcfpi_pf_pf_->InitializeStrucPresMap(poro_->FluidStructureCoupling().SlaveDofMap(),poro_->FluidStructureCoupling().PermMasterDofMap());

  //safety check
  if (!mcfpi_ps_ps_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    dserror("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (psps)!");
  if (!mcfpi_ps_pf_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    dserror("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (pspf)!");
  if (!mcfpi_pf_ps_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    dserror("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (pfps)!");
  if (!mcfpi_pf_pf_->IDispnp()->Map().SameAs(*GetMapExtractor(0)->Map(1)))
    dserror("XFPCoupling_Manager: Maps of Condition and Mesh Coupling do not fit (pfpf)!");

  // storage of the resulting Robin-type structural forces from the old timestep
  // Recovering of Lagrange multiplier happens on fluid field
  lambda_ps_ = Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->Map(1),true));
  lambda_pf_ = Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->Map(1),true));
}

void XFEM::XFPCoupling_Manager::SetCouplingStates()
{
  //1 Set Displacement on both mesh couplings ... we get them from the structure field!
  InsertVector(0,poro_->StructureField()->Dispnp(),0,mcfpi_ps_ps_->IDispnp(),Coupling_Comm_Manager::full_to_partial);
  InsertVector(0,poro_->StructureField()->Dispnp(),0,mcfpi_ps_pf_->IDispnp(),Coupling_Comm_Manager::full_to_partial);
  InsertVector(0,poro_->StructureField()->Dispnp(),0,mcfpi_pf_ps_->IDispnp(),Coupling_Comm_Manager::full_to_partial);
  InsertVector(0,poro_->StructureField()->Dispnp(),0,mcfpi_pf_pf_->IDispnp(),Coupling_Comm_Manager::full_to_partial);

  mcfpi_ps_ps_->SetFullDispnp(poro_->StructureField()->Dispnp(),poro_->FluidField()->ExtractPressurePart(poro_->FluidField()->Velnp()));
  mcfpi_ps_pf_->SetFullDispnp(poro_->StructureField()->Dispnp(),poro_->FluidField()->ExtractPressurePart(poro_->FluidField()->Velnp()));
  mcfpi_pf_ps_->SetFullDispnp(poro_->StructureField()->Dispnp(),poro_->FluidField()->ExtractPressurePart(poro_->FluidField()->Velnp()));
  mcfpi_pf_pf_->SetFullDispnp(poro_->StructureField()->Dispnp(),poro_->FluidField()->ExtractPressurePart(poro_->FluidField()->Velnp()));


  //2 Set Structural Velocity onto ps mesh coupling
  InsertVector(0,poro_->StructureField()->Velnp(),0,mcfpi_ps_ps_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
  InsertVector(0,poro_->StructureField()->Velnp(),0,mcfpi_pf_ps_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
//  poro_->StructureField()->Velnp()->Print(std::cout);

  //  InsertVector(1,poro_->FluidField()->GridVel(),0,mcfpi_ps_ps_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
//  InsertVector(1,poro_->FluidField()->GridVel(),0,mcfpi_pf_ps_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
  //3 Set Fluid Velocity onto pf mesh coupling
  InsertVector(1,poro_->FluidField()->Velnp(),0,mcfpi_ps_pf_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
  InsertVector(1,poro_->FluidField()->Velnp(),0,mcfpi_pf_pf_->IVelnp(),Coupling_Comm_Manager::full_to_partial);
}

void XFEM::XFPCoupling_Manager::AddCouplingMatrix(LINALG::BlockSparseMatrixBase& systemmatrix, double scaling)
{
  if (idx_.size() == 2) //assum that the poro field is not split and we just have a blockmatrix P/F
  {
    const double scaling_disp_vel = 1/((1-poro_->StructureField()->TimIntParam())*poro_->StructureField()->Dt());
    const double dt = poro_->FluidField()->Dt();

    LINALG::SparseMatrix& C_ss_block = (systemmatrix)(idx_[0],idx_[0]);
    LINALG::SparseMatrix& C_fs_block = (systemmatrix)(idx_[1],idx_[0]);
    LINALG::SparseMatrix& C_sf_block = (systemmatrix)(idx_[0],idx_[1]);

    //1// Add Blocks f-ps(2), ps-f(3), ps-ps(4)
    C_ss_block.Add(*xfluid_->C_ss_Matrix(cond_name_ps_ps_),false,scaling*scaling_disp_vel,1.0);
    C_sf_block.Add(*xfluid_->C_sx_Matrix(cond_name_ps_ps_),false,scaling,1.0);
    C_fs_block.Add(*xfluid_->C_xs_Matrix(cond_name_ps_ps_),false,scaling*scaling_disp_vel,1.0);

    //2// Add Blocks f-pf(5), ps-pf(6)
    Teuchos::RCP<LINALG::SparseMatrix> C_ps_pf = Teuchos::rcp(new LINALG::SparseMatrix(xfluid_->C_ss_Matrix(cond_name_ps_pf_)->RowMap(),81,false));
    InsertMatrix(-1,0,*xfluid_->C_ss_Matrix(cond_name_ps_pf_),1,*C_ps_pf,Coupling_Comm_Manager::col,1,true,false);
    C_ps_pf->Complete(*GetMapExtractor(1)->Map(1),*GetMapExtractor(0)->Map(1));
    Teuchos::RCP<LINALG::SparseMatrix> C_f_pf = Teuchos::rcp(new LINALG::SparseMatrix(xfluid_->C_xs_Matrix(cond_name_ps_pf_)->RowMap(),81,false));
    InsertMatrix(-1,0,*xfluid_->C_xs_Matrix(cond_name_ps_pf_),1,*C_f_pf,Coupling_Comm_Manager::col,1,true,false);
    C_f_pf->Complete(*GetMapExtractor(1)->Map(1),C_fs_block.RangeMap());
    C_fs_block.Add(*C_f_pf,false,scaling,1.0);
    C_ss_block.Add(*C_ps_pf,false,scaling,1.0);

    //3// Add Blocks pf-f(7), pf-ps(8)
    Teuchos::RCP<LINALG::SparseMatrix> C_pf_ps = Teuchos::rcp(new LINALG::SparseMatrix(*GetMapExtractor(1)->Map(1),81,false));
    Teuchos::RCP<LINALG::SparseMatrix> C_pf_f = Teuchos::rcp(new LINALG::SparseMatrix(*GetMapExtractor(1)->Map(1),81,false));
    InsertMatrix(-1,0,*xfluid_->C_ss_Matrix(cond_name_pf_ps_),1,*C_pf_ps,Coupling_Comm_Manager::row,1,true,false);
    C_pf_ps->Complete(*GetMapExtractor(0)->Map(1),*GetMapExtractor(1)->Map(1));
    InsertMatrix(-1,0,*xfluid_->C_sx_Matrix(cond_name_pf_ps_),1,*C_pf_f,Coupling_Comm_Manager::row,1,true,false);
    C_pf_f->Complete(*xfluid_->DofRowMap(),*GetMapExtractor(1)->Map(1));
    C_ss_block.Add(*C_pf_ps,false,scaling*scaling_disp_vel*dt,1.0);
    C_sf_block.Add(*C_pf_f,false,scaling*dt,1.0);

    //4// Add Block pf-pf(9)
    Teuchos::RCP<LINALG::SparseMatrix> C_pf_pf = Teuchos::rcp(new LINALG::SparseMatrix(*GetMapExtractor(1)->Map(1),81,false));
    InsertMatrix(-1,0,*xfluid_->C_ss_Matrix(cond_name_pf_pf_),1,*C_pf_pf,Coupling_Comm_Manager::row_and_col);
    C_pf_pf->Complete(*GetMapExtractor(1)->Map(1),*GetMapExtractor(1)->Map(1));
    C_ss_block.Add(*C_pf_pf,false,scaling*dt,1.0);

  }
  else
    dserror("XFPCoupling_Manager::AddCouplingMatrix: Not implemented for number of blocks = %d", idx_.size());
}

///*-----------------------------------------------------------------------------------------*
//| Add the coupling rhs                                                        ager 06/2016 |
//*-----------------------------------------------------------------------------------------*/
void XFEM::XFPCoupling_Manager::AddCouplingRHS(Teuchos::RCP<Epetra_Vector> rhs,const LINALG::MultiMapExtractor& me, double scaling)
{
  if (idx_.size() == 2) //assum that the poro field is not split and we just have a blockmatrix P/F
  {

    const double dt = poro_->FluidField()->Dt();

      Teuchos::RCP<const Epetra_Vector> rhs_C_ps_ps = xfluid_->RHS_s_Vec(cond_name_ps_ps_);
      Teuchos::RCP<const Epetra_Vector> rhs_C_ps_pf = xfluid_->RHS_s_Vec(cond_name_ps_pf_);
      Teuchos::RCP<const Epetra_Vector> rhs_C_pf_ps = xfluid_->RHS_s_Vec(cond_name_pf_ps_);
      Teuchos::RCP<const Epetra_Vector> rhs_C_pf_pf = xfluid_->RHS_s_Vec(cond_name_pf_pf_);

      Teuchos::RCP<Epetra_Vector> prhs = Teuchos::rcp(new Epetra_Vector(*me.Map(idx_[0]),true));

      InsertVector(0,rhs_C_ps_ps,0,prhs,Coupling_Comm_Manager::partial_to_global,true,scaling);
      InsertVector(0,rhs_C_ps_pf,0,prhs,Coupling_Comm_Manager::partial_to_global,true,scaling);

      InsertVector(0,rhs_C_pf_ps,1,prhs,Coupling_Comm_Manager::partial_to_global,true,scaling*dt);
      InsertVector(0,rhs_C_pf_pf,1,prhs,Coupling_Comm_Manager::partial_to_global,true,scaling*dt);

      //Add lambda contribution
      if (lambda_ps_ != Teuchos::null && lambda_pf_ != Teuchos::null)
      {
        /*----------------------------------------------------------------------*/
        // get time integration parameters of structure and fluid time integrators
        // to enable consistent time integration among the fields
        /*----------------------------------------------------------------------*/

        /*----------------------------------------------------------------------*/
        // this is the interpolation weight for quantities from last time step
        // alpha_f for genalpha and (1-theta) for OST (weighting of the old time step n for displacements)
        // TimeIntegration for poro needs to be consistent!
        const double stiparam = poro_->StructureField()->TimIntParam();    // (1-theta) for OST and alpha_f for Genalpha

        // scale factor for the structure system matrix w.r.t the new time step
        const double scaling_S = 1.0/(1.0-stiparam);  // 1/(1-alpha_F) = 1/weight^S_np

       InsertVector(0,lambda_ps_,0,prhs,Coupling_Comm_Manager::partial_to_global,true,stiparam*scaling_S);
       InsertVector(0,lambda_pf_,1,prhs,Coupling_Comm_Manager::partial_to_global,true,stiparam*scaling_S);
      }

      me.AddVector(prhs,idx_[0],rhs);
  }
  else
    dserror("XFPCoupling_Manager::AddCouplingRHS: Not implemented for number of blocks = %d", idx_.size());
}

/*----------------------------------------------------------------------*/
/* Store the Coupling RHS of the Old Timestep in lambda     ager 06/2016 |
 *----------------------------------------------------------------------*/
 void XFEM::XFPCoupling_Manager::Update(double scaling)
 {
   /*----------------------------------------------------------------------*/
   // we directly store the fluid-unscaled rhs_C_s residual contribution from the fluid solver which corresponds to the actual acting forces

   // scaling for the structural residual is done when it is added to the global residual vector
   // get the coupling rhs from the xfluid, this vector is based on the boundary dis which is part of the structure dis
   lambda_ps_->Update(scaling,*xfluid_->RHS_s_Vec(cond_name_ps_ps_),0.0);
   lambda_ps_->Update(scaling,*xfluid_->RHS_s_Vec(cond_name_ps_pf_),1.0);

   const double dt = poro_->FluidField()->Dt();
   lambda_pf_->Update(scaling*dt,*xfluid_->RHS_s_Vec(cond_name_pf_ps_),0.0);
   lambda_pf_->Update(scaling*dt,*xfluid_->RHS_s_Vec(cond_name_pf_pf_),1.0);
   return;
}

/*----------------------------------------------------------------------*/
/* Write Output                                             ager 06/2016 |
*-----------------------------------------------------------------------*/
void XFEM::XFPCoupling_Manager::Output(IO::DiscretizationWriter& writer)
{
  //--------------------------------
  // output for Lagrange multiplier field (ie forces onto the structure, Robin-type forces
  // consisting of fluid forces and the Nitsche penalty term contribution)
  //--------------------------------
  Teuchos::RCP<Epetra_Vector> lambdafull =  Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(),true));
  InsertVector(0,lambda_ps_,0,lambdafull,Coupling_Comm_Manager::partial_to_full);
  writer.WriteVector("fpilambda_ps", lambdafull);

  lambdafull =  Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(),true));
  InsertVector(0,lambda_pf_,0,lambdafull,Coupling_Comm_Manager::partial_to_full);
  writer.WriteVector("fpilambda_pf", lambdafull);
 return;
}
/*----------------------------------------------------------------------*/
/* Read Restart on the interface                            ager 06/2016 |
*-----------------------------------------------------------------------*/
void XFEM::XFPCoupling_Manager::ReadRestart(IO::DiscretizationReader& reader)
{
  Teuchos::RCP<Epetra_Vector> lambdafull =  Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(),true));
  reader.ReadVector(lambdafull, "fpilambda_ps");
  InsertVector(0,lambdafull,0,lambda_ps_,Coupling_Comm_Manager::full_to_partial);

  lambdafull =  Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->FullMap(),true));
  reader.ReadVector(lambdafull, "fpilambda_pf");
  InsertVector(0,lambdafull,0,lambda_pf_,Coupling_Comm_Manager::full_to_partial);
  return;
}

/*-----------------------------------------------------------------------------------------*
| Get Timeface on the interface (for OST this is 1/(theta dt))                ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
double XFEM::XFPCoupling_Manager::GetInterfaceTimefac()
{
  dserror("Check if you really want this!");
  /*
   * Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double dt = xfluid_->Dt();
  if (interface_second_order_)
    return 2./dt;
  else
    return 1./dt;
}
