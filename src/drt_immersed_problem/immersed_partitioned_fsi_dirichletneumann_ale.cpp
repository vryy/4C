/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann_ale.cpp

\brief partitioned immersed fsi algorithm for neumann-neumann like coupling (volume force coupling)

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_dirichletneumann_ale.H"

#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_field_wrapper.H"
#include "../drt_adapter/ad_fld_wrapper.H"
#include "../drt_adapter/ad_ale_wrapper.H"

#include "../drt_adapter/ad_fld_fluid_ale_immersed.H"

#include "../drt_structure/stru_aux.H"


IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ImmersedPartitionedFSIDirichletNeumannALE(const Epetra_Comm& comm)
  : ImmersedPartitionedFSIDirichletNeumann(comm),
    idispnp_(Teuchos::null),
    ivelnp_(Teuchos::null)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  SetupCoupling(fsidyn,comm);

  /// Immersed+FSI interface vector of new state to use in residual calculation
  combined_newstate_ = Teuchos::rcp(new Epetra_Vector(*(immersedstructure_->CombinedInterface()->FullMap()),true));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetupCoupling(const Teuchos::ParameterList& fsidyn ,const Epetra_Comm& comm)
{
  std::cout<<" Setup ALE coupling at FSI Interface ..."<<std::endl;
  FSI::Partitioned::SetupCoupling(fsidyn,comm);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetStatesFluidOP()
{
  // for FluidOP
  structdis_->SetState(0,"displacement",immersedstructure_->Dispnp());
  structdis_->SetState(0,"velocity",immersedstructure_->Velnp());
  fluiddis_ ->SetState(0,"dispnp",MBFluidField()->FluidField()->Dispnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetStatesVelocityCorrection()
{
  fluiddis_->SetState(0,"velnp",MBFluidField()->FluidField()->Velnp());
  fluiddis_->SetState(0,"dispnp",MBFluidField()->FluidField()->Dispnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SetStatesStructOP()
{
  // for StructOP
  fluiddis_->SetState(0,"velnp",MBFluidField()->FluidField()->Velnp());
  fluiddis_ ->SetState(0,"dispnp",MBFluidField()->FluidField()->Dispnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::InitialGuess()
{

  if(displacementcoupling_)
    return immersedstructure_->PredictFullInterfaceDispnp();
  else
  {

    //struct_bdry_traction_old_ = Teuchos::rcp(new Epetra_Vector(*(immersedstructure_->DofRowMap()),true));
    Teuchos::RCP<Epetra_Vector> combined_traction = Teuchos::rcp(new Epetra_Vector(*(immersedstructure_->CombinedInterface()->FullMap()),true));

    // insert Immersed and FSI vector to combined vector (vector of whole interface)
    immersedstructure_->CombinedInterface()->InsertOtherVector(immersedstructure_->Interface()->ExtractIMMERSEDCondVector(struct_bdry_traction_),combined_traction);
    immersedstructure_->CombinedInterface()->InsertCondVector(FluidToStruct(MBFluidField()->ExtractInterfaceForces()),combined_traction);

    return combined_traction;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::FluidOp(Teuchos::RCP<Epetra_Vector> bforce,
    const FillType fillFlag)
{
  // the displacement -> velocity conversion at the ALE FSI interface
  idispnp_ = immersedstructure_->ExtractInterfaceDispnp();
  ivelnp_ = InterfaceVelocity(idispnp_);

  // immersed part
  IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FluidOp(bforce,fillFlag);

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::SolveFluid()
{
  MBFluidField()->NonlinearSolve(StructToFluid(idispnp_),StructToFluid(ivelnp_));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::StructOp(Teuchos::RCP<Epetra_Vector> struct_bdry_traction,
    const FillType fillFlag)
{
  if(fillFlag==User)
  {
    dserror("fillFlag==User not implemented");
    return Teuchos::null;
  }
  else
  {
    // immersed part
    IMMERSED::ImmersedPartitionedFSIDirichletNeumann::StructOp(struct_bdry_traction,fillFlag);

    // add FSI part
    immersedstructure_->CombinedInterface()->AddCondVector(FluidToStruct(MBFluidField()->ExtractInterfaceForces()),struct_bdry_traction);

    // solve the structure
    IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SolveStruct();

    return ExtractInterfaceDispnp();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::UpdateCurrentPositionsFluidNodes()
{
  // get displacement state
    Teuchos::RCP<const Epetra_Vector> displacements = fluid_->FluidField()->Dispnp();

   // update positions
   for (int lid = 0; lid < fluiddis_->NumMyColNodes(); ++lid)
   {
     const DRT::Node* node = fluiddis_->lColNode(lid);
     LINALG::Matrix<3,1> currpos;
     std::vector<int> dofstoextract(4);
     std::vector<double> mydisp(4);

     // get the current displacement
     fluiddis_->Dof(node,0,dofstoextract);
     DRT::UTILS::ExtractMyValues(*displacements,mydisp,dofstoextract);

     currpos(0) = node->X()[0]+mydisp.at(0);
     currpos(1) = node->X()[1]+mydisp.at(1);
     currpos(2) = node->X()[2]+mydisp.at(2);

     currpositions_fluid_[node->Id()] = currpos;
   }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ExtractInterfaceDispnp()
{
  return immersedstructure_->ExtractFullInterfaceDispnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> full_traction_vec)
{
  immersedstructure_->ApplyImmersedInterfaceForces(FluidToStruct(MBFluidField()->ExtractInterfaceForces()),immersedstructure_->CombinedInterface()->ExtractOtherVector(*full_traction_vec));

  return;
}

void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ExtractPreviousInterfaceSolution()
{
  FSI::Partitioned::ExtractPreviousInterfaceSolution();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::AddDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidAleImmersed >(MBFluidField())->AddDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::RemoveDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidAleImmersed >(MBFluidField())->RemoveDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::CalcResidual(Epetra_Vector &F, const Teuchos::RCP<Epetra_Vector> newstate , const Teuchos::RCP<Epetra_Vector> oldstate)
{
  // insert Immersed and FSI vector to combined vector (vector of whole interface)
  immersedstructure_->CombinedInterface()->InsertOtherVector(immersedstructure_->Interface()->ExtractIMMERSEDCondVector(newstate),combined_newstate_);
  immersedstructure_->CombinedInterface()->InsertCondVector(immersedstructure_->Interface()->ExtractFSICondVector(newstate),combined_newstate_);

  int err = F.Update(1.0, *combined_newstate_, -1.0, *oldstate, 0.0);


  return err;
}
