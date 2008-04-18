/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_xfem.cpp

\brief 

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "adapter_fluid_xfem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(const Teuchos::ParameterList& prbdyn,
                                          std::string condname)
  : fluid_(prbdyn,false),
    boundary_(prbdyn,false)
{
//  icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
//                                   FluidField().Interface(),
//                                  *AleField().Discretization(),
//                                   AleField().Interface(),
//                                   condname);

  //FSI::Coupling& coupfa = FluidAleFieldCoupling();

  // the fluid-ale coupling always matches
  //const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  //const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

//  coupfa_.SetupCoupling(*FluidField().Discretization(),
//                        *AleField().Discretization(),
//                        *fluidnodemap,
//                        *alenodemap);

  //FluidField().SetMeshMap(coupfa_.MasterDofMap());

  // the ale matrix is build just once
  //AleField().BuildSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::Discretization()
{
  return FluidField().Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const LINALG::MapExtractor& ADAPTER::FluidXFEM::Interface() const
{
  return FluidField().Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::PrepareTimeStep()
{
  FluidField().PrepareTimeStep();
  //AleField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update()
{
  FluidField().Update();
  //AleField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output()
{
  FluidField().Output();
  //AleField().Output();

  FluidField().LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFEM::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  //AleField().ReadRestart(step);
  return FluidField().Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                          Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    //AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));
    FluidField().ApplyInterfaceVelocities(ivel);
  }

  //if (FluidField().FreeSurface().Relevant())
  //{
  //  Teuchos::RCP<const Epetra_Vector> dispnp = FluidField().Dispnp();
  //  Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField().FreeSurface().ExtractCondVector(dispnp);
  //  AleField().ApplyFreeSurfaceDisplacements(fscoupfa_.MasterToSlave(fsdispnp));
  //}

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  //AleField().Solve();
  //Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  //FluidField().ApplyMeshDisplacement(fluiddisp);
  FluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                                      double dt)
{
  // Here we have a mesh position independent of the
  // given trial vector, but still the grid velocity depends on the
  // trial vector only.

  // grid velocity
  //AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));

  //AleField().Solve();
  //Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  //fluiddisp->Scale(1./dt);

  //FluidField().ApplyMeshVelocity(fluiddisp);

  // grid position is done inside RelaxationSolve

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1./dt);

  return FluidField().RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceForces()
{
  return FluidField().ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::IntegrateInterfaceShape()
{
  return FluidField().IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidXFEM::CreateFieldTest()
{
  return FluidField().CreateFieldTest();
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
//{
//  return coupfa_.SlaveToMaster(iv);
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const
//{
//  return coupfa_.SlaveToMaster(iv);
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
//{
//  return icoupfa_.MasterToSlave(iv);
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const
//{
//  return icoupfa_.MasterToSlave(iv);
//}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ADAPTER::XFluidFSIBoundary::XFluidFSIBoundary(RCP<DRT::Discretization> actdis,
                              Teuchos::RCP<LINALG::Solver> solver,
                              Teuchos::RCP<ParameterList> params,
                              Teuchos::RCP<IO::DiscretizationWriter> output,
                              int aletype,
                              bool dirichletcond)
  : discret_(actdis),
    solver_ (solver),
    params_ (params),
    output_ (output),
    step_(0),
    time_(0.0),
    aletype_(aletype),
    sysmat_(null),
    restartstep_(0),
    uprestart_(params->get("write restart every", -1))
{
  numstep_ = params_->get<int>("numstep");
  maxtime_ = params_->get<double>("maxtime");
  dt_      = params_->get<double>("dt");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispn_          = LINALG::CreateVector(*dofrowmap,true);
  dispnp_         = LINALG::CreateVector(*dofrowmap,true);
  residual_       = LINALG::CreateVector(*dofrowmap,true);
  dirichtoggle_   = LINALG::CreateVector(*dofrowmap,true);

  UTILS::SetupNDimExtractor(*actdis,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*actdis,"FREESURFCoupling",freesurface_);

  // set fixed nodes (conditions != 0 are not supported right now)
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,dirichtoggle_);

  if (dirichletcond)
  {
    // for partitioned FSI the interface becomes a Dirichlet boundary

    Teuchos::RCP<Epetra_Vector> idisp = LINALG::CreateVector(*interface_.CondMap(),false);
    idisp->PutScalar(1.0);
    interface_.InsertCondVector(idisp,dirichtoggle_);
  }

  if (dirichletcond and freesurface_.Relevant())
  {
    // for partitioned solves the free surface becomes a Dirichlet boundary

    Teuchos::RCP<Epetra_Vector> idisp = LINALG::CreateVector(*freesurface_.CondMap(),false);
    idisp->PutScalar(1.0);
    freesurface_.InsertCondVector(idisp,dirichtoggle_);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::BuildSystemMatrix(bool full)
{
  // build linear matrix once and for all
  if (full)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,false,true));
  }
  else
  {
    sysmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(interface_,interface_,81,false,true));
  }

  EvaluateElements();
  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const LINALG::SparseMatrix* ADAPTER::XFluidFSIBoundary::InteriorMatrixBlock() const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* bm =
    dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&*sysmat_);
  if (bm!=NULL)
  {
    return &bm->Matrix(0,0);
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const LINALG::SparseMatrix* ADAPTER::XFluidFSIBoundary::InterfaceMatrixBlock() const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* bm =
    dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&*sysmat_);
  if (bm!=NULL)
  {
    return &bm->Matrix(0,1);
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp) const
{
  // We save the current solution here. This will not change the
  // result of our element call, but the next time somebody asks us we
  // know the displacements.
  //
  // Note: What we get here is the sum of all increments in this time
  // step, not just the latest increment. Be careful.

  if (ddisp!=Teuchos::null)
  {
    // Dirichlet boundaries != 0 are not supported.

    dispnp_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::Solve()
{
  if (aletype_==ALE_DYNAMIC::springs)
    EvaluateElements();

  // set fixed nodes
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,dirichtoggle_);

  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);

  solver_->Solve(sysmat_->EpetraOperator(),dispnp_,residual_,true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::Update()
{
  dispn_->Update(1.0,*dispnp_,0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::Output()
{
  // We do not need any output -- the fluid writes its
  // displacements itself. But we need restart.

  restartstep_ += 1;

  output_->NewStep    (step_,time_);
  output_->WriteVector("dispnp", dispnp_);

  if (restartstep_ == uprestart_)
  {
    restartstep_ = 0;

    // add restart data
    output_->WriteVector("dispn", dispn_);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::Integrate()
{
  while (step_ < numstep_-1 and time_ <= maxtime_)
  {
    PrepareTimeStep();
    Solve();
    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::EvaluateElements()
{
  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // set vector values needed by elements
  discret_->ClearState();

  // action for elements
  if (aletype_==ALE_DYNAMIC::classic_lin)
  {
    eleparams.set("action", "calc_ale_lin_stiff");
  }
  else if (aletype_==ALE_DYNAMIC::springs)
  {
    discret_->SetState("dispnp", dispnp_);
    eleparams.set("action", "calc_ale_spring");
  }
  else
  {
    dserror("unsupported ale type");
  }

  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  sysmat_->Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_.InsertCondVector(idisp,dispnp_);

  // apply displacements to the rhs as well
  interface_.InsertCondVector(idisp,residual_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSIBoundary::ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp)
{
  freesurface_.InsertCondVector(fsdisp,dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSIBoundary::ExtractDisplacement() const
{
  // We know that the ale dofs are coupled with their original map. So
  // we just return them here.
  return dispnp_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSIBoundary::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(dispnp_);
// }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double ADAPTER::XFluidFSIBoundary::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(dispnp_, "dispnp");
  return time_;
}

#endif
