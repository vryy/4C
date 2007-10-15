
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_fluid.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



FSI::Fluid::Fluid(RefCountPtr<DRT::Discretization> dis,
                  RefCountPtr<LINALG::Solver> solver,
                  RefCountPtr<ParameterList> params,
                  RefCountPtr<IO::DiscretizationWriter> output)
  : FluidImplicitTimeInt(dis,*solver,*params,*output,true),
    solver_(solver),
    params_(params),
    output_(output)
{
  stepmax_    = params_->get<int>   ("max number timesteps");
  maxtime_    = params_->get<double>("total time");
  theta_      = params_->get<double>("theta");
  timealgo_   = params_->get<FLUID_TIMEINTTYPE>("time int algo");
  dtp_ = dta_ = params_->get<double>("time step size");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  relax_ = LINALG::CreateVector(*dofrowmap,true);
  griddisp_ = LINALG::CreateVector(*dofrowmap,true);
}


int FSI::Fluid::Itemax() const
{
  return params_->get<int>("max nonlin iter steps");
}


void FSI::Fluid::SetItemax(int itemax)
{
  params_->set<int>("max nonlin iter steps", itemax);
}


void FSI::Fluid::SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im)
{
  ivelmap_ = im;
  extractor_ = rcp(new Epetra_Export(residual_->Map(), *ivelmap_));
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::Fluid::ExtractInterfaceForces()
{
  Teuchos::RefCountPtr<Epetra_Vector> iforce = rcp(new Epetra_Vector(*ivelmap_));
  int err = iforce->Export(*trueresidual_,*extractor_,Insert);
  if (err)
    dserror("Export using exporter returned err=%d",err);
  return iforce;
}


void FSI::Fluid::ApplyInterfaceVelocities(Teuchos::RefCountPtr<Epetra_Vector> ivel)
{
  int err = velnp_->Import(*ivel,*extractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  err = dirichtoggle_->Import(*ivel,*extractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);
}


void FSI::Fluid::SetMeshMap(Teuchos::RefCountPtr<Epetra_Map> mm)
{
  meshmap_ = mm;
  meshextractor_ = rcp(new Epetra_Export(residual_->Map(), *meshmap_));
}


void FSI::Fluid::ApplyMeshDisplacement(Teuchos::RefCountPtr<Epetra_Vector> fluiddisp)
{
  int err = dispnp_->Import(*fluiddisp,*meshextractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);

  // new grid velocity
  // There are other choices how to approximate that. Which one to
  // chose?
  gridv_->Update(1/dta_, *dispnp_, -1/dta_, *dispn_, 0.0);
}


void FSI::Fluid::ApplyMeshVelocity(Teuchos::RefCountPtr<Epetra_Vector> gridvel)
{
  //gridv_->PutScalar(0);
  int err = gridv_->Import(*gridvel,*meshextractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::Fluid::RelaxationSolve(Teuchos::RefCountPtr<Epetra_Vector> ivel)
{
  //
  // This method is really stupid, but simple. We calculate the fluid
  // elements twice here. First because we need the global matrix to
  // do a linear solve. We do not have any RHS other than the one from
  // the Dirichlet condition at the FSI interface.
  //
  // After that we recalculate the matrix so that we can get the
  // reaction forces at the FSI interface easily.
  //
  // We do more work that required, but we do not need any special
  // element code to perform the steepest descent calculation. This is
  // quite a benefit as the special code in the old discretization was
  // a real nightmare.
  //

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  relax_->PutScalar(0.0);
  int err = relax_->Import(*ivel,*extractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);

  // set the grid displacement independent of the trial value at the
  // interface
  griddisp_->Update(1., *dispnp_, -1., *dispn_, 0.);

  // dirichtoggle_ has already been set up

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  ParameterList eleparams;
  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("total time",time_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("using stationary formulation",false);
  eleparams.set("include reactive terms for linearisation",newton_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist"  ,zeros_);
#if 0
  discret_->SetState("dispnp", zeros_);
  discret_->SetState("gridv", zeros_);
#else
  discret_->SetState("dispnp", griddisp_);
  discret_->SetState("gridv", gridv_);
#endif

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // finalize the system matrix
  LINALG::Complete(*sysmat_);

  // No, we do not want to have any rhs. There cannot be any.
  residual_->PutScalar(0.0);

  //--------- Apply dirichlet boundary conditions to system of equations
  //          residual discplacements are supposed to be zero at
  //          boundary conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,relax_,dirichtoggle_);

  //-------solve for residual displacements to correct incremental displacements
  solver_->Solve(sysmat_,incvel_,residual_,true,true);

  // and now we need the reaction forces

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("using stationary formulation",false);

  // set vector values needed by elements
  discret_->ClearState();
  //discret_->SetState("velnp",incvel_);
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist"  ,zeros_);
#if 0
  discret_->SetState("dispnp", zeros_);
  discret_->SetState("gridv", zeros_);
#else
  discret_->SetState("dispnp", griddisp_);
  discret_->SetState("gridv", gridv_);
#endif

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  double density = eleparams.get("density", 0.0);

//   double norm;
//   trueresidual_->Norm2(&norm);
//   if (trueresidual_->Map().Comm().MyPID()==0)
//     cout << "==> residual norm = " << norm << " <==\n";

  // finalize the system matrix
  LINALG::Complete(*sysmat_);

  if (sysmat_->Apply(*incvel_, *trueresidual_)!=0)
    dserror("sysmat_->Apply() failed");
  trueresidual_->Scale(-density/dta_/theta_);

  return ExtractInterfaceForces();
}


RefCountPtr<Epetra_Vector> FSI::Fluid::IntegrateInterfaceShape()
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set("action","integrate_Shapefunction");

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  RefCountPtr<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

  // call loop over elements
  discret_->ClearState();
  discret_->SetState("dispnp", dispnp_);
  discret_->EvaluateCondition(eleparams,*integratedshapefunc,"FSICoupling");
  discret_->ClearState();

  Teuchos::RefCountPtr<Epetra_Vector> ishape = rcp(new Epetra_Vector(*ivelmap_));
  int err = ishape->Export(*integratedshapefunc,*extractor_,Insert);
  if (err)
    dserror("Export using exporter returned err=%d",err);
  return ishape;
}

#endif
#endif
