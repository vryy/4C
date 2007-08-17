
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


Teuchos::RefCountPtr<Epetra_Vector> FSI::Fluid::RelaxationSolve(Teuchos::RefCountPtr<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  relax_->PutScalar(0.0);
  int err = relax_->Import(*ivel,*extractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);

  // dirichtoggle_ has already been set up

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  ParameterList eleparams;
  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("time constant for integration",theta_*dta_);
  eleparams.set("using stationary formulation",false);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p at time n+1 (trial)",velnp_);
  discret_->SetState("old solution data for rhs"  ,zeros_);
  discret_->SetState("dispnp", zeros_);
  discret_->SetState("gridv", zeros_);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // finalize the system matrix
  LINALG::Complete(*sysmat_);

  //--------- Apply dirichlet boundary conditions to system of equations
  //          residual discplacements are supposed to be zero at
  //          boundary conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,relax_,dirichtoggle_);

  //-------solve for residual displacements to correct incremental displacements
  solver_->Solve(sysmat_,incvel_,residual_,true,true);

  // and now we need the residuum

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("time constant for integration",theta_*dta_);
  eleparams.set("using stationary formulation",false);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p at time n+1 (trial)",incvel_);
  discret_->SetState("old solution data for rhs"  ,zeros_);
  discret_->SetState("dispnp", zeros_);
  discret_->SetState("gridv", zeros_);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  double density = eleparams.get("density", 0.0);
  trueresidual_->Update(density/dta_/theta_,*residual_,0.0);

  return ExtractInterfaceForces();
}


#endif
#endif
