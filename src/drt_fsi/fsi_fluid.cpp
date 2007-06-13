
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
                  RefCountPtr<DiscretizationWriter> output)
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


#endif
#endif
