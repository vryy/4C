/*----------------------------------------------------------------------------*/
/*! \file
\brief Basic HDG weakly compressible time-integration scheme

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "fluid_timint_hdg_weak_comp.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_ele/fluid_ele_hdg_weak_comp.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/fluid_weakly_compressible.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                  laspina 08/19 |
 *----------------------------------------------------------------------*/
FLD::TimIntHDGWeakComp::TimIntHDGWeakComp(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntGenAlpha(actdis, solver, params, output, alefluid),
      timealgoset_(INPAR::FLUID::timeint_afgenalpha),
      firstAssembly_(false)
{
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                  laspina 08/19 |
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::Init()
{
  DRT::DiscretizationHDG* hdgdis = dynamic_cast<DRT::DiscretizationHDG*>(discret_.get());
  if (hdgdis == NULL) dserror("Did not receive an HDG discretization");

  // get number of spatial dimensions
  const unsigned int nsd = params_->get<int>("number of velocity degrees of freedom");

  // initialize density/momentum splitting
  std::vector<int> dof_all;
  std::set<int> dofset_r;
  std::set<int> dofset_w;

  // fill dofset
  for (int i = 0; i < hdgdis->NumMyRowFaces(); ++i)
  {
    dof_all = hdgdis->Dof(0, hdgdis->lRowFace(i));

    for (unsigned int j_r = 0; j_r < (dof_all.size() / (1 + nsd)); ++j_r)
      dofset_r.insert(dof_all[j_r]);

    for (unsigned int j_w = (dof_all.size() / (1 + nsd)); j_w < dof_all.size(); ++j_w)
      dofset_w.insert(dof_all[j_w]);
  }

  // define density dof map
  std::vector<int> dofmapvec_r;
  dofmapvec_r.reserve(dofset_r.size());
  dofmapvec_r.assign(dofset_r.begin(), dofset_r.end());
  dofset_r.clear();
  Teuchos::RCP<Epetra_Map> dofmap_r =
      Teuchos::rcp(new Epetra_Map(-1, dofmapvec_r.size(), &dofmapvec_r[0], 0, hdgdis->Comm()));

  // define momentum dof map
  std::vector<int> dofmapvec_w;
  dofmapvec_w.reserve(dofset_w.size());
  dofmapvec_w.assign(dofset_w.begin(), dofset_w.end());
  dofset_w.clear();
  Teuchos::RCP<Epetra_Map> dofmap_w =
      Teuchos::rcp(new Epetra_Map(-1, dofmapvec_w.size(), &dofmapvec_w[0], 0, hdgdis->Comm()));

  // build density/momentum (actually velocity/pressure) splitter
  velpressplitter_->Setup(*hdgdis->DofRowMap(), dofmap_r, dofmap_w);

  // implement ost and bdf2 through gen-alpha facilities
  if (timealgo_ == INPAR::FLUID::timeint_bdf2)
  {
    alphaM_ = 1.5;
    alphaF_ = 1.0;
    gamma_ = 1.0;
  }
  else if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)
  {
    alphaM_ = 1.0;
    alphaF_ = 1.0;
    gamma_ = params_->get<double>("theta");
  }
  else if (timealgo_ == INPAR::FLUID::timeint_stationary)
  {
    // mimic backward Euler neglecting inertial terms
    alphaM_ = 1.0;
    alphaF_ = 1.0;
    gamma_ = 1.0;
  }

  timealgoset_ = timealgo_;
  timealgo_ = INPAR::FLUID::timeint_afgenalpha;

  // call Init()-functions of base classes
  // note: this order is important
  FLD::TimIntGenAlpha::Init();
}



/*----------------------------------------------------------------------*
| calculate pseudo-theta for startalgo_, modified for HDG laspina 08/19 |
*-----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::SetTheta()
{
  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  // starting algorithm
  if (startalgo_ || (step_ <= 2 && timealgoset_ == INPAR::FLUID::timeint_bdf2))
  {
    // use backward-Euler-type parameter combination
    if (step_ <= numstasteps_ || (step_ <= 1 && timealgoset_ == INPAR::FLUID::timeint_bdf2))
    {
      if (myrank_ == 0)
      {
        std::cout << "Starting algorithm for Af_GenAlpha active. "
                  << "Performing step " << step_ << " of " << numstasteps_
                  << " Backward Euler starting steps" << std::endl;
      }
      alphaM_ = 1.0;
      alphaF_ = 1.0;
      gamma_ = 1.0;
    }
    else
    {
      // recall original user wish
      if (timealgoset_ == INPAR::FLUID::timeint_one_step_theta)
      {
        alphaM_ = alphaF_ = 1.0;
        gamma_ = params_->get<double>("theta");
      }
      else if (timealgoset_ == INPAR::FLUID::timeint_bdf2)
      {
        alphaF_ = gamma_ = 1.0;
        alphaM_ = 3. / 2.;
      }
      else
      {
        alphaM_ = params_->get<double>("alpha_M");
        alphaF_ = params_->get<double>("alpha_F");
        gamma_ = params_->get<double>("gamma");
      }

      // do not enter starting algorithm section in the future
      startalgo_ = false;
    }
  }

  // compute "pseudo-theta" for af-generalized-alpha scheme
  theta_ = alphaF_ * gamma_ / alphaM_;
}


/*----------------------------------------------------------------------*
 * Explicit predictor                                     laspina 08/19 |
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::ExplicitPredictor()
{
  if (predictor_ == "steady_state")
  {
    // this has already been done in TimeUpdate()
  }
  else if (predictor_ == "zero_acceleration")
  {
    velnp_->Update(1.0, *veln_, (1.0 - theta_) * dta_, *accn_, 0.0);
    intvelnp_->Update(1.0, *intveln_, (1.0 - theta_) * dta_, *intaccn_, 0.0);
  }
  else if (predictor_ == "constant_acceleration")
  {
    velnp_->Update(1.0, *veln_, dta_, *accn_, 0.0);
    intvelnp_->Update(1.0, *intveln_, dta_, *intaccn_, 0.0);
  }
  else if (predictor_ == "constant_increment")
  {
    velnp_->Update(2.0, *veln_, -1.0, *velnm_, 0.0);
    intvelnp_->Update(2.0, *intveln_, -1.0, *intvelnm_, 0.0);
  }
  else if (predictor_ == "explicit_second_order_midpoint")
  {
    velnp_->Update(1.0, *velnm_, 2.0 * dta_, *accn_, 0.0);
    intvelnp_->Update(1.0, *intvelnm_, 2.0 * dta_, *intaccn_, 0.0);
  }
  else
    dserror("Unknown fluid predictor %s", predictor_.c_str());
}

/*----------------------------------------------------------------------*
| set old part of right hand side                         laspina 08/19 |
*-----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::SetOldPartOfRighthandside()
{
  hist_->PutScalar(0.0);

  // This code is entered at the beginning of the nonlinear iteration, so
  // store that the assembly to be done next is going to be the first one
  // (without combined vector update) for HDG.
  firstAssembly_ = true;
}



/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha              laspina 08/19 |
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::GenAlphaUpdateAcceleration()
{
  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (i)
  //    acc    = acc * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // as opposed to standard fluid, update all variables

  // compute factors
  const double fact1 = 1.0 / (gamma_ * dta_);
  const double fact2 = 1.0 - (1.0 / gamma_);

  accnp_->Update(fact2, *accn_, 0.0);
  accnp_->Update(fact1, *velnp_, -fact1, *veln_, 1.0);

  intaccnp_->Update(fact2, *intaccn_, 0.0);
  intaccnp_->Update(fact1, *intvelnp_, -fact1, *intveln_, 1.0);
}



/*----------------------------------------------------------------------*
 | compute values at intermediate time steps              laspina 08/19 |
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::GenAlphaIntermediateValues()
{
  // set intermediate values for accelerations
  //
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  accam_->Update((alphaM_), *accnp_, (1.0 - alphaM_), *accn_, 0.0);
  intaccam_->Update((alphaM_), *intaccnp_, (1.0 - alphaM_), *intaccn_, 0.0);

  // set intermediate values for mixed variable, density and momentum
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  velaf_->Update((alphaF_), *velnp_, (1.0 - alphaF_), *veln_, 0.0);
  intvelaf_->Update((alphaF_), *intvelnp_, (1.0 - alphaF_), *intveln_, 0.0);
}


/*----------------------------------------------------------------------*
| set HDG state vectors                                   laspina 08/19 |
*-----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::SetStateTimInt()
{
  discret_->SetState(0, "velaf", velaf_);
  discret_->SetState(1, "intvelaf", intvelaf_);
  discret_->SetState(1, "intaccam", intaccam_);
  discret_->SetState(1, "intvelnp", intvelnp_);
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                   laspina 08/19 |
*-----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  eleparams.set<bool>("needslocalupdate", !firstAssembly_);
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                   laspina 08/19 |
*-----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::ClearStateAssembleMatAndRHS()
{
  if (!firstAssembly_)
  {
    // Wrote into the state vector during element calls, need to transfer the
    // data back before it disappears when clearing the state (at least for nproc>1)
    const Epetra_Vector& intvelnpGhosted = *discret_->GetState(1, "intvelnp");
    for (int i = 0; i < intvelnp_->MyLength(); ++i)
      (*intvelnp_)[i] = intvelnpGhosted[intvelnpGhosted.Map().LID(intvelnp_->Map().GID(i))];
  }
  firstAssembly_ = false;
  FluidImplicitTimeInt::ClearStateAssembleMatAndRHS();
}


/*----------------------------------------------------------------------*
| update within iteration                                 laspina 08/19 |
*-----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment)
{
  // call element routine to update local solution
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::update_local_solution);

  // location array
  DRT::Element::LocationArray la(2);

  // interior dofs map
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);

  // dummy variables
  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;

  // initialize elemental local increments
  Epetra_SerialDenseVector elemintinc;

  // initialize increments of local variables
  Teuchos::RCP<Epetra_Vector> intvelincnp = LINALG::CreateVector(*intdofrowmap, true);

  // set state
  SetStateTimInt();
  discret_->SetState(0, "globaltraceinc", increment);

  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // get element
    DRT::Element* ele = discret_->lColElement(el);
    ele->LocationVector(*discret_, la, false);

    // evaluate interior local increments
    ele->Evaluate(params, *discret_, la[0].lm_, dummyMat, dummyMat, elemintinc, dummyVec, dummyVec);

    // fill the interior increment vector for all the discretization
    if (ele->Owner() == discret_->Comm().MyPID())
    {
      std::vector<int> localDofs = discret_->Dof(1, ele);
      for (unsigned int i = 0; i < localDofs.size(); ++i)
        localDofs[i] = intdofrowmap->LID(localDofs[i]);
      intvelincnp->ReplaceMyValues(localDofs.size(), elemintinc.A(), &localDofs[0]);
    }
  }

  // update interior values by adding increments
  intvelnp_->Update(1.0, *intvelincnp, 1.0);

  // set new state
  SetStateTimInt();

  // call base function
  FluidImplicitTimeInt::IterUpdate(increment);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::TimeUpdate()
{
  FluidImplicitTimeInt::TimeUpdate();

  // local solution of this step become most recent
  // local solution of the last step
  intvelnm_->Update(1.0, *intveln_, 0.0);
  intveln_->Update(1.0, *intvelnp_, 0.0);

  intaccnm_->Update(1.0, *intaccn_, 0.0);
  intaccn_->Update(1.0, *intaccnp_, 0.0);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::UpdateGridv()
{
  if (timealgoset_ == INPAR::FLUID::timeint_afgenalpha ||
      timealgoset_ == INPAR::FLUID::timeint_npgenalpha ||
      timealgoset_ == INPAR::FLUID::timeint_bdf2)  // 2nd order methods
  {
    if (step_ <= 1)
    {
      // use BDF1 in 1st step
      gridv_->Update(1.0 / dta_, *dispnp_, -1.0 / dta_, *dispn_, 0.0);
    }
    else
    {
      // use BDF2 after 1st step
      gridv_->Update(1.5 / dta_, *dispnp_, -2.0 / dta_, *dispn_, 0.0);
      gridv_->Update(0.5 / dta_, *dispnm_, 1.0);
    }
  }
  else if (timealgoset_ == INPAR::FLUID::timeint_one_step_theta)  // 1st order methods
  {
    // use BDF1
    gridv_->Update(1.0 / dta_, *dispnp_, -1.0 / dta_, *dispn_, 0.0);
  }
  else if (timealgoset_ == INPAR::FLUID::timeint_stationary)
  {
    gridv_->Scale(0.0);
  }
}



/*----------------------------------------------------------------------*
 |  set initial flow field                                 laspina 08/19|
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield, const int startfuncno)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;
  initParams.set<int>("action", FLD::project_fluid_field);
  initParams.set("startfuncno", startfuncno);
  initParams.set<int>("initfield", initfield);

  // loop over all elements on the processor
  DRT::Element::LocationArray la(2);
  double error = 0;
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    DRT::Element* ele = discret_->lColElement(el);

    ele->LocationVector(*discret_, la, false);
    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1, ele)) elevec2.Shape(discret_->NumDof(1, ele), 1);

    ele->Evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);

    // now fill the ele vector into the discretization
    for (unsigned int i = 0; i < la[0].lm_.size(); ++i)
    {
      const int lid = dofrowmap->LID(la[0].lm_[i]);
      if (lid >= 0)
      {
        if ((*velnp_)[lid] != 0) error += std::abs((*velnp_)[lid] - elevec1(i));
        (*velnp_)[lid] = elevec1(i);
        (*veln_)[lid] = elevec1(i);
        (*velnm_)[lid] = elevec1(i);
      }
    }

    if (ele->Owner() == discret_->Comm().MyPID())
    {
      std::vector<int> localDofs = discret_->Dof(1, ele);
      dsassert(localDofs.size() == static_cast<std::size_t>(elevec2.M()), "Internal error");
      for (unsigned int i = 0; i < localDofs.size(); ++i)
        localDofs[i] = intdofrowmap->LID(localDofs[i]);
      intvelnp_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
      intveln_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
      intvelnm_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
    }
  }

  double globerror = 0;
  discret_->Comm().SumAll(&error, &globerror, 1);
  if (discret_->Comm().MyPID() == 0)
    std::cout << "Error project when setting face twice: " << globerror << std::endl;
}



/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions laspina 08/19|
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> FLD::TimIntHDGWeakComp::EvaluateErrorComparedToAnalyticalSol()
{
  // HDG needs one more state vector for the interior solution (i.e., the actual solution)
  INPAR::FLUID::CalcError calcerr =
      DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_, "calculate error");

  switch (calcerr)
  {
    case INPAR::FLUID::no_error_calculation:
    {
      return Teuchos::null;
      break;
    }
    case INPAR::FLUID::byfunct:
    {
      discret_->SetState(1, "intvelnp", intvelnp_);

      // std::vector containing
      // [0]: absolute L2 mixed variable error
      // [1]: absolute L2 density error
      // [2]: absolute L2 momentum error
      Teuchos::RCP<std::vector<double>> abserror = Teuchos::rcp(new std::vector<double>(3));

      // create the parameters for the discretization
      Teuchos::ParameterList eleparams;

      // action for elements
      eleparams.set<int>("action", FLD::calc_fluid_error);
      eleparams.set<int>("Physical Type", physicaltype_);
      eleparams.set<int>("calculate error", calcerr);

      // get function number
      const int errorfunctno = params_->get<int>("error function number", -1);
      eleparams.set<int>("error function number", errorfunctno);

      // set scheme-specific element parameters and vector values
      SetStateTimInt();

      if (alefluid_) discret_->SetState(2, "dispnp", dispnp_);

      // get (squared) error values
      // 0: delta mixed variable for L2-error norm
      // 1: delta density for L2-error norm
      // 2: delta momentum for L2-error norm
      // (3: analytical mixed variable for L2 norm)
      // (4: analytical density for L2 norm)
      // (5: analytical momentum for L2 norm)
      Teuchos::RCP<Epetra_SerialDenseVector> errors =
          Teuchos::rcp(new Epetra_SerialDenseVector(3 + 3));

      // call loop over elements (assemble nothing)
      discret_->EvaluateScalars(eleparams, errors);
      discret_->ClearState();

      // evaluate absolute L2 error
      (*abserror)[0] = sqrt((*errors)[0]);
      (*abserror)[1] = sqrt((*errors)[1]);
      (*abserror)[2] = sqrt((*errors)[2]);

      if (myrank_ == 0)
      {
        {
          std::cout.precision(8);
          std::cout << std::endl;
          std::cout << "---- Error norm for analytical solution -------------------" << std::endl;
          std::cout << "| absolute L_2 mixed variable error norm:   " << (*abserror)[0]
                    << std::endl;
          std::cout << "| absolute L_2 density        error norm:   " << (*abserror)[1]
                    << std::endl;
          std::cout << "| absolute L_2 momentum       error norm:   " << (*abserror)[2]
                    << std::endl;
          std::cout << "-----------------------------------------------------------" << std::endl;
          std::cout << std::endl;
        }

        // print last error in a seperate file

        // append error of the last time step to the error file
        if ((step_ == stepmax_) or (time_ == maxtime_))  // write results to file
        {
          std::ostringstream temp;
          const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
          const std::string fname = simulation + ".abserror";

          std::ofstream f;
          f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
          f << "#| " << simulation << "\n";
          f << "#| Step | Time | abs. L2-error mixed variable | abs. L2-error density | abs. "
               "L2-error momentum |\n";
          f << step_ << " " << time_ << " " << (*abserror)[0] << " " << (*abserror)[1] << " "
            << (*abserror)[2] << "\n";
          f.flush();
          f.close();
        }

        std::ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation + "_time.abserror";

        if (step_ == 1)
        {
          std::ofstream f;
          f.open(fname.c_str());
          f << "#| Step | Time | abs. L2-error mixed variable | abs. L2-error density | abs. "
               "L2-error momentum |\n";
          f << std::setprecision(10) << step_ << " " << std::setw(1) << std::setprecision(5)
            << time_ << std::setw(1) << std::setprecision(6) << " " << (*abserror)[0]
            << std::setw(1) << std::setprecision(6) << " " << (*abserror)[1] << std::setprecision(6)
            << " " << (*abserror)[2] << "\n";
          f.flush();
          f.close();
        }
        else
        {
          std::ofstream f;
          f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
          f << std::setprecision(10) << step_ << " " << std::setw(3) << std::setprecision(5)
            << time_ << std::setw(1) << std::setprecision(6) << " " << (*abserror)[0]
            << std::setw(1) << std::setprecision(6) << " " << (*abserror)[1] << std::setprecision(6)
            << " " << (*abserror)[2] << "\n";
          f.flush();
          f.close();
        }
      }
      return abserror;
    }
    break;
    default:
      dserror("Cannot calculate error. Unknown type of analytical test problem");
      break;
  }

  return Teuchos::null;
}



/*------------------------------------------------------------------------------------------------*
 |
 *------------------------------------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::Reset(bool completeReset, int numsteps, int iter)
{
  FluidImplicitTimeInt::Reset(completeReset, numsteps, iter);
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  intvelnp_ = LINALG::CreateVector(*intdofrowmap, true);
  intvelaf_ = LINALG::CreateVector(*intdofrowmap, true);
  intvelnm_ = LINALG::CreateVector(*intdofrowmap, true);
  intveln_ = LINALG::CreateVector(*intdofrowmap, true);
  intaccnp_ = LINALG::CreateVector(*intdofrowmap, true);
  intaccam_ = LINALG::CreateVector(*intdofrowmap, true);
  intaccnm_ = LINALG::CreateVector(*intdofrowmap, true);
  intaccn_ = LINALG::CreateVector(*intdofrowmap, true);
  if (discret_->Comm().MyPID() == 0)
    std::cout << "Number of degrees of freedom in HDG system: "
              << discret_->DofRowMap(0)->NumGlobalElements() << std::endl;
}



namespace
{
  // internal helper function for output
  void getNodeVectorsHDGWeakComp(DRT::Discretization& dis,
      const Teuchos::RCP<Epetra_Vector>& interiorValues,
      const Teuchos::RCP<Epetra_Vector>& traceValues, const int ndim,
      Teuchos::RCP<Epetra_MultiVector>& mixedvar, Teuchos::RCP<Epetra_Vector>& density,
      Teuchos::RCP<Epetra_MultiVector>& momentum, Teuchos::RCP<Epetra_Vector>& traceden,
      Teuchos::RCP<Epetra_MultiVector>& tracemom)
  {
    const int msd = ndim * (ndim + 1.0) / 2.0;

    // create dofsets for mixed variable, density and momentum at nodes
    if (density.get() == NULL || density->GlobalLength() != dis.NumGlobalNodes())
    {
      mixedvar.reset(new Epetra_MultiVector(*dis.NodeRowMap(), msd));
      density.reset(new Epetra_Vector(*dis.NodeRowMap()));
      momentum.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
    }
    traceden.reset(new Epetra_Vector(density->Map()));
    tracemom.reset(new Epetra_MultiVector(momentum->Map(), ndim));

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::interpolate_hdg_to_node);
    dis.SetState(1, "intvelnp", interiorValues);
    dis.SetState(0, "velnp", traceValues);
    std::vector<int> dummy;
    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(dis.NumMyRowNodes());
    mixedvar->PutScalar(0.);
    density->PutScalar(0.);
    momentum->PutScalar(0.);

    for (int el = 0; el < dis.NumMyColElements(); ++el)
    {
      DRT::Element* ele = dis.lColElement(el);
      if (interpolVec.M() == 0) interpolVec.Resize(ele->NumNode() * (msd + 1 + ndim + 1 + ndim));

      ele->Evaluate(params, dis, dummy, dummyMat, dummyMat, interpolVec, dummyVec, dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i = 0; i < ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = dis.NodeRowMap()->LID(node->Id());
        if (localIndex < 0) continue;
        touchCount[localIndex]++;
        for (int m = 0; m < msd; ++m)
          (*mixedvar)[m][localIndex] += interpolVec(i + m * ele->NumNode());
        (*density)[localIndex] += interpolVec(i + msd * ele->NumNode());
        for (int d = 0; d < ndim; ++d)
          (*momentum)[d][localIndex] += interpolVec(i + (msd + 1 + d) * ele->NumNode());
        (*traceden)[localIndex] += interpolVec(i + (msd + 1 + ndim) * ele->NumNode());
        for (int d = 0; d < ndim; ++d)
          (*tracemom)[d][localIndex] += interpolVec(i + (msd + 1 + ndim + 1 + d) * ele->NumNode());
      }
    }

    for (int i = 0; i < density->MyLength(); ++i)
    {
      for (int m = 0; m < msd; ++m) (*mixedvar)[m][i] /= touchCount[i];
      (*density)[i] /= touchCount[i];
      for (int d = 0; d < ndim; ++d) (*momentum)[d][i] /= touchCount[i];
      (*traceden)[i] /= touchCount[i];
      for (int d = 0; d < ndim; ++d) (*tracemom)[d][i] /= touchCount[i];
    }
    dis.ClearState();
  }
}  // namespace



/*----------------------------------------------------------------------*
 | output of solution vector to binio                      laspina 08/19|
 *----------------------------------------------------------------------*/
void FLD::TimIntHDGWeakComp::Output()
{
  // output of solution, currently only small subset of functionality
  if (step_ % upres_ == 0)
  {
    // step number and time
    output_->NewStep(step_, time_);

    // get number of spatial dimensions
    const unsigned int nsd = params_->get<int>("number of velocity degrees of freedom");

    // initialize trace variables
    Teuchos::RCP<Epetra_Vector> traceDen;
    Teuchos::RCP<Epetra_MultiVector> traceMom;

    // get node vectors
    getNodeVectorsHDGWeakComp(*discret_, intvelnp_, velnp_, nsd, interpolatedMixedVar_,
        interpolatedDensity_, interpolatedMomentum_, traceDen, traceMom);

    // get weakly compressible material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(
        INPAR::MAT::m_fluid_weakly_compressible);
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::WeaklyCompressibleFluid* actmat =
        static_cast<const MAT::PAR::WeaklyCompressibleFluid*>(mat);

    // evaluate derived variables
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity;
    Teuchos::RCP<Epetra_Vector> interpolatedPressure;
    interpolatedVelocity.reset(new Epetra_MultiVector(interpolatedMomentum_->Map(), nsd));
    interpolatedPressure.reset(new Epetra_Vector(interpolatedDensity_->Map()));
    for (int i = 0; i < interpolatedDensity_->MyLength(); ++i)
    {
      for (unsigned int d = 0; d < nsd; ++d)
        (*interpolatedVelocity)[d][i] = (*interpolatedMomentum_)[d][i] / (*interpolatedDensity_)[i];

      (*interpolatedPressure)[i] =
          actmat->refpressure_ +
          1.0 / actmat->comprcoeff_ * ((*interpolatedDensity_)[i] - actmat->refdensity_);
    }

    // write solution variables
    output_->WriteVector("Mixedvar", interpolatedMixedVar_, IO::nodevector);
    output_->WriteVector("Density", interpolatedDensity_, IO::nodevector);
    output_->WriteVector("Momentum", interpolatedMomentum_, IO::nodevector);
    output_->WriteVector("Trace_density", traceDen, IO::nodevector);
    output_->WriteVector("Trace_momentum", traceMom, IO::nodevector);

    // write derived variables
    output_->WriteVector("Velocity", interpolatedVelocity, IO::nodevector);
    output_->WriteVector("Pressure", interpolatedPressure, IO::nodevector);

    // write ALE variables
    if (alefluid_)
    {
      Teuchos::RCP<Epetra_MultiVector> AleDisplacement;
      AleDisplacement.reset(new Epetra_MultiVector(interpolatedMomentum_->Map(), nsd));
      for (int i = 0; i < interpolatedDensity_->MyLength(); ++i)
        for (unsigned int d = 0; d < nsd; ++d) (*AleDisplacement)[d][i] = (*dispnp_)[(i * nsd) + d];

      output_->WriteVector("Ale_displacement", AleDisplacement, IO::nodevector);
    }

    if (step_ == upres_ or step_ == 0) output_->WriteElementData(true);
  }
}
