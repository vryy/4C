/*----------------------------------------------------------------------*/
/*!
\file strutimint.cpp
\brief Time integration for spatially discretised structural dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint.H"

/*----------------------------------------------------------------------*/
/* Map stress input string to enum */
enum StruTimInt::StressEnum StruTimInt::MapStressStringToEnum
(
  const std::string name  //!< identifier
)
{
  if ( (name == "cauchy") or (name == "Cauchy") )
  {
    return stress_cauchy;
  }
  else if ( (name == "2pk") or (name == "2PK")
            or (name == "Yes") or (name == "yes") or (name == "YES") )
  {
    return stress_pk2;
  }
  else if ( (name == "No") or (name == "NO") or (name == "No") )
  {
    return stress_none;
  }
  else
  {
    dserror("Cannot handle (output) stress type %s", name.c_str());
    return stress_none;
  }
}

/*----------------------------------------------------------------------*/
/* Map strain input string to enum */
enum StruTimInt::StrainEnum StruTimInt::MapStrainStringToEnum
(
  const std::string name  //!< identifier
)
{
  if ( (name == "ea") or (name == "EA") )
  {
    return strain_ea;
  }
  else if ( (name == "gl") or (name == "GL")
            or (name == "Yes") or (name == "yes") or (name == "YES") )
  {
    return strain_gl;
  }
  else if ( (name == "No") or (name == "NO") or (name == "No") )
  {
    return strain_none;
  }
  else
  {
    dserror("Cannot handle (output) strain type %s", name.c_str());
    return strain_none;
  }
}

/*----------------------------------------------------------------------*/
/* constructor */
StruTimInt::StruTimInt
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  DRT::Discretization& actdis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
: discret_(actdis),
  myrank_(actdis.Comm().MyPID()),
  dofrowmap_(actdis.Filled() ? actdis.DofRowMap() : NULL),
  solver_(solver),
  solveradapttol_(Teuchos::getIntegralValue<int>(sdynparams,"ADAPTCONV")==1),
  solveradaptolbetter_(sdynparams.get<double>("ADAPTCONV_BETTER")),
  output_(output),
  printscreen_(true),  // ADD INPUT PARAMETER
  errfile_(xparams.get<FILE*>("err file")), 
  printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
  printiter_(true),  // ADD INPUT PARAMETER
  writerestartevery_(sdynparams.get<int>("RESTARTEVRY")),
  writestate_((bool) Teuchos::getIntegralValue<int>(ioparams,"STRUCT_DISP")),
  writestateevery_(sdynparams.get<int>("RESEVRYDISP")),
  writestrevery_(sdynparams.get<int>("RESEVRYSTRS")),
  writestress_(MapStressStringToEnum(ioparams.get<string>("STRUCT_STRESS"))),
  writestrain_(MapStrainStringToEnum(ioparams.get<string>("STRUCT_STRAIN"))),
  damping_((bool) Teuchos::getIntegralValue<int>(sdynparams,"DAMPING")),
  dampk_(sdynparams.get<double>("K_DAMP")),
  dampm_(sdynparams.get<double>("M_DAMP")),
  time_(0.0),  // HERE SHOULD BE SOMETHING LIKE (sdynparams.get<double>("TIMEINIT"))
  timen_(0.0),
  dt_(sdynparams.get<double>("TIMESTEP")),
  timemax_(sdynparams.get<double>("MAXTIME")),
  stepmax_(sdynparams.get<int>("NUMSTEP")),
  step_(0),
  stepn_(0),
  dirichtoggle_(Teuchos::null),
  invtoggle_(Teuchos::null),
  zeros_(Teuchos::null),
  state_(Teuchos::null),
  staten_(Teuchos::null),
  stiff_(Teuchos::null),
  mass_(Teuchos::null),
  damp_(Teuchos::null)
{
  // welcome user
  if (myrank_ == 0)
  {
    std::cout << "Welcome to Structural Time Integration " << std::endl;
    std::cout << "     __o__                          __o__       " << std::endl;
    std::cout << "__  /-----\\__                  __  /-----\\__           " << std::endl;
    std::cout << "\\ \\/       \\ \\    |       \\    \\ \\/       \\ \\          " << std::endl;
    std::cout << " \\ |  tea  | |    |-------->    \\ |  tea  | |          " << std::endl;
    std::cout << "  \\|       |_/    |       /      \\|       |_/          " << std::endl;
    std::cout << "    \\_____/   ._                   \\_____/   ._ _|_ /| " << std::endl;
    std::cout << "              | |                            | | |   | " << std::endl;
    std::cout << std::endl;
  }

  // check wether discretisation has been completed
  if (not discret_.Filled())
  {
    dserror("Discretisation is not complete!");
  }

  // time state 
  timen_ = time_;  // set target time to initial time
  step_ = 0;  // time step

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //if (not discret_.Filled())
  //{
  //  discret_.FillComplete();
  //}
  //dofrowmap_ = discret_.DofRowMap();

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap_, true);

  // Dirichlet vector
  // vector of full length; for each component
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap_, true);
  // set Dirichlet toggle vector
  {
    Teuchos::ParameterList p;
    p.set("total time", timen_);
    discret_.EvaluateDirichlet(p, zeros_, null, null, dirichtoggle_);
    zeros_->PutScalar(0.0); // just in case of change
  }
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  invtoggle_ = LINALG::CreateVector(*dofrowmap_, false);
  // compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0, *dirichtoggle_, 1.0);

  // create 1 state vector for last step at t_n
  state_ = Teuchos::rcp(new StruTimIntState(0,0,dofrowmap_,true));
  // associate convenience pointers
  dis_ = state_->Dis(0);
  vel_ = state_->Vel(0);
  acc_ = state_->Acc(0);

  // create 1 state vector for last step at t_{n+1}
  staten_ = Teuchos::rcp(new StruTimIntState(1,1,dofrowmap_,true));
  // associate convenience pointers
  disn_ = staten_->Dis(1);
  veln_ = staten_->Vel(1);
  accn_ = staten_->Acc(1);

  // create empty matrices
  stiff_ = Teuchos::rcp(
    new LINALG::SparseMatrix(*dofrowmap_, 81, true, false)
  );
  mass_ = Teuchos::rcp(
    new LINALG::SparseMatrix(*dofrowmap_, 81, true, false)
  );
  if (damping_)
  {
    damp_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*dofrowmap_, 81, true, false)
    );
  }

  // determine mass, damping and initial accelerations
  DetermineMassDampConsistAccel();

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* equilibrate system at initial state
 * and identify consistent accelerations */
void StruTimInt::DetermineMassDampConsistAccel()
{
  // temporary force vectors in this routine
  Teuchos::RCP<Epetra_Vector> fext 
    = LINALG::CreateVector(*dofrowmap_, true); // external force
  Teuchos::RCP<Epetra_Vector> fint 
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  // overwrite initial state vectors with DirichletBCs
  ApplyDirichletBC(time_, dis_, vel_, acc_);

  // get external force
  ApplyForceExternal(time_, dis_, fext);
  
  // initialise matrices
  stiff_->Zero();
  mass_->Zero();

  // get initial internal force and stiffness and mass
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "calc_struct_nlnstiffmass");
    // other parameters that might be needed by the elements
    p.set("total time", time_);
    p.set("delta time", dt_);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement", zeros_);
    discret_.SetState("displacement", dis_);
    //discret_.SetState("velocity",vel_); // not used at the moment
    discret_.Evaluate(p, stiff_, mass_, fint, null, null);
    discret_.ClearState();
  }

  // finish mass matrix
  mass_->Complete();

  // close stiffness matrix
  stiff_->Complete();

  // build Rayleigh damping matrix if desired
  if (damping_)
  {
    damp_->Add(*stiff_, false, dampk_, 0.0);
    damp_->Add(*mass_, false, dampm_, 1.0);
    damp_->Complete();
  }

  // calculate consistent initial accelerations
  // WE MISS:
  //   - surface stress forces
  //   - potential forces
  {
    Teuchos::RCP<Epetra_Vector> rhs 
      = LINALG::CreateVector(*dofrowmap_, true);
    if (damping_)
    {
      damp_->Multiply(false, *vel_, *rhs);
    }
    rhs->Update(-1.0, *fint, 1.0, *fext, -1.0);
    Epetra_Vector rhscopy = Epetra_Vector(*rhs);
    rhs->Multiply(1.0, *invtoggle_, rhscopy, 0.0);
    solver_.Solve(mass_->EpetraMatrix(), acc_, rhs, true, true);
  }

  // leave this
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate Dirichlet BC at t_{n+1} */
void StruTimInt::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> acc
)
{
  // apply DBCs
  // needed parameters
  ParameterList p;
  p.set("total time", time);  // target time
  
  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_.ClearState();
  discret_.EvaluateDirichlet(p, dis, vel, acc, dirichtoggle_);
  discret_.ClearState();

  // compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0, *dirichtoggle_, 1.0);

  // ciao
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate external forces at t_{n+1} */
void StruTimInt::ApplyForceExternal
(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
  Teuchos::RCP<Epetra_Vector>& fext  //!< external force
)
{
  ParameterList p;
  // action for elements
  p.set("action", "calc_struct_eleload");
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_.ClearState();
  discret_.SetState("displacement", dis);
  // get load vector
  discret_.EvaluateNeumann(p, *fext);
  discret_.ClearState();

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force, its stiffness at state */
void StruTimInt::ApplyForceStiffInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // residual displacements
  Teuchos::RCP<Epetra_Vector> fint,  // internal force
  Teuchos::RCP<LINALG::SparseMatrix> stiff  // stiffness matrix
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  const std::string action = "calc_struct_nlnstiff";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  // set vector values needed by elements
  discret_.ClearState();
  discret_.SetState("residual displacement", disi);
  discret_.SetState("displacement", dis);
  //discret_.SetState("velocity", veln_); // not used at the moment
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_.Evaluate(p, stiff, null, fint, null, null);
  discret_.ClearState();
  
  // that's it
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force */
void StruTimInt::ApplyForceInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // incremental displacements
  Teuchos::RCP<Epetra_Vector> fint  // internal force
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  const std::string action = "calc_struct_internalforce";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  // set vector values needed by elements
  discret_.ClearState();
  discret_.SetState("residual displacement", disi);  // these are incremental
  discret_.SetState("displacement", dis);
  //discret_.SetState("velocity", veln_); // not used at the moment
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_.Evaluate(p, null, null, fint, null, null);
  discret_.ClearState();
  
  // where the fun starts
  return;
}

/*----------------------------------------------------------------------*/
/* integrate */
void StruTimInt::Integrate()
{
  // set target time and step
  timen_ = time_ + dt_;
  stepn_ = step_ + 1;

  // time loop
  while ( (timen_ <= timemax_) and (stepn_ <= stepmax_) )
  {
    // integrate time step
    // after this step we hold disn_, etc
    IntegrateStep();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    UpdateStep();

    // update time and step
    time_ = timen_;
    step_ = stepn_;
    // 
    timen_ += dt_;
    stepn_ += 1;

    // print info about finished time step
    PrintStep();

    // write output
    Output();
  }

  // that's it
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
