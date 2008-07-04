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
    staten_(Teuchos::null)
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

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* integrate */
void StruTimIntImpl::Integrate()
{
  // set target time and step
  timen_ = time_ + dt_;
  stepn_ = step_ + 1;

  // time loop
  while ( (timen_ <= timemax_) and (step_ <= stepmax_) )
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

#endif  // #ifdef CCADISCRET
