/*----------------------------------------------------------------------*/
/*!
\file strtimint.cpp
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
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "strtimint_mstep.H"
#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* Map stress input string to enum */
enum STR::TimInt::StressEnum STR::TimInt::MapStressStringToEnum
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
enum STR::TimInt::StrainEnum STR::TimInt::MapStrainStringToEnum
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
/* provide string for identifying enum */
std::string STR::TimInt::MapNameEnumToString
(
  const enum NameEnum term  //!< the enum
)
{
  switch (term)
  {
  case name_statics :
    return "Statics";
    break;
  case name_genalpha :
    return "GenAlpha";
    break;
  case name_onesteptheta :
    return "OneStepTheta";
    break;
  case name_gemm :
    return "GEMM";
    break;
  case name_ab2 :
    return "AdamsBashforth2";
    break;
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }
}

/*----------------------------------------------------------------------*/
/* print tea time logo */
void STR::TimInt::Logo()
{
  std::cout << "Welcome to Structural Time Integration " << std::endl;
  std::cout << "     __o__                          __o__" << std::endl;
  std::cout << "__  /-----\\__                  __  /-----\\__" << std::endl;
  std::cout << "\\ \\/       \\ \\    |       \\    \\ \\/       \\ \\" << std::endl;
  std::cout << " \\ |  tea  | |    |-------->    \\ |  tea  | |" << std::endl;
  std::cout << "  \\|       |_/    |       /      \\|       |_/" << std::endl;
  std::cout << "    \\_____/   ._                   \\_____/   ._ _|_ /|" << std::endl;
  std::cout << "              | |                            | | |   |" << std::endl;
  std::cout << std::endl;
}

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimInt::TimInt
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: discret_(actdis),
  myrank_(actdis->Comm().MyPID()),
  dofrowmap_(actdis->Filled() ? actdis->DofRowMap() : NULL),
  solver_(solver),
  solveradapttol_(Teuchos::getIntegralValue<int>(sdynparams,"ADAPTCONV")==1),
  solveradaptolbetter_(sdynparams.get<double>("ADAPTCONV_BETTER")),
  dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor())),
  output_(output),
  printlogo_(true),  // DON'T EVEN DARE TO SET THIS TO FALSE
  printscreen_(true),  // ADD INPUT PARAMETER
  errfile_(xparams.get<FILE*>("err file")),
  printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
  printiter_(true),  // ADD INPUT PARAMETER
  writerestartevery_(sdynparams.get<int>("RESTARTEVRY")),
  writestate_((bool) Teuchos::getIntegralValue<int>(ioparams,"STRUCT_DISP")),
  writestateevery_(sdynparams.get<int>("RESEVRYDISP")),
  writestrevery_(sdynparams.get<int>("RESEVRYSTRS")),
  writestress_(MapStressStringToEnum(ioparams.get<std::string>("STRUCT_STRESS"))),
  writestrain_(MapStrainStringToEnum(ioparams.get<std::string>("STRUCT_STRAIN"))),
  writeenergyevery_(sdynparams.get<int>("RESEVRYERGY")),
  energyfile_(NULL),
  damping_(Teuchos::getIntegralValue<INPAR::STR::DampKind>(sdynparams,"DAMPING")),
  dampk_(sdynparams.get<double>("K_DAMP")),
  dampm_(sdynparams.get<double>("M_DAMP")),
  conman_(Teuchos::null),
  consolv_(Teuchos::null),
  surfstressman_(Teuchos::null),
  potman_(Teuchos::null),
  time_(Teuchos::null),
  timen_(0.0),
  dt_(Teuchos::null),
  timemax_(sdynparams.get<double>("MAXTIME")),
  stepmax_(sdynparams.get<int>("NUMSTEP")),
  step_(Teuchos::null),
  stepn_(0),
  zeros_(Teuchos::null),
  dis_(Teuchos::null),
  vel_(Teuchos::null),
  acc_(Teuchos::null),
  disn_(Teuchos::null),
  veln_(Teuchos::null),
  accn_(Teuchos::null),
  stiff_(Teuchos::null),
  mass_(Teuchos::null),
  damp_(Teuchos::null)
{
  // welcome user
  if ( (printlogo_) and (myrank_ == 0) )
  {
    Logo();
  }

  // check wether discretisation has been completed
  if (not discret_->Filled())
  {
    dserror("Discretisation is not complete!");
  }


  // time state
  time_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, 0.0));  // HERE SHOULD BE SOMETHING LIKE (sdynparams.get<double>("TIMEINIT")) -- you may not believe it, but time does not always start at zero
  dt_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, sdynparams.get<double>("TIMESTEP")));
  step_ = 0;
  timen_ = (*time_)[0] + (*dt_)[0];  // set target time to initial time plus step size
  stepn_ = step_ + 1;

  // output file for energy
  if ( (writeenergyevery_ != 0) and (myrank_ == 0) )
    AttachEnergyFile();

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap_, true);

  // Map containing Dirichlet DOFs
  {
    Teuchos::ParameterList p;
    p.set("total time", timen_);
    discret_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null, 
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }
  
  // displacements D_{n}
  dis_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new TimIntMStep<Epetra_Vector>(0, 0, dofrowmap_, true));

  // displacements D_{n+1} at t_{n+1}
  disn_ = LINALG::CreateVector(*dofrowmap_, true);
  // velocities V_{n+1} at t_{n+1}
  veln_ = LINALG::CreateVector(*dofrowmap_, true);
  // accelerations A_{n+1} at t_{n+1}
  accn_ = LINALG::CreateVector(*dofrowmap_, true);

  // create empty matrices
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));
  mass_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));
  if (damping_ == INPAR::STR::damp_rayleigh)
    damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81, true, true));

  // initialize constraint manager
  conman_ = Teuchos::rcp(new UTILS::ConstrManager(discret_,
                                                  (*dis_)(0),
                                                  sdynparams));
  // initialize constraint solver iff constraints are defined
  if (conman_->HaveConstraint())
  {
    consolv_ = Teuchos::rcp(new UTILS::ConstraintSolver(discret_,
                                                        *solver_,
                                                        dbcmaps_,
                                                        sdynparams)); 
  }
  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  dofrowmap_ = discret_->DofRowMap();

  // Check for surface stress conditions due to interfacial phenomena
  {
    vector<DRT::Condition*> surfstresscond(0);
    discret_->GetCondition("SurfaceStress", surfstresscond);
    if (surfstresscond.size())
    {
      surfstressman_ = rcp(new UTILS::SurfStressManager(*discret_));
    }
  }

  // Check for potential conditions
  {
    vector<DRT::Condition*> potentialcond(0);
    discret_->GetCondition("Potential", potentialcond);
    if (potentialcond.size())
    {
      potman_ = rcp(new UTILS::PotentialManager(Discretization(), *discret_));
    }
  }

  // stay with us
  return;
}

/*----------------------------------------------------------------------*/
/* equilibrate system at initial state
 * and identify consistent accelerations */
void STR::TimInt::DetermineMassDampConsistAccel()
{
  // temporary force vectors in this routine
  Teuchos::RCP<Epetra_Vector> fext
    = LINALG::CreateVector(*dofrowmap_, true); // external force
  Teuchos::RCP<Epetra_Vector> fint
    = LINALG::CreateVector(*dofrowmap_, true); // internal force

  // overwrite initial state vectors with DirichletBCs
  ApplyDirichletBC((*time_)[0], (*dis_)(0), (*vel_)(0), (*acc_)(0), false);

  // get external force
  ApplyForceExternal((*time_)[0], (*dis_)(0), (*vel_)(0), fext);

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
    p.set("total time", (*time_)[0]);
    p.set("delta time", (*dt_)[0]);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement", zeros_);
    discret_->SetState("displacement", (*dis_)(0));
    if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", (*vel_)(0));
    discret_->Evaluate(p, stiff_, mass_, fint, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // finish mass matrix
  mass_->Complete();

  // close stiffness matrix
  stiff_->Complete();

  // build Rayleigh damping matrix if desired
  if (damping_ == INPAR::STR::damp_rayleigh)
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
    if (damping_ == INPAR::STR::damp_rayleigh)
    {
      damp_->Multiply(false, (*vel_)[0], *rhs);
    }
    rhs->Update(-1.0, *fint, 1.0, *fext, -1.0);
    // blank RHS on DBC DOFs
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), rhs);
    solver_->Solve(mass_->EpetraMatrix(), (*acc_)(0), rhs, true, true);
  }

  // We need to reset the stiffness matrix because its graph (topology)
  // is not finished yet in case of constraints and posssibly other side
  // effects (basically managers).
  stiff_->Reset();

  // leave this hell
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate Dirichlet BC at t_{n+1} */
void STR::TimInt::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> acc,
  bool recreatemap
)
{
  // apply DBCs
  // needed parameters
  ParameterList p;
  p.set("total time", time);  // target time

  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_->ClearState();
  if (recreatemap)
  {
    discret_->EvaluateDirichlet(p, dis, vel, acc,
                                Teuchos::null, dbcmaps_);
  }
  else
  {
    discret_->EvaluateDirichlet(p, dis, vel, acc,
                                Teuchos::null, Teuchos::null);
  }
  discret_->ClearState();

  // ciao
  return;
}

/*----------------------------------------------------------------------*/
/* Update time and step counter */
void STR::TimInt::UpdateStepTime()
{
  // update time and step
  time_->UpdateSteps(timen_);  // t_{n} := t_{n+1}, etc
  step_ = stepn_;  // n := n+1
  // 
  timen_ += (*dt_)[0];
  stepn_ += 1;

  // new deal
  return;
}

/*----------------------------------------------------------------------*/
/* Reset configuration after time step */
void STR::TimInt::ResetStep()
{
  // reset state vectors
  disn_->Update(1.0, (*dis_)[0], 0.0);
  veln_->Update(1.0, (*vel_)[0], 0.0);
  accn_->Update(1.0, (*acc_)[0], 0.0);

  // reset anything that needs to be reset at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("action", "calc_struct_reset_istep");
    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  // I am gone
  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart values */
void STR::TimInt::ReadRestart
(
  const int step
)
{
  IO::DiscretizationReader reader(discret_, step);
  step_ = reader.ReadInt("step");
  if (step_ != step) dserror("Time step on file not equal to given step");
  time_ = Teuchos::rcp(new TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
  ReadRestartState();
  ReadRestartForce();
  ReadRestartConstraint();
  ReadRestartSurfstress();
  ReadRestartMultiScale();
  // fix pointer to #dofrowmap_, which has not really changed, but is
  // located at different place
  dofrowmap_ = discret_->DofRowMap();
}

/*----------------------------------------------------------------------*/
/* Read and set restart state */
void STR::TimInt::ReadRestartState()
{
  IO::DiscretizationReader reader(discret_, step_);
  reader.ReadVector(disn_, "displacement");
  dis_->UpdateSteps(*disn_);
  reader.ReadVector(veln_, "velocity");
  vel_->UpdateSteps(*veln_);
  reader.ReadVector(accn_, "acceleration");
  acc_->UpdateSteps(*accn_);
  reader.ReadMesh(step_);
  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart values for constraints */
void STR::TimInt::ReadRestartConstraint()
{
  if (conman_->HaveConstraint())
  {
    IO::DiscretizationReader reader(discret_, step_);
    double uzawatemp = reader.ReadDouble("uzawaparameter");
    consolv_->SetUzawaParameter(uzawatemp);
    Teuchos::RCP<Epetra_Map> constrmap=conman_->GetConstraintMap();
    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*constrmap, true);
    reader.ReadVector(tempvec, "lagrmultiplier");
    conman_->SetLagrMultVector(tempvec);
    reader.ReadVector(tempvec, "refconval");
    conman_->SetRefBaseValues(tempvec, (*time_)[0]);
  }
}


/*----------------------------------------------------------------------*/
/* Read and set restart values for constraints */
void STR::TimInt::ReadRestartSurfstress()
{
  if (surfstressman_ != Teuchos::null)
  {
    IO::DiscretizationReader reader(discret_, step_);
    Teuchos::RCP<Epetra_Map> surfmap = surfstressman_->GetSurfRowmap();
    RCP<Epetra_Vector> A = LINALG::CreateVector(*surfmap,true);
    RCP<Epetra_Vector> con = LINALG::CreateVector(*surfmap,true);
    RCP<Epetra_Vector> gamma = LINALG::CreateVector(*surfmap,true);
    reader.ReadVector(A, "A");
    reader.ReadVector(con, "con");
    reader.ReadVector(gamma, "gamma");
    surfstressman_->SetHistory(A, con, gamma);
  }
}

/*----------------------------------------------------------------------*/
/* Read and set restart values for multi-scale */
void STR::TimInt::ReadRestartMultiScale()
{
  if (DRT::Problem::Instance()->ProblemType() == "struct_multi")
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action", "multi_readrestart");
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }
}

/*----------------------------------------------------------------------*/
/* output to file
 * originally by mwgee 03/07 */
void STR::TimInt::OutputStep()
{
  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if (writerestartevery_ and (step_%writerestartevery_ == 0) )
  {
    OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if ( writestate_
       and writestateevery_ and (step_%writestateevery_ == 0)
       and (not datawritten) )
  {
    OutputState(datawritten);
  }

  // output stress & strain
  if ( writestrevery_
       and ( (writestress_ != stress_none)
             or (writestrain_ != strain_none) )
       and (step_%writestrevery_ == 0) )
  {
    OutputStressStrain(datawritten);
  }

  // output energy
  if ( writeenergyevery_ and (step_%writeenergyevery_ == 0) )
  {
    OutputEnergy();
  }

  // what's next?
  return;
}

/*----------------------------------------------------------------------*/
/* write restart
 * originally by mwgee 03/07 */
void STR::TimInt::OutputRestart
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // write restart output, please
  output_->WriteMesh(step_, (*time_)[0]);
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("velocity", (*vel_)(0));
  output_->WriteVector("acceleration", (*acc_)(0));
  output_->WriteVector("fexternal", Fext());

  // surface stress
  if (surfstressman_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Map> surfrowmap
      = surfstressman_->GetSurfRowmap();
    Teuchos::RCP<Epetra_Vector> A
      = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
    Teuchos::RCP<Epetra_Vector> con
      = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
    Teuchos::RCP<Epetra_Vector> gamma
      = Teuchos::rcp(new Epetra_Vector(*surfrowmap, true));
    surfstressman_->GetHistory(A,con,gamma);
    output_->WriteVector("A", A);
    output_->WriteVector("con", con);
    output_->WriteVector("gamma", gamma);
  }

  // constraints
  if (conman_->HaveConstraint())
  {
    output_->WriteDouble("uzawaparameter",
                          consolv_->GetUzawaParameter());
    output_->WriteVector("lagrmultiplier",
                          conman_->GetLagrMultVector());
    output_->WriteVector("refconval",
                          conman_->GetRefBaseValues());
  }

  // info dedicated to user's eyes staring at standard out
  if ( (myrank_ == 0) and printscreen_)
  {
    printf("====== Restart written in step %d\n", step_);
    fflush(stdout);
  }

  // info dedicated to processor error file
  if (printerrfile_)
  {
    fprintf(errfile_, "====== Restart written in step %d\n", step_);
    fflush(errfile_);
  }

  // we will say what we did
  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities and accelerations
 * originally by mwgee 03/07 */
void STR::TimInt::OutputState
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // write now
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("velocity", (*vel_)(0));
  output_->WriteVector("acceleration", (*acc_)(0));
  output_->WriteVector("fexternal", Fext());
  output_->WriteElementData();

  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* stress calculation and output
 * originally by lw */
void STR::TimInt::OutputStressStrain
(
  bool& datawritten
)
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  p.set("action", "calc_struct_stress");
  // other parameters that might be needed by the elements
  p.set("total time", (*time_)[0]);
  p.set("delta time", (*dt_)[0]);

  // stress
  if (writestress_ == stress_cauchy)
  {
    // output of Cauchy stresses instead of 2PK stresses
    p.set("cauchy", true);
  }
  else
  {
    // this will produce 2nd PK stress ????
    p.set("cauchy", false);
  }
  Teuchos::RCP<std::vector<char> > stressdata
    = Teuchos::rcp(new std::vector<char>());
  p.set("stress", stressdata);

  // strain
  if (writestrain_ == strain_ea)
  {
    p.set("iostrain", "euler_almansi");
  }
  else if (writestrain_ == strain_gl)
  {
    // WILL THIS CAUSE TROUBLE ????
    // THIS STRING DOES NOT EXIST IN SO3
    p.set("iostrain", "green_lagrange");
  }
  else
  {
    p.set("iostrain", "none");
  }
  Teuchos::RCP<std::vector<char> > straindata
    = Teuchos::rcp(new std::vector<char>());
  p.set("strain", straindata);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("residual displacement", zeros_);
  discret_->SetState("displacement", (*dis_)(0));
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // Make new step
  if (not datawritten)
  {
    output_->NewStep(step_, (*time_)[0]);
  }
  datawritten = true;

  // write stress
  if (writestress_ != stress_none)
  {
    std::string stresstext = "";
    if (writestress_ == stress_cauchy)
    {
      stresstext = "gauss_cauchy_stresses_xyz";
    }
    else if (writestress_ == stress_pk2)
    {
      stresstext = "gauss_2PK_stresses_xyz";
    }
    output_->WriteVector(stresstext, *stressdata,
                         *(discret_->ElementColMap()));
  }

  // write strain
  if (writestrain_ != strain_none)
  {
    std::string straintext = "";
    if (writestrain_ == strain_ea)
    {
      straintext = "gauss_EA_strains_xyz";
    }
    else
    {
      straintext = "gauss_GL_strains_xyz";
    }
    output_->WriteVector(straintext, *straindata,
                         *(discret_->ElementColMap()));
  }

  // leave me alone
  return;
}

/*----------------------------------------------------------------------*/
/* output system energies */
void STR::TimInt::OutputEnergy()
{
  // internal/strain energy
  double intergy = 0.0;  // total internal energy
  {
    ParameterList p;
    // other parameters needed by the elements
    p.set("action", "calc_struct_energy");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("displacement", (*dis_)(0));
    // get energies
    Teuchos::RCP<Epetra_SerialDenseVector> energies
      = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(p, energies);
    discret_->ClearState();
    intergy = (*energies)(0);
  }

  // global calculation of kinetic energy
  double kinergy = 0.0;  // total kinetic energy
  {
    Teuchos::RCP<Epetra_Vector> linmom
      = LINALG::CreateVector(*dofrowmap_, true);
    mass_->Multiply(false, (*vel_)[0], *linmom);
    linmom->Dot((*vel_)[0], &kinergy);
    kinergy *= 0.5;
  }

  // external energy
  double extergy = 0.0;  // total external energy
  {
    // WARNING: This will only work with dead loads!!!
    Teuchos::RCP<Epetra_Vector> fext = Fext();
    fext->Dot((*dis_)[0], &extergy);
  }

  // total energy
  double totergy = kinergy + intergy - extergy;

  // the output
  if (myrank_ == 0)
  {
    *energyfile_ << " " << std::setw(9) << step_
                 << std::scientific  << std::setprecision(16)
                 << " " << (*time_)[0]
                 << " " << totergy
                 << " " << kinergy
                 << " " << intergy
                 << " " << extergy
                 << std::endl;
  }

  // in god we trust
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate external forces at t_{n+1} */
void STR::TimInt::ApplyForceExternal
(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
  const Teuchos::RCP<Epetra_Vector> vel,  //!< velocity state
  Teuchos::RCP<Epetra_Vector>& fext  //!< external force
)
{
  ParameterList p;
  // other parameters needed by the elements
  p.set("total time", time);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("displacement", dis);
  if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", vel);
  // get load vector
  discret_->EvaluateNeumann(p, *fext);
  discret_->ClearState();

  // go away
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force, its stiffness at state */
void STR::TimInt::ApplyForceStiffInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // residual displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
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
  discret_->ClearState();
  discret_->SetState("residual displacement", disi);
  discret_->SetState("displacement", dis);
  if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // that's it
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force */
void STR::TimInt::ApplyForceInternal
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state
  const Teuchos::RCP<Epetra_Vector> disi,  // incremental displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
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
  discret_->ClearState();
  discret_->SetState("residual displacement", disi);  // these are incremental
  discret_->SetState("displacement", dis);
  if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     fint, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // where the fun starts
  return;
}

/*----------------------------------------------------------------------*/
/* integrate */
void STR::TimInt::Integrate()
{
  // set target time and step
  timen_ = (*time_)[0] + (*dt_)[0];
  stepn_ = step_ + 1;

  // time loop
  while ( (timen_ <= timemax_) and (stepn_ <= stepmax_) )
  {
    // integrate time step
    // after this step we hold disn_, etc
    IntegrateStep();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    UpdateStepState();

    // update time and step
    UpdateStepTime();

    // print info about finished time step
    PrintStep();

    // write output
    OutputStep();
  }

  // that's it
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
