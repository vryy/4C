/*----------------------------------------------------------------------*/
/*!
\file strtimint_expl.cpp
\brief Explicit time integration for spatial discretised 
       structural dynamics

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
#include <sstream>

#include "strtimint.H"
#include "strtimint_expl.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntExpl:: TimIntExpl
(
  const Teuchos::ParameterList& ioparams,  //!< ioflags
  const Teuchos::ParameterList& sdynparams,  //!< input parameters
  const Teuchos::ParameterList& xparams,  //!< extra flags
  Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
  Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
  Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
)
: TimInt
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    output
  )
{
  // explicit time integrators cannot handle constraints
  if (conman_->HaveConstraint())
  {
    dserror("Explicit TIS cannot handle constraints");
  }

  // get away
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntExpl::PrintStep()
{
  // print out
  if ( (myrank_ == 0) and printscreen_ )
  {
    PrintStepText(stdout);
  }

  if (printerrfile_)
  {
    PrintStepText(errfile_);
  }

  // fall asleep
  return;
}

/*----------------------------------------------------------------------*/
/* print step summary */
void STR::TimIntExpl::PrintStepText
(
  FILE* ofile
)
{
  fprintf(ofile,
          "Finalised: step %6d"
          " | nstep %6d"
          " | time %-14.8E"
          " | dt %-14.8E"
          " | numiter %3d\n",
          step_, stepmax_, (*time_)[0], (*dt_)[0], 0);
  // print a beautiful line made exactly of 80 dashes
  fprintf(ofile,
          "--------------------------------------------------------------"
          "------------------\n");
  // do it, print now!
  fflush(ofile);

  // fall asleep
  return;
}

/*----------------------------------------------------------------------*/
/* output to file
 * originally by mwgee 03/07 */
void STR::TimIntExpl::OutputStep()
{
  // this flag is passed
  bool datawritten = false;

  // output restart (try this first)
  OutputRestart(datawritten);

  // output results (not necessary if restart in same step)
  OutputState(datawritten);

  // output stress & strain
  OutputStressStrain(datawritten);

  // what's next?
  return;
}

/*----------------------------------------------------------------------*/
/* write restart
 * originally mwgee 03/07 */
void STR::TimIntExpl::OutputRestart
(
  bool& datawritten
)
{
  // write restart step
  if (writerestartevery_ and (step_%writerestartevery_ == 0) )
  {
    // Yes, we are going to write...
    datawritten = true;

    // write restart output, please
    output_->WriteMesh(step_, (*time_)[0]);
    output_->NewStep(step_, (*time_)[0]);
    output_->WriteVector("displacement", (*dis_)(0));
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
    //output_->WriteVector("fexternal", fext_);  // CURRENTLY NOT AVAILABLE THINK OF SCENARIO

/*
    // surface stress
    if (surfstressman_ != Teuchos::null)
    {
      RCP<Epetra_Map> surfrowmap = surfstressman_->GetSurfRowmap();
      RCP<Epetra_Vector> A = rcp(new Epetra_Vector(*surfrowmap, true));
      RCP<Epetra_Vector> con = rcp(new Epetra_Vector(*surfrowmap, true));
      surfstressman_->GetHistory(A,con);
      output_->WriteVector("Aold", A);
      output_->WriteVector("conquot", con);
    }
    
    // potential forces
    if (potman_ != Teuchos::null)
    {
      RCP<Epetra_Map> surfrowmap = potman_->GetSurfRowmap();
      RCP<Epetra_Vector> A = rcp(new Epetra_Vector(*surfrowmap, true));
      potman_->GetHistory(A);
      output_->WriteVector("Aold", A);
    }

    // constraints
    if (conman_->HaveConstraint())
    {
      output_->WriteDouble("uzawaparameter",
                          uzawasolv_->GetUzawaParameter());
    }

    // info dedicated to user's eyes staring at standard out
    if ( (myrank_ == 0) and printscreen_)
    { 
      printf("====== Restart written in step %d\n", step_);
      fflush(stdout);
    }
*/

    // info dedicated to processor error file
    if (printerrfile_)
    {
      fprintf(errfile_, "====== Restart written in step %d\n", step_);
      fflush(errfile_);
    }
  }

  // we will say what we did
  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities and accelerations
 * originally mwgee 03/07 */
void STR::TimIntExpl::OutputState
(
  bool& datawritten
)
{
  if ( writestate_ 
       and writestateevery_ and (step_%writestateevery_ == 0)
       and (not datawritten) )
  {
    // Yes, we are going to write...
    datawritten = true;

    // write now
    output_->NewStep(step_, (*time_)[0]);
    output_->WriteVector("displacement", (*dis_)(0));
    output_->WriteVector("velocity", (*vel_)(0));
    output_->WriteVector("acceleration", (*acc_)(0));
    //output_->WriteVector("fexternal",fext_);  // CURRENTLY NOT AVAILABLE
    output_->WriteElementData();
  }

  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* stress output
 * originally by lw */
void STR::TimIntExpl::OutputStressStrain
(
  bool& datawritten
)
{
  // do stress calculation and output
  if ( writestrevery_
       and ( (writestress_ != stress_none)
             or (writestrain_ != strain_none) )
       and (step_%writestrevery_ == 0) )
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
    discret_->Evaluate(p, null, null, null, null, null);
    discret_->ClearState();

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
  }

  // leave me alone
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
