/*----------------------------------------------------------------------*/
/*! \file

\brief Structural time integration with one-step-theta for cell migration

\level 1

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/

#include "str_timint_ost_immersed.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_sparsematrix.H"

/*======================================================================*/
/* constructor */
STR::TimIntOneStepThetaImmersed::TimIntOneStepThetaImmersed(
    const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<LINALG::Solver> contactsolver, Teuchos::RCP<IO::DiscretizationWriter> output)
    : TimIntOneStepTheta(
          timeparams, ioparams, sdynparams, xparams, actdis, solver, contactsolver, output)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
  return;
}


/*----------------------------------------------------------------------------------------------*
 * Initialize this class                                                            rauch 05/17 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntOneStepThetaImmersed::Init(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver)
{
  // call Init() in base class
  STR::TimIntOneStepTheta::Init(timeparams, sdynparams, xparams, actdis, solver);

  return;
}


/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 05/17 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntOneStepThetaImmersed::Setup()
{
  // call Setup() in base class
  STR::TimIntOneStepTheta::Setup();

  return;
}


/*-------------------------------------------------------------------------------------------*
 * Create matrices when setting up time integrator
 *-------------------------------------------------------------------------------------------*/
void STR::TimIntOneStepThetaImmersed::createFields()
{
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMapView(), true);

  // Map containing Dirichlet DOFs
  {
    Teuchos::ParameterList p;
    p.set("total time", timen_);
    discret_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);  // just in case of change
  }

  // create empty matrices
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, false, true));
  mass_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, false, true));
  if (damping_ != INPAR::STR::damp_none)
    damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, false, true));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::TimIntOneStepThetaImmersed::DetermineMass()
{
  Teuchos::RCP<Epetra_Vector> acc_aux = LINALG::CreateVector(*DofRowMapView(), true);
  Teuchos::RCP<LINALG::SparseOperator> stiff_aux =
      Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, false, true));

  // initialise matrices
  mass_->Zero();

  // compute new inner radius
  discret_->ClearState();
  discret_->SetState(0, "displacement", (*dis_)(0));

  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // action for elements
  if (lumpmass_ == false) p.set("action", "calc_struct_nlnstiffmass");
  // lumping the mass matrix
  else
    p.set("action", "calc_struct_nlnstifflmass");
  // other parameters that might be needed by the elements
  p.set("total time", (*time_)[0]);
  p.set("delta time", (*dt_)[0]);

  // set vector values needed by elements
  discret_->ClearState();
  // extended SetState(0,...) in case of multiple dofsets (e.g. TSI)
  discret_->SetState(0, "residual displacement", zeros_);
  discret_->SetState(0, "displacement", (*dis_)(0));
  discret_->SetState(0, "velocity", (*vel_)(0));

  // The acceleration is only used as a dummy here and should not be applied inside an element,
  // since this is not the consistent initial acceleration vector which will be determined later on
  discret_->SetState(0, "acceleration", acc_aux);

  discret_->Evaluate(p, stiff_aux, mass_, acc_aux, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // finish mass matrix
  mass_->Complete();

  // calculate intertial force
  mass_->Multiply(false, *acct_, *finertt_);

  return;
}


/*----------------------------------------------------------------------*/
/* check whether the initial conditions are fulfilled */
void STR::TimIntOneStepThetaImmersed::NonlinearMassSanityCheck(
    Teuchos::RCP<const Epetra_Vector> fext, Teuchos::RCP<const Epetra_Vector> dis,
    Teuchos::RCP<const Epetra_Vector> vel, Teuchos::RCP<const Epetra_Vector> acc,
    const Teuchos::ParameterList* sdynparams) const
{
  double fextnorm;
  fext->Norm2(&fextnorm);

  if (fextnorm > 1.0e-14)
  {
    dserror(
        "Initial configuration does not fulfill equilibrium, check your "
        "initial external forces, velocities and accelerations!!!");
  }

  return;
}
