/*----------------------------------------------------------------------*/
/*! \file
\brief Static Prestress analysis

\level 2

*/
/*----------------------------------------------------------------------*/

/* headers */
#include "strtimint_prestress.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"
#include "../drt_lib/prestress_service.H"
#include "../drt_constraint/springdashpot_manager.H"

/*======================================================================*/
/* constructor */
STR::TimIntPrestress::TimIntPrestress(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioparams, const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams, const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<LINALG::Solver>& contactsolver,
    const Teuchos::RCP<IO::DiscretizationWriter>& output)
    : TimIntStatics(
          timeparams, ioparams, sdynparams, xparams, actdis, solver, contactsolver, output)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntPrestress::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // MULF
  if (::UTILS::PRESTRESS::IsMulf())
  {
    if (::UTILS::PRESTRESS::IsMulfActive((*time_)[0]))
    {
      if (!discret_->Comm().MyPID()) IO::cout << "====== Entering MULF update" << IO::endl;
      // action for elements
      p.set("action", "calc_struct_prestress_update");
      discret_->ClearState();
      discret_->SetState(0, "residual displacement", zeros_);
    }
    else
    {
      // action for elements
      p.set("action", "calc_struct_update_istep");
      discret_->ClearState();
    }
  }

  // INVERSE DESIGN
  else if (::UTILS::PRESTRESS::IsInverseDesign())
  {
    if (::UTILS::PRESTRESS::IsInverseDesignActive((*time_)[0]))
    {
      if (!discret_->Comm().MyPID()) IO::cout << "====== Entering INVERSEDESIGN update" << IO::endl;
      // action for elements
      p.set("action", "calc_struct_inversedesign_update");
    }
    else
    {
      // action for elements
      p.set("action", "calc_struct_update_istep");
      discret_->ClearState();
    }
  }

  // params for both MULF and ID
  p.set("total time", (*time_)[0]);
  p.set("delta time", (*dt_)[0]);

  // go to elements
  discret_->SetState("displacement", (*dis_)(0));
  discret_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);


  if (::UTILS::PRESTRESS::IsInverseDesignActive((*time_)[0]) &&
      !::UTILS::PRESTRESS::IsInverseDesignActive(timen_))
  {
    // switch in id mode:
    dis_->UpdateSteps(*zeros_);
    vel_->UpdateSteps(*zeros_);  // this simply copies zero vectors
    acc_->UpdateSteps(*zeros_);  // this simply copies zero vectors
    if (!discret_->Comm().MyPID()) IO::cout << "XXXXXX Entering INVERSEDESIGN SWITCH" << IO::endl;
    // action for elements
    p.set("action", "calc_struct_inversedesign_switch");
    p.set("total time", timen_);
    discret_->Evaluate(
        p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }

  if (::UTILS::PRESTRESS::IsMulfActive((*time_)[0]))
  {
    // prestressing for spring in spring dashpot - corresponds to storage of deformation gradient in
    // material law (mhv 12/2015) pass current displacement state to spring at end of MULF step
    if (springman_->HaveSpringDashpot())
    {
      springman_->ResetPrestress(disn_);
    }
    // only for MULF prestressing mode:
    dis_->UpdateSteps(*zeros_);
    vel_->UpdateSteps(*zeros_);  // this simply copies zero vectors
    acc_->UpdateSteps(*zeros_);  // this simply copies zero vectors
  }
}

/*----------------------------------------------------------------------*/
