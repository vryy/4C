/*----------------------------------------------------------------------*/
/*!
\brief Static Prestress analysis

\level 2

\maintainer Fabian Braeu
*/
/*----------------------------------------------------------------------*/

/* headers */
#include "strtimint_prestress.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"
#include "../drt_constraint/springdashpot_manager.H"

/*======================================================================*/
/* constructor */
STR::TimIntPrestress::TimIntPrestress(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioparams, const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<LINALG::Solver> contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter> output)
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

/*----------------------------------------------------------------------------------------------*
 * Initialize this class                                                            rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntPrestress::Init(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver)
{
  // call Init() in base class
  STR::TimIntStatics::Init(timeparams, sdynparams, xparams, actdis, solver);

  return;
}

/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntPrestress::Setup()
{
  // call Setup() in base class
  STR::TimIntStatics::Setup();

  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntPrestress::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // which prestress type?
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  INPAR::STR::PreStress pstype =
      DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
  double pstime = sdyn.get<double>("PRESTRESSTIME");

  // MULF
  if (pstype == INPAR::STR::prestress_mulf)
  {
    if ((*time_)[0] <= pstime)
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
  else if (pstype == INPAR::STR::prestress_id)
  {
    if ((*time_)[0] <= pstime)
    {
      if (!discret_->Comm().MyPID()) IO::cout << "====== Entering INVERSEDESIGN update" << IO::endl;
      // action for elements
      p.set("action", "calc_struct_inversedesign_update");
    }
    else if ((*time_)[0] > pstime)
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


  if (pstype == INPAR::STR::prestress_id && (*time_)[0] <= pstime && timen_ > pstime)
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

  if (pstype == INPAR::STR::prestress_mulf && (*time_)[0] <= pstime)
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
/* write internal and external forces for restart */
/* only necessary for dynamic restart! */
void STR::TimIntPrestress::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  output->WriteVector("fexternal", fextn_);
  output->WriteVector("fint", fintn_);
  output->WriteVector("finert", zeros_);
  return;
}

/*----------------------------------------------------------------------*/
