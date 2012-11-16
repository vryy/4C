/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.cpp
\brief Static Prestress analysis

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_prestress.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_pstream.H"

/*======================================================================*/
/* constructor */
STR::TimIntPrestress::TimIntPrestress
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<LINALG::Solver> contactsolver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimIntStatics
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output
  )
{
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntPrestress::UpdateStepElement()
{
  // create the parameters for the discretization
  ParameterList p;
  
  // which prestress type?
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
  double pstime = sdyn.get<double>("PRESTRESSTIME");

  // MULF
  if (pstype == INPAR::STR::prestress_mulf)
  {
    if ( (*time_)[0] <= pstime)
    {
      if (!discret_->Comm().MyPID()) IO::cout << "====== Entering MULF update" << IO::endl;
      // action for elements
      p.set("action", "calc_struct_prestress_update");
      discret_->ClearState();
      discret_->SetState(0,"residual displacement",zeros_);
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
    if ( (*time_)[0] <= pstime)
    {
      if (!discret_->Comm().MyPID()) IO::cout << "====== Entering INVERSEDESIGN update" << IO::endl;
      // action for elements
      p.set("action","calc_struct_inversedesign_update");
    }
    else if ( (*time_)[0] > pstime)
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
    p.set("action","calc_struct_inversedesign_switch");
    p.set("total time",timen_ );
    discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  if (pstype == INPAR::STR::prestress_mulf &&  (*time_)[0] <= pstime)
  {
    // only for MULF prestressing mode:
    dis_->UpdateSteps(*zeros_); 
    vel_->UpdateSteps(*zeros_);  // this simply copies zero vectors
    acc_->UpdateSteps(*zeros_);  // this simply copies zero vectors
  }
  
}

/*----------------------------------------------------------------------*/

