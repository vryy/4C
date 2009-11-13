/*----------------------------------------------------------------------*/
/*!
\file adapter_thermo_timint.cpp

\brief Thermo field adapter

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                             bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 |  headers                                                 bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "adapter_thermo_timint.H"

#include "../drt_thermo/thrtimint_impl.H"
#include "../drt_thermo/thrtimint_statics.H"
#include "../drt_thermo/thrtimint_ost.H"
#include "../drt_thermo/thrtimint_genalpha.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_thermo/thr_resulttest.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for ThermoBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |  constructor                                             bborn 08/09 |
 *----------------------------------------------------------------------*/
ADAPTER::ThermoTimInt::ThermoTimInt(
  Teuchos::RCP<Teuchos::ParameterList> ioparams,
  Teuchos::RCP<Teuchos::ParameterList> tdynparams,
  Teuchos::RCP<Teuchos::ParameterList> xparams,
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
  )
: thermo_(Create(*ioparams, *tdynparams, *xparams, discret, solver, output)),
  discret_(discret),
  ioparams_(ioparams),
  tdynparams_(tdynparams),
  xparams_(xparams),
  solver_(solver),
  output_(output)
{
  // make sure
  if (thermo_ == Teuchos::null)
    dserror("Failed to create thermal integrator");

  // initialise temperature increments to 0 (in words zero)
  tempinc_ = Teuchos::rcp(new Epetra_Vector(*(DofRowMap()),true));

  // good bye
  return;
}

/*----------------------------------------------------------------------*
 | create implicit marching time integrator                 bborn 08/09 |
 | originally included in the file strtimint_create.cpp                 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<THR::TimIntImpl> ADAPTER::ThermoTimInt::Create(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& tdyn,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization>& actdis,
  Teuchos::RCP<LINALG::Solver>& solver,
  Teuchos::RCP<IO::DiscretizationWriter>& output
  )
{
  Teuchos::RCP<THR::TimIntImpl> tti = Teuchos::null;

  // create specific time integrator
  switch (Teuchos::getIntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"))
  {
  // Static analysis
  case INPAR::THR::dyna_statics :
  {
    tti = Teuchos::rcp(new THR::TimIntStatics(ioflags, tdyn, xparams,
                                              actdis, solver, output));
    break;
  }

  // One-step-theta (OST) time integration
  case INPAR::THR::dyna_onesteptheta :
  {
    tti = Teuchos::rcp(new THR::TimIntOneStepTheta(ioflags, tdyn, xparams,
                                                   actdis, solver, output));
    break;
  }

  /*
  // Generalised energy-momentum method (GEMM) time integration
   case INPAR::THR::dyna_gemm :
   {
     tti = Teuchos::rcp(new THR::TimIntGEMM(ioflags, tdyn, xparams,
                                                    actdis, solver, output));
     break;
    }
*/
  // Generalized alpha time integration
  case INPAR::THR::dyna_genalpha :
  {
    tti = Teuchos::rcp(new THR::TimIntGenAlpha(ioflags, tdyn, xparams,
                                               actdis, solver, output));
    break;
  }

  // Everything else
  default :
  {
    // do nothing
    break;
  }
  } // switch

  // return the integrator
  return tti;
}

/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimInt::InitialGuess()
{
  return thermo_->TempRes();
}

/*----------------------------------------------------------------------*
 | right-hand side alias the dynamic force residual         bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimInt::RHS()
{
  // this expects a _negative_ (Newton-ready) residual with blanked
  // Dirichlet DOFs. We did it in #Evaluate.
  return thermo_->ForceRes();
}

/*----------------------------------------------------------------------*
 | get current temperature T_{n+1}                          bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimInt::Tempnp()
{
  return thermo_->TempNew();
}

/*----------------------------------------------------------------------*
 | get last converged temperature T_{n}                     bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimInt::Tempn()
{
  return thermo_->Temp();
}

/*----------------------------------------------------------------------*
 | non-overlapping DOF map                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::ThermoTimInt::DofRowMap()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*
 | tangent, i.e. force residual R_{n+1} differentiated      bborn 08/09 |
 | by temperatures T_{n+1}                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::ThermoTimInt::SystemMatrix()
{
  return thermo_->Tang();
}

/*----------------------------------------------------------------------*
 | get discretisation                                       bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::ThermoTimInt::Discretization()
{
  return thermo_->Discretization();
}

/*----------------------------------------------------------------------*
 | External force F_{ext,n+1}                               bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimInt::FExtnp()
{
  return thermo_->FextNew();
}

/*----------------------------------------------------------------------*
 | prepare time step                                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimInt::PrepareTimeStep()
{
  // Note: MFSI requires a constant predictor. Otherwise the fields will get
  // out of sync.

  // predict
  thermo_->Predict();

  // initialise incremental temperatures
  tempinc_->PutScalar(0.0);

}

/*----------------------------------------------------------------------*
 | build linear system tangent matrix and rhs/force residual bborn 08/09 |
 | Monolithic TSI accesses the linearised thermo problem                |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimInt::Evaluate(
  Teuchos::RCP<const Epetra_Vector> temp
  )
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (temp != Teuchos::null)
  {
    // residual temperatures (or iteration increments or iteratively incremental temperatures)
    Teuchos::RCP<Epetra_Vector> tempi = Teuchos::rcp(new Epetra_Vector(*temp));
    tempi->Update(-1.0, *tempinc_, 1.0);

    // update incremental temperature member to provided step increments
    // shortly: tempinc_^<i> := temp^<i+1>
    tempinc_->Update(1.0, *temp, 0.0);

    // do thermal update with provided residual temperatures
    thermo_->UpdateIterIncrementally(tempi);
  }
  else
  {
    thermo_->UpdateIterIncrementally(Teuchos::null);
  }

  // builds tangent, residual and applies DBC
  thermo_->EvaluateRhsTangResidual();
  thermo_->PrepareSystemForNewtonSolve();
}

/*----------------------------------------------------------------------*
 | update time step                                         bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimInt::Update()
{
  thermo_->UpdateStepState();
  thermo_->UpdateStepTime();
  return;
}

/*----------------------------------------------------------------------*
 | output                                                   bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimInt::Output()
{
  thermo_->OutputStep();
}

/*----------------------------------------------------------------------*
 | domain map                                               bborn 08/09 |
 *----------------------------------------------------------------------*/
const Epetra_Map& ADAPTER::ThermoTimInt::DomainMap()
{
  return thermo_->GetDomainMap();
}

/*----------------------------------------------------------------------*
 | read restart                                             bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimInt::ReadRestart(const int step)
{
  thermo_->ReadRestart(step);
}

/*----------------------------------------------------------------------*
 | find iteratively solution                                bborn 08/09 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimInt::Solve()
{
  thermo_->Solve();
}

/*----------------------------------------------------------------------*
 | thermal result test                                      bborn 08/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ThermoTimInt::CreateFieldTest()
{
  return Teuchos::rcp(new THR::ResultTest(*thermo_));
}

/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
