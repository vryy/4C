/*----------------------------------------------------------------------*/
/*!
\file adapter_thermo_timint_expl.cpp

\brief Thermo field adapter

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 01/12 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 |  headers                                                  dano 01/12 |
 *----------------------------------------------------------------------*/
#include "adapter_thermo_timint_expl.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_thermo/thr_resulttest.H"
#include "../drt_thermo/thrtimint_expleuler.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for ThermoBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |  constructor                                              dano 01/12 |
 *----------------------------------------------------------------------*/
ADAPTER::ThermoTimIntExpl::ThermoTimIntExpl(
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
  // DofRowMap for multiple dofsets
  tempinc_ = Teuchos::rcp(new Epetra_Vector(*(DofRowMap(0)),true));

  // good bye
  return;
}


/*----------------------------------------------------------------------*
 | create implicit marching time integrator                  dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<THR::TimIntExpl> ADAPTER::ThermoTimIntExpl::Create(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& tdyn,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization>& actdis,
  Teuchos::RCP<LINALG::Solver>& solver,
  Teuchos::RCP<IO::DiscretizationWriter>& output
  )
{
  // build the thermal time integration tti
  Teuchos::RCP<THR::TimIntExpl> tti = Teuchos::null;

  // create specific time integrator
  switch (DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"))
  {
    // explicit analysis using forward Euler
    case INPAR::THR::dyna_expleuler :
    {
      tti = Teuchos::rcp(new THR::TimIntExplEuler(ioflags, tdyn, xparams, actdis, solver, output));
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
 |                                                           dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimIntExpl::InitialGuess()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 | right-hand side alias the dynamic force residual          dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimIntExpl::RHS()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 | get current temperature T_{n+1}                           dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimIntExpl::Tempnp()
{
  return thermo_->TempNew();
}


/*----------------------------------------------------------------------*
 | get last converged temperature T_{n}                      dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimIntExpl::Tempn()
{
  return thermo_->Temp();
}


/*----------------------------------------------------------------------*
 | non-overlapping DOF map for multiple dofsets              dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::ThermoTimIntExpl::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------*
 | non-overlapping DOF map                                   dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::ThermoTimIntExpl::DofRowMap()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}


/*----------------------------------------------------------------------*
 | tangent, i.e. force residual R_{n+1} differentiated       dano 01/12 |
 | by temperatures T_{n+1}                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::ThermoTimIntExpl::SystemMatrix()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 | recalculate thermal matrices for tsi simulations          dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::TSIMatrix()
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*
 | get discretisation                                        dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::ThermoTimIntExpl::Discretization()
{
  return thermo_->Discretization();
}


/*----------------------------------------------------------------------*
 | External force F_{ext,n+1}                                dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::ThermoTimIntExpl::FExtnp()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 | Thermal contact manager                                   mgit 07/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<THR::ThermoContactMan> ADAPTER::ThermoTimIntExpl::ThermoContactManager()
{
  return thermo_->ThermoContactManager();
}

/*----------------------------------------------------------------------*
 | prepare time step                                        dano 01/12|
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::PrepareTimeStep()
{
  //do nothing: there is no predictor for explicit time integration
  return;
}


/*----------------------------------------------------------------------*
 | build linear system tangent matrix, rhs/force residual   dano 01/12|
 | Monolithic TSI accesses the linearised thermo problem                |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::Evaluate(Teuchos::RCP<const Epetra_Vector> temp)
{
  dserror("not implemented");
}



/*----------------------------------------------------------------------*
 | update time step                                          dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::Update()
{
  // update temperature and temperature rate
  // after this call we will have tempn_ == temp_ (temp_{n+1} == temp_n), etc.
  thermo_->UpdateStepState();
  // update time and step
  thermo_->UpdateStepTime();
  // currently nothing, can include history dependency of materials
  thermo_->UpdateStepElement();
  return;
}


/*----------------------------------------------------------------------*
 | print step summary                                        dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::PrintStep()
{
  thermo_->PrintStep();
  return;
}


/*----------------------------------------------------------------------*
 | output                                                    dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::Output()
{
  thermo_->OutputStep();
}


/*----------------------------------------------------------------------*
 | domain map                                                dano 01/12 |
 *----------------------------------------------------------------------*/
const Epetra_Map& ADAPTER::ThermoTimIntExpl::DomainMap()
{
  return thermo_->GetDomainMap();
}


/*----------------------------------------------------------------------*
 | read restart                                              dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::ReadRestart(const int step)
{
  thermo_->ReadRestart(step);
}

/*----------------------------------------------------------------------*
 | prepare thermal contact                                   mgit 06/11 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::PrepareThermoContact(Teuchos::RCP<MORTAR::ManagerBase> cmtman,
                                                  Teuchos::RCP<DRT::Discretization> discretstruct)
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*
 | find iteratively solution                                 dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::Solve()
{
  thermo_->IntegrateStep();
}

/*----------------------------------------------------------------------*
 | thermal result test                                       dano 01/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ThermoTimIntExpl::CreateFieldTest()
{
  return Teuchos::rcp(new THR::ResultTest(*thermo_));
}

/*----------------------------------------------------------------------*
 | extract temperature T_{n} for TSI                         dano 01/12 |
 | (named like in FSI)                                                  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ThermoTimIntExpl::ExtractTempn()
{
  dserror("not implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 | extract current temperature T_{n+1} for TSI               dano 01/12 |
 | (named like in FSI)                                                  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ThermoTimIntExpl::ExtractTempnp()
{
  dserror("not implemented");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 | update Newton step                                        dano 02/11 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::UpdateNewton(
  Teuchos::RCP<const Epetra_Vector> temp
  )
{
  dserror("not implemented");
}

/*----------------------------------------------------------------------*
 | build linear system tangent matrix, rhs/force residual    dano 02/11 |
 | Monolithic TSI accesses the linearised thermo problem                |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::Evaluate()
{
  dserror("not implemented");
}


/*----------------------------------------------------------------------*
 | apply current displacements and velocities needed in TSI  dano 01/12 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::ApplyStructVariables(
  Teuchos::RCP<const Epetra_Vector> disp,
  Teuchos::RCP<const Epetra_Vector> vel
  )
{
  // pass current displacements and velocities to the thermo field
  thermo_->ApplyStructVariables(disp,vel);
}

/*----------------------------------------------------------------------*
 | evaluate the residual forces, use the last converged      dano 12/10 |
 | solution for predictor                                               |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::PreparePartitionStep()
{
  dserror("not yet checked in ADAPTER::ThermoTimIntExpl::PreparePartitionStep()!");
//  thermo_->PreparePartitionStep();
}


/*----------------------------------------------------------------------*
 | external interface loads (heat fluxes) are applied to                |
 | the thermo field                                         ghamm 12/10 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::ApplyInterfaceForces
(
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  thermo_->SetForceInterface(ithermoload);
}


/*----------------------------------------------------------------------*
 | map extractor is set to communicate external loads       ghamm 12/10 |
 *----------------------------------------------------------------------*/
void ADAPTER::ThermoTimIntExpl::SetSurfaceTFSI
(
  Teuchos::RCP<const LINALG::MapExtractor> tfsisurface  //!< the TFSI surface
)
{
  thermo_->SetSurfaceTFSI(tfsisurface);
}


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
