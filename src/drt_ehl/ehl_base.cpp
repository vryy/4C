/*--------------------------------------------------------------------------*/
/*!
\file ehl_base.cpp

\brief base class for all elastohydrodynamic lubrication (lubrication structure interaction) algorithms

\level 3

<pre>
\maintainer Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "ehl_base.H"

#include "ehl_partitioned.H"
#include "ehl_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_lubrication.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset_aux_proxy.H"

#include "../drt_lubrication/lubrication_timint_implicit.H"

#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 | constructor                                     (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Base::Base(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string lubrication_disname):
    AlgorithmBase(comm, globaltimeparams),
    structure_(Teuchos::null),
    lubrication_(Teuchos::null),
    fieldcoupling_(DRT::INPUT::IntegralValue<INPAR::EHL::FieldCoupling>(DRT::Problem::Instance()->ElastoHydroDynamicParams(),"FIELDCOUPLING"))
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // get the solver number used for Lubrication solver
  const int linsolvernumber = lubricationparams.get<int>("LINEAR_SOLVER");

  //2.- Setup discretizations and coupling.
  SetupDiscretizations(comm,struct_disname, lubrication_disname);
  SetupFieldCoupling(struct_disname, lubrication_disname);

  //3.- Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis(struct_disname);

  // todo: what to do in the moving grid case?
  bool isale = false;

  // determine which time params to use to build the single fields
  // in case of time stepping time params have to be read from single field sections
  // in case of equal timestep size for all fields the time params are controlled solely
  // by the problem section (e.g. ehl or cell dynamic)
  const Teuchos::ParameterList* structtimeparams = &globaltimeparams;
  const Teuchos::ParameterList* lubricationtimeparams = &globaltimeparams;
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ElastoHydroDynamicParams(),"DIFFTIMESTEPSIZE"))
  {
    structtimeparams = &structparams;
    lubricationtimeparams = &lubricationparams;
  }

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(*structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::Structure>(structure->StructureField());
  structure_->Setup();
  lubrication_ =  Teuchos::rcp(new ADAPTER::LubricationBaseAlgorithm());
  lubrication_->Setup(*lubricationtimeparams,lubricationparams,problem->SolverParams(linsolvernumber),lubrication_disname,isale);

}

/*----------------------------------------------------------------------*
 | read restart information for given time step   (public) wirtz 12/15  |
 *----------------------------------------------------------------------*/
void EHL::Base::ReadRestart( int restart )
{

  if (restart)
  {
    lubrication_->LubricationField()->ReadRestart(restart);
    structure_->ReadRestart(restart);
    SetTimeStep(structure_->TimeOld(), restart);
  }

  return;
}

/*----------------------------------------------------------------------*
 | calculate velocities by a FD approximation               wirtz 12/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> EHL::Base::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> dispnp
  )
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector( *(structure_->Dispn()) ) );
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1./Dt(), *dispnp, -1./Dt());

  return vel;
}  // CalcVelocity()

/*----------------------------------------------------------------------*
 | read restart information for given time        (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::ReadRestartfromTime( double restarttime )
{
  if ( restarttime > 0.0 )
  {

    int restartstructure = EHL::Utils::CheckTimeStepping(structure_->Dt(), restarttime);
    int restartlubrication = EHL::Utils::CheckTimeStepping(lubrication_->LubricationField()->Dt(), restarttime);

    lubrication_->LubricationField()->ReadRestart(restartlubrication);
    structure_->ReadRestart(restartstructure);
    SetTimeStep(structure_->TimeOld(), restartstructure);

  }

  return;
}

/*----------------------------------------------------------------------*
 | test results (if necessary)                     (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(lubrication_->CreateLubricationFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string lubrication_disname)
{
  // Scheme   : the structure discretization is received from the input. Then, an ale-lubrication disc. is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> lubricationdis = problem->GetDis(lubrication_disname);
  if(!structdis->Filled())
    structdis->FillComplete();
  if(!lubricationdis->Filled())
    lubricationdis->FillComplete();

// todo: check if the auxiliary dofsets (proxy) are the things we need

  //first call FillComplete for single discretizations.
  //This way the physical dofs are numbered successively
  structdis->FillComplete();
  lubricationdis->FillComplete();

  //build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_lubrication = lubricationdis->NumDof(0,
      lubricationdis->lRowNode(0));
  const int ndofperelement_lubrication = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux_lubrication =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(ndofpernode_lubrication,ndofperelement_lubrication, 0, true));
  if ( structdis->AddDofSet(dofsetaux_lubrication)!= 1 )
    dserror("unexpected dof sets in structure field");

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux_struct =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct,ndofperelement_struct, 0, true));
  if ( lubricationdis->AddDofSet(dofsetaux_struct)!= 1 )
    dserror("unexpected dof sets in lubrication field");

  //call AssignDegreesOfFreedom also for auxiliary dofsets
  //note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. lubrication dofs
  // 3. structure auxiliary dofs
  // 4. lubrication auxiliary dofs
  structdis->FillComplete(true, false, false);
  lubricationdis->FillComplete(true, false, false);

}

// todo: the following member functions are not yet used and not implemented

/*----------------------------------------------------------------------*
 | set structure solution on lubrication field              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                                       Teuchos::RCP<const Epetra_Vector> vel )
{
  // todo: here, we have to think about what quantity is relavant for the lubrication framework to get from the structure framework
  SetMeshDisp(disp);
  SetVelocityFields(vel);
}

/*----------------------------------------------------------------------*
 | set lubrication solution on structure field              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetLubricationSolution( Teuchos::RCP<const Epetra_Vector> quantity )
{
  // todo: here, we have to think about what quantity is relavant for the structure framework to get from the lubrication framework
}


/*----------------------------------------------------------------------*
 | set structure velocity fields on lubrication field       wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetVelocityFields( Teuchos::RCP<const Epetra_Vector> quantity)
{
  switch(fieldcoupling_)
  {
  case INPAR::EHL::coupling_none:
    break;
  default:
    dserror("unknown field coupling type in SetVelocityFields()");
    break;
  }
}

/*----------------------------------------------------------------------*
 | set structure mesh displacement on lubrication field     wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp )
{
  switch(fieldcoupling_)
  {
  case INPAR::EHL::coupling_none:
    break;
  default:
    dserror("unknown field coupling type in SetMeshDisp()");
    break;
  }
}

/*----------------------------------------------------------------------*
 | setup adapters for EHL on boundary                       wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetupFieldCoupling(const std::string struct_disname, const std::string lubrication_disname)
{

  // todo: here, we need to setup a coupling adapter

}

/*----------------------------------------------------------------------*
 | update (protected)                                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Base::Update()
{
  StructureField()->Update();
  lubrication_->LubricationField()->Update();

  return;
}

/*----------------------------------------------------------------------*
 | output (protected)                                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Base::Output(bool forced_writerestart)
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.

  //===========================
  // output for structurefield:
  //===========================
//  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());
  StructureField()->Output(forced_writerestart);

  //=============================
  // output for lubricationfield:
  //=============================
//  ApplyStructCouplingState(StructureField()->Dispnp(),StructureField()->Velnp());
  lubrication_->LubricationField()->Output(forced_writerestart);

  //reset states
  StructureField()->Discretization()->ClearState(true);
  lubrication_->LubricationField()->Discretization()->ClearState(true);
}  // Output()
