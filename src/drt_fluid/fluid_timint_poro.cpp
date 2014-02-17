/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_poro.cpp
\brief TimIntPoro

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_poro.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntPoro::TimIntPoro(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid)
{

  Teuchos::ParameterList *  stabparams;
  if (discret_->Name()=="fluid" and DRT::Problem::Instance()->ProblemType()==prb_fpsi)
    stabparams=&(params_->sublist("RESIDUAL-BASED STABILIZATION"));
  else if (discret_->Name()=="porofluid")
    stabparams=&(params_->sublist("POROUS-FLOW STABILIZATION"));
  else
    stabparams=&(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if(stabparams->get<std::string>("STABTYPE")=="residual_based")
    if(stabparams->get<std::string>("TDS") == "time_dependent")
      dserror("TDS is not implemented for Poro yet. An error will occur in FluidImplicitTimeInt::TimeUpdate().");

  if (alefluid_)
  {
    if (  (physicaltype_ == INPAR::FLUID::poro
        or physicaltype_ == INPAR::FLUID::poro_p1
        or physicaltype_ == INPAR::FLUID::poro_p2)
        and discret_->Name()=="porofluid" )
      //gridvn_ can also be moved to poro class?
      gridvn_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  }

  //set some poro-specific parameters
  if (  (physicaltype_ == INPAR::FLUID::poro
      or physicaltype_ == INPAR::FLUID::poro_p1
      or physicaltype_ == INPAR::FLUID::poro_p2)
      and discret_->Name()=="porofluid")
    SetElementCustomParameter();
  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntPoro::~TimIntPoro()
{
  return;
}

/*----------------------------------------------------------------------*
| set params in constructor                                    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetElementGeneralFluidParameter()
{

  //set some poro-specific parameters only in specific poro cases
  if (  not ((physicaltype_ == INPAR::FLUID::poro
      or physicaltype_ == INPAR::FLUID::poro_p1
      or physicaltype_ == INPAR::FLUID::poro_p2)
      and discret_->Name()=="porofluid"))
  {
    FluidImplicitTimeInt::SetElementGeneralFluidParameter();
  }
  return;
}

/*----------------------------------------------------------------------*
| set params in constructor                                    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetElementTurbulenceParameter()
{

  //set some poro-specific parameters only in specific poro cases
  if (  not ((physicaltype_ == INPAR::FLUID::poro
      or physicaltype_ == INPAR::FLUID::poro_p1
      or physicaltype_ == INPAR::FLUID::poro_p2)
      and discret_->Name()=="porofluid"))
  {
    FluidImplicitTimeInt::SetElementTurbulenceParameter();
  }
  return;
}

void FLD::TimIntPoro::AssembleMatAndRHS()
{
  FluidImplicitTimeInt::AssembleMatAndRHS();
  PoroIntUpdate();

  return;
}

void FLD::TimIntPoro::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);

  FluidImplicitTimeInt::ReadRestart(step);

  if(alefluid_)
  {
    if((physicaltype_ == INPAR::FLUID::poro or physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2) and discret_->Name()=="porofluid" )
      reader.ReadVector(gridv_,"gridv");
    if((physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2) and discret_->Name()=="porofluid")
      reader.ReadVector(gridvn_,"gridvn");
  }
  return;
}

// -------------------------------------------------------------------
// set poro parameters                               vuong  11/2012
// -------------------------------------------------------------------
void FLD::TimIntPoro::SetElementCustomParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_poro_parameter);

  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // set poro specific element parameters
  eleparams.set<bool>("conti partial integration",params_->get<bool>("conti partial integration"));

  // parameter for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");
  eleparams.sublist("EDGE-BASED STABILIZATION") = params_->sublist("EDGE-BASED STABILIZATION");
  eleparams.sublist("POROUS-FLOW STABILIZATION") = params_->sublist("POROUS-FLOW STABILIZATION");

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
 |  set initial field for porosity                                      |
 *----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetInitialPorosityField(
    const INPAR::POROELAST::InitialField init,
    const int startfuncno)
{
  std::cout<<"FLD::FluidImplicitTimeInt::SetInitialPorosityField()"<<std::endl;

  switch(init)
  {
  case INPAR::POROELAST::initfield_field_by_function:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = discret_->Dof(lnode);

      int numdofs = nodedofset.size();
      double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(0,lnode->X(),time_,NULL);

      // check whether there are invalid values of porosity
      if (initialval < EPS15) dserror("zero or negative initial porosity");
      if (initialval >= 1) dserror("initial porosity greater or equal than 1");
      for (int k=0;k< numdofs;++k)
      {
        const int dofgid = nodedofset[k];
        int doflid = dofrowmap->LID(dofgid);
        // evaluate component k of spatial function
        int err = initporosityfield_->ReplaceMyValues(1,&initialval,&doflid);
        if (err != 0) dserror("dof not on proc");

      }
    }

    break;
  }
  default:
    dserror("Unknown option for initial field: %d", init);
    break;
  } // switch(init)

  return;
} // TimIntPoro::SetInitialPorosityField

/*----------------------------------------------------------------------*
| add some functionality to UpdateIterIncrementally            bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::UpdateIterIncrementally(
  Teuchos::RCP<const Epetra_Vector> vel)  //!< input residual velocities

{
  FluidImplicitTimeInt::UpdateIterIncrementally(vel);
  // set the new solution we just got
  if (vel != Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(
        *(discret_->DofRowMap(0)), true);


    if ((physicaltype_ == INPAR::FLUID::poro or physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2) and discret_->Name()=="porofluid" )
    {
      //only one step theta
      // new end-point accelerations
      aux->Update(1.0 / (theta_ * dta_), *velnp_, -1.0 / (theta_ * dta_),
          *(*veln_)(0), 0.0);
      aux->Update(-(1.0 - theta_) / theta_, *(*accn_)(0), 1.0);
      // put only to free/non-DBC DOFs
      dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(accnp_), aux);
      *accnp_ = *aux;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 | overloading function                                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntPoro::Output()
{

  FluidImplicitTimeInt::Output();
  // output of solution
  if (step_%upres_ == 0)
  {
    if((physicaltype_ == INPAR::FLUID::poro or physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2) and discret_->Name()=="porofluid" )
    {
      Teuchos::RCP<Epetra_Vector>  convel= Teuchos::rcp(new Epetra_Vector(*velnp_));
      convel->Update(-1.0,*gridv_,1.0);
      output_->WriteVector("convel", convel);
      output_->WriteVector("gridv", gridv_);
      if(physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2)
        output_->WriteVector("gridvn", gridvn_);
    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ > 0 && step_%uprestart_ == 0)
  {
    if(physicaltype_ == INPAR::FLUID::poro or physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2)
      output_->WriteVector("gridv", gridv_);
    if(physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2)
      output_->WriteVector("gridvn", gridvn_);
  }
  return;
} // TimIntPoro::Output

/*----------------------------------------------------------------------*
| set params in constructor                                    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  if (discret_->Name()=="porofluid")
  {
      std::string ptype = DRT::Problem::Instance()->FluidDynamicParams().get<std::string>("PHYSICAL_TYPE");
      if(ptype == (std::string)"Poro_P1")
        physicaltype_=INPAR::FLUID::poro_p1;
      else if(ptype == (std::string)"Poro_P2")
        physicaltype_=INPAR::FLUID::poro_p2;
      else if(ptype == (std::string)"Poro")
        physicaltype_=INPAR::FLUID::poro;
  }
  else if (discret_->Name()=="fluid" and DRT::Problem::Instance()->ProblemType() == prb_fpsi)
  {
    physicaltype_ = INPAR::FLUID::incompressible;
  }
  else
    dserror("Poro but neither 'porofluid' nor 'fluid and prb_fpsi'");

  eleparams.set<int>("physical type",physicaltype_);

  if (alefluid_)
  {
    if (   physicaltype_ == INPAR::FLUID::poro
        or physicaltype_ == INPAR::FLUID::poro_p1
        or physicaltype_ == INPAR::FLUID::poro_p2)
    {
      //just for poroelasticity
      discret_->SetState("dispn", dispn_);
      discret_->SetState("accnp", accnp_);
      discret_->SetState("accn", accn_);
      discret_->SetState("gridvn", gridvn_);

      eleparams.set("total time", time_);
      eleparams.set("delta time", dta_);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
| update sysmat after AssembleMatAndRHS                        bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::PoroIntUpdate()
{
  if (   physicaltype_ == INPAR::FLUID::poro
        or physicaltype_ == INPAR::FLUID::poro_p1
        or physicaltype_ == INPAR::FLUID::poro_p2)
  {
    sysmat_->UnComplete();

    std::string condname = "PoroPartInt";
    std::vector<DRT::Condition*> poroPartInt;
    discret_->GetCondition(condname,poroPartInt);
    if(poroPartInt.size())
    {
      Teuchos::ParameterList eleparams;

      // set action for elements
      eleparams.set<int>("action",FLD::poro_boundary);
      eleparams.set("total time", time_);
      eleparams.set("delta time", dta_);
      eleparams.set<POROELAST::coupltype>("coupling",POROELAST::fluidfluid);
      eleparams.set<int>("physical type",physicaltype_);

      discret_->ClearState();
      discret_->SetState("dispnp", dispnp_);
      discret_->SetState("gridv", gridv_);
      discret_->SetState("velnp",velnp_);
      discret_->SetState("scaaf",scaaf_);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
      discret_->ClearState();
    }

    condname = "PoroPresInt";
    std::vector<DRT::Condition*> poroPresInt;
    discret_->GetCondition(condname,poroPresInt);
    if(poroPresInt.size())
    {
      Teuchos::ParameterList eleparams;

      // set action for elements
      eleparams.set<int>("action",FLD::poro_prescoupl);
      eleparams.set<POROELAST::coupltype>("coupling",POROELAST::fluidfluid);
      eleparams.set<int>("physical type",physicaltype_);

      discret_->ClearState();
      discret_->SetState("dispnp", dispnp_);
      discret_->SetState("gridv", gridv_);
      discret_->SetState("velnp",velnp_);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
      discret_->ClearState();
    }
    sysmat_->Complete();
  }

  return;
}

/*----------------------------------------------------------------------*
| Calculate acceleration for poro                              bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::TimIntCalculateAcceleration()
{
  Teuchos::RCP<Epetra_Vector> onlyaccn = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> onlyvelnm = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> onlyveln = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> onlyvelnp = Teuchos::null;

  if (not (    physicaltype_ == INPAR::FLUID::poro
            or physicaltype_ == INPAR::FLUID::poro_p1
            or physicaltype_ == INPAR::FLUID::poro_p2
          )
            or
          (    DRT::Problem::Instance()->ProblemType()==prb_fpsi
               and discret_->Name()=="fluid"
          )
     ) //standard case
  {
    onlyaccn  = velpressplitter_.ExtractOtherVector(accn_);
    onlyaccnp = velpressplitter_.ExtractOtherVector(accnp_);
    onlyvelnm = velpressplitter_.ExtractOtherVector(velnm_);
    onlyveln  = velpressplitter_.ExtractOtherVector(veln_);
    onlyvelnp = velpressplitter_.ExtractOtherVector(velnp_);
  }
  else //poroelasticity case
  {
    onlyaccn = accn_;
    onlyaccnp = accnp_;
    onlyvelnm = velnm_;
    onlyveln = veln_;
    onlyvelnp = velnp_;
  }

  CalculateAcceleration(onlyvelnp,
                        onlyveln ,
                        onlyvelnm,
                        onlyaccn ,
                        onlyaccnp);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*accnp_);

  if ( physicaltype_ == INPAR::FLUID::poro_p1 or physicaltype_ == INPAR::FLUID::poro_p2)
    gridvn_ ->Update(1.0,*gridv_,0.0);
  return;
}
