/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
TOPOPT::Optimizer::Optimizer(
    RCP<DRT::Discretization> discret,
    const ParameterList& params
) :
discret_(discret),
params_(params)
{
  const Teuchos::ParameterList& optimizer_params = params_.sublist("TOPOLOGY OPTIMIZER");

  // topology density fields
  dens_i_ = rcp(new Epetra_Vector(*discret_->NodeRowMap(),false));
  dens_ip_ = rcp(new Epetra_Vector(*discret_->NodeRowMap(),false));

  // TODO segmentation fault???
  // fluid fields for optimization
  vel_ = rcp(new map<int,Teuchos::RCP<Epetra_Vector> >);
  adjointvel_ = rcp(new map<int,Teuchos::RCP<Epetra_Vector> >);

  // set initial density field if present
  SetInitialDensityField(DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialDensityField>(optimizer_params,"INITIALFIELD"),
      optimizer_params.get<int>("INITFUNCNO"));

  // set parameters for the objective
  dissipation_ = (DRT::INPUT::IntegralValue<int>(params_,"OBJECTIVE_DISSIPATION")==1);
  pressure_inlet_ = (DRT::INPUT::IntegralValue<int>(params_,"OBJECTIVE_INLET_PRESSURE")==1);
  pressure_drop_ = (DRT::INPUT::IntegralValue<int>(params_,"OBJECTIVE_PRESSURE_DROP")==1);

  dissipation_fac_ = params_.get<double>("DISSIPATION_FAC");
  pressure_inlet_fac_ = params_.get<double>("PRESSURE_INLET_FAC");
  pressure_drop_fac_ = params_.get<double>("PRESSURE_DROP_FAC");

  // general data
  if (DRT::INPUT::IntegralValue<int>(params_,"IS_STATIONARY")==1)
    num_timesteps_ = 1;
  else
    num_timesteps_ = params_.get<int>("NUMSTEP");


  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double TOPOPT::Optimizer::ComputeObjectiveValue(
    Teuchos::RCP<Epetra_Vector> porosity
)
{
  double value = 0.0;

  DataComplete();

  if (dissipation_)
    value += EvaluateObjective(porosity,"dissipation");

  if (pressure_inlet_)
    value += EvaluateObjective(porosity,"pressure inlet");

  if (pressure_drop_)
    value += EvaluateObjective(porosity,"pressure drop");

  return value;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::SetInitialDensityField(
    const INPAR::TOPOPT::InitialDensityField initfield,
    const int startfuncno
)
{
  switch (initfield)
  {
  case INPAR::TOPOPT::initdensfield_zero_field:
  {
    dens_ip_->PutScalar(0.0);
    break;
  }
  case INPAR::TOPOPT::initdensfield_field_by_function:
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // evaluate component k of spatial function
      double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(0,lnode->X(),0.0,NULL); // scalar
      dens_ip_->ReplaceMyValues(1,&initialval,&lnodeid); // lnodeid = ldofid
    }
    break;
  }
  default:
    dserror("unknown initial field");
  }

  dens_i_->Update(1.0,*dens_ip_,0.0);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportFluidData(
    RCP<Epetra_Vector> vel,
    int step)
{
  vel_->insert(pair<int,RCP<Epetra_Vector> >(step,vel));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportAdjointFluidData(
    RCP<Epetra_Vector> vel,
    int step)
{
  adjointvel_->insert(pair<int,RCP<Epetra_Vector> >(step,vel));
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::DataComplete() const
{
  if (num_timesteps_!=vel_->size())
    return false;

  if (num_timesteps_!=adjointvel_->size())
    return false;

  return true;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* TOPOPT::Optimizer::RowMap()
{
  return discret_->NodeRowMap();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* TOPOPT::Optimizer::ColMap()
{
  return discret_->NodeColMap();
}

