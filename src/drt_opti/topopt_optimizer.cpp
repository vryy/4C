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
#include "topopt_fluidAdjoint3_impl_parameter.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
TOPOPT::Optimizer::Optimizer(
    RCP<DRT::Discretization> optidis,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    const ParameterList& params
) :
optidis_(optidis),
fluiddis_(fluiddis),
params_(params)
{
  const Teuchos::ParameterList& optimizer_params = params_.sublist("TOPOLOGY OPTIMIZER");

  // fluid fields for optimization
  fluidvel_ = rcp(new map<int,Teuchos::RCP<Epetra_Vector> >);
  adjointvel_ = rcp(new map<int,Teuchos::RCP<Epetra_Vector> >);

  // topology density fields
  dens_i_ = rcp(new Epetra_Vector(*optidis_->NodeRowMap(),false));
  dens_ip_ = rcp(new Epetra_Vector(*optidis_->NodeRowMap(),false));

  // value of the objective function
  obj_value_ = 0.0;
  // gradient of the objective function
  obj_grad_ = rcp(new Epetra_Vector(*optidis_->NodeRowMap()));

  // set initial density field if present
  SetInitialDensityField(DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialDensityField>(optimizer_params,"INITIALFIELD"),
      optimizer_params.get<int>("INITFUNCNO"));

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeObjective()
{
  // check if all data is present
  DataComplete();

  // Transform fluid and adjoint field so that it is better readable at element level
  TransformFlowFields();

  Teuchos::ParameterList params;

  params.set("action","compute_objective");

  params.set("fluidvel",fluidvel_);
  params.set("adjointvel",adjointvel_);
  params.set("fluiddis",fluiddis_);

  optidis_->ClearState();

  optidis_->SetState("density",dens_ip_);
  optidis_->Evaluate(params,Teuchos::null,obj_grad_);

  optidis_->ClearState();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeGradient()
{
  TEUCHOS_FUNC_TIME_MONITOR(" evaluate objective gradient");
  obj_grad_->PutScalar(0.0);

  DataComplete();

  // Transform fluid and adjoint field so that it is better readable at element level
  TransformFlowFields();

  Teuchos::ParameterList params;

  params.set("action","compute_gradient");

  params.set("objective_value",obj_value_);

  params.set("fluidvel",fluidvel_);
  params.set("adjointvel",adjointvel_);
  params.set("fluiddis",fluiddis_);

  optidis_->ClearState();

  optidis_->SetState("density",dens_ip_);
  optidis_->Evaluate(params,Teuchos::null,Teuchos::null);

  optidis_->ClearState();
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
    for(int lnodeid=0;lnodeid<optidis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = optidis_->lRowNode(lnodeid);

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
void TOPOPT::Optimizer::ImportFlowParams(
    Teuchos::RCP<Teuchos::ParameterList>& fluidParams
)
{
  fluidParams_ = fluidParams;

  // set the general parameter on element level
  Teuchos::ParameterList opti_ele_params;
  opti_ele_params.set("action","set_general_optimization_parameter");

  // flow parameter
  opti_ele_params.set("timealgo",fluidParams_->get<int>("time int algo"));
  opti_ele_params.set("dt",fluidParams_->get<double>("time step size"));
  opti_ele_params.set("maxtimesteps",fluidParams_->get<int>("max number timesteps"));

  // objective parameter
  opti_ele_params.set("dissipation",fluidParams_->get<bool>("OBJECTIVE_DISSIPATION"));
  opti_ele_params.set("pres_drop",fluidParams_->get<bool>("OBJECTIVE_PRESSURE_DROP"));

  opti_ele_params.set("dissipation_fac",fluidParams_->get<double>("DISSIPATION_FAC"));
  opti_ele_params.set("pres_drop_fac",fluidParams_->get<double>("PRESSURE_DROP_FAC"));

  optidis_->Evaluate(opti_ele_params,Teuchos::null,Teuchos::null);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportFluidData(
    RCP<Epetra_Vector> vel,
    int step
)
{
  if (fluidvel_->find(step)!=fluidvel_->end())
    *(fluidvel_->find(step)->second) = *vel;
  else
  {
    // a copy is required here, otherwise the values are changed within the fluid time integration
    RCP<Epetra_Vector> new_vel = rcp(new Epetra_Vector(*vel)); // copy
    fluidvel_->insert(pair<int,RCP<Epetra_Vector> >(step,new_vel));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportAdjointFluidData(
    RCP<Epetra_Vector> vel,
    int step)
{
  if (adjointvel_->find(step)!=adjointvel_->end())
    *(adjointvel_->find(step)->second) = *vel;
  else
  {
    // a copy is required here, otherwise the values are changed within the adjoint time integration
    RCP<Epetra_Vector> new_vel = rcp(new Epetra_Vector(*vel));
    adjointvel_->insert(pair<int,RCP<Epetra_Vector> >(step,new_vel));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::DataComplete() const
{
  // n timesteps
  // -> solutions at time t^0,t^1,...,t^n
  // -> n+1 solutions

  size_t num_sols = 0;
  if (fluidParams_->get<int>("time int algo")==INPAR::FLUID::timeint_stationary)
    num_sols = 1;
  else
    num_sols = fluidParams_->get<int>("max number timesteps")+1;

  if (num_sols!=fluidvel_->size())
    dserror("fluid field and time step numbers do not fit: n_f = %i, n_t = %i",fluidvel_->size(),num_sols);

  if (num_sols!=adjointvel_->size())
    dserror("adjoint field and time step numbers do not fit: n_a = %i, n_t = %i",adjointvel_->size(),num_sols);

  return true;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::TransformFlowFields()
{
  const Epetra_Map* colmap = fluiddis_->DofColMap();

  // if maps have been mapped before, dont change them
  if (fluidvel_->begin()->second->Map().PointSameAs(*colmap)) return;

  RCP<Epetra_Vector> vec = rcp(new Epetra_Vector(*colmap,false));

  for (map<int,RCP<Epetra_Vector> >::iterator i=fluidvel_->begin();
      i!=fluidvel_->end();i++)
  {
    // export vector from row to column map
    LINALG::Export(*i->second,*vec);

    // set new vector
    i->second = vec;
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=adjointvel_->begin();
      i!=adjointvel_->end();i++)
  {
    // export vector from row to column map
    LINALG::Export(*i->second,*vec);

    // set new vector
    i->second = vec;
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* TOPOPT::Optimizer::RowMap()
{
  return optidis_->NodeRowMap();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* TOPOPT::Optimizer::ColMap()
{
  return optidis_->NodeColMap();
}

