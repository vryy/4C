/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.cpp

\brief optimizer of the topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"
#include "opti_GCMMA.H"
#include "opti_resulttest.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/optimization_density.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
TOPOPT::Optimizer::Optimizer(
    RCP<DRT::Discretization> optidis,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    const ParameterList& params,
    Teuchos::RCP<IO::DiscretizationWriter>& output
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
  dens_ = rcp(new Epetra_Vector(*optidis_->NodeRowMap(),false));

  // value of the objective function
  obj_value_ = 0.0;
  // gradient of the objective function
  obj_grad_ = rcp(new Epetra_Vector(*optidis_->NodeRowMap()));

  /// number of constraints
  switch (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(optimizer_params,"TESTCASE"))
  {
  case INPAR::TOPOPT::optitest_no:
  case INPAR::TOPOPT::optitest_snake_one_constr:
  {
    num_constr_ = 1;
    break;
  }
  case INPAR::TOPOPT::optitest_snake_multiple_constr:
  {
    num_constr_ = 5;
    break;
  }
  }

  // value of the constraint(s);
  constr_ = rcp(new Epetra_SerialDenseVector(num_constr_));
  // gradient of the constraint(s)
  constr_deriv_ = rcp(new Epetra_MultiVector(*optidis_->NodeRowMap(),num_constr_,false));

  // set initial density field if present
  SetInitialDensityField(DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialDensityField>(optimizer_params,"INITIALFIELD"),
      optimizer_params.get<int>("INITFUNCNO"));

  optimizer_ = rcp(new OPTI::GCMMA(
      optidis_,
      optimizer_params,
      dens_,
      num_constr_,
      Teuchos::null,
      Teuchos::null,
      output
  ));
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeValues()
{
  obj_value_ = 0.0; // initialize with zero

  // initialize
  double* constr = constr_->Values();
  for (int i=0;i<num_constr_;i++)
  {
    *constr = 0.0;
    constr++;
  }

  // check if all data is present
  DataComplete(false);

  // Transform fluid and adjoint field so that it is better readable at element level
  TransformFlowFields();

  Teuchos::ParameterList params;

  params.set("action","compute_values");

  params.set("objective_value",obj_value_);
  params.set("constraint_values",constr_);

  params.set("fluidvel",fluidvel_);
  params.set("fluiddis",fluiddis_);

  optidis_->ClearState();

  optidis_->SetState("density",dens_);
  optidis_->Evaluate(params,Teuchos::null,Teuchos::null);

  optidis_->ClearState();

  // extract objective value from parameter list
  obj_value_ = params.get<double>("objective_value");
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeGradients()
{
  TEUCHOS_FUNC_TIME_MONITOR(" evaluate objective gradient");
  obj_grad_->PutScalar(0.0);
  constr_deriv_->PutScalar(0.0);

  DataComplete(true);

  // Transform fluid and adjoint field so that it is better readable at element level
  TransformFlowFields();

  Teuchos::ParameterList params;

  params.set("action","compute_gradients");

  params.set("constraints_derivations",constr_deriv_);

  params.set("fluidvel",fluidvel_);
  params.set("adjointvel",adjointvel_);
  params.set("fluiddis",fluiddis_);

  optidis_->ClearState();

  optidis_->SetState("density",dens_);
  optidis_->Evaluate(params,Teuchos::null,obj_grad_);

  optidis_->ClearState();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::Output()
{
  optimizer_->Output();
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
    dens_->PutScalar(0.0);
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
      int err = dens_->ReplaceMyValues(1,&initialval,&lnodeid); // lnodeid = ldofid
      if (err != 0) dserror("dof not on proc");
    }
    break;
  }
  default:
    dserror("unknown initial field");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ImportFlowParams(
    Teuchos::RCP<Teuchos::ParameterList>& fluidParams
)
{
  // save the fluid parameter
  fluidParams_ = fluidParams;

  // set the general parameter on element level
  Teuchos::ParameterList opti_ele_params;
  opti_ele_params.set("action","set_general_optimization_parameter");

  // set material parameters of fluid and optimization material
  {
    // get the material
    const int nummat = DRT::Problem::Instance()->Materials()->Num();
    for (int id = 1; id-1 < nummat; ++id)
    {
      Teuchos::RCP<const MAT::PAR::Material> imat = DRT::Problem::Instance()->Materials()->ById(id);

      if (imat == Teuchos::null)
        dserror("Could not find material Id %d", id);
      else
      {
        switch (imat->Type())
        {
        case INPAR::MAT::m_fluid:
        {
          const MAT::PAR::Parameter* matparam = imat->Parameter();
          const MAT::PAR::NewtonianFluid* mat = static_cast<const MAT::PAR::NewtonianFluid* >(matparam);

          opti_ele_params.set("density",mat->density_);
          opti_ele_params.set("viscosity",mat->viscosity_);
          break;
        }
        case INPAR::MAT::m_opti_dens:
        {
          const MAT::PAR::Parameter* matparam = imat->Parameter();
          const MAT::PAR::TopOptDens* mat = static_cast<const MAT::PAR::TopOptDens* >(matparam);

          opti_ele_params.set("min_poro",mat->poro_bd_down_);
          opti_ele_params.set("max_poro",mat->poro_bd_up_);
          opti_ele_params.set("smear_fac",mat->smear_fac_);
          break;
        }
        default:
          dserror("unknown material %s",imat->Name().c_str());
        }
      }
    }
  }

  // flow parameter
  opti_ele_params.set("timealgo",fluidParams_->get<int>("time int algo"));
  opti_ele_params.set("dt",fluidParams_->get<double>("time step size"));
  opti_ele_params.set("maxtimesteps",fluidParams_->get<int>("max number timesteps"));
  opti_ele_params.set("theta",fluidParams_->get<double>("theta"));
  opti_ele_params.set("theta_pre",fluidParams_->get<double>("theta_pre"));
  opti_ele_params.set("theta_div",fluidParams_->get<double>("theta_div"));

  // objective parameter
  opti_ele_params.set("dissipation",fluidParams_->get<bool>("OBJECTIVE_DISSIPATION"));
  opti_ele_params.set("pres_drop",fluidParams_->get<bool>("OBJECTIVE_PRESSURE_DROP"));

  opti_ele_params.set("dissipation_fac",fluidParams_->get<double>("DISSIPATION_FAC"));
  opti_ele_params.set("pres_drop_fac",fluidParams_->get<double>("PRESSURE_DROP_FAC"));

  opti_ele_params.set("vol_bd",params_.sublist("TOPOLOGY OPTIMIZER").get<double>("VOLUME_BOUNDARY"));

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
void TOPOPT::Optimizer::Iterate(
    bool& doGradient
)
{
  if (doGradient)
  {
    dens_ = optimizer_->Iterate(
        obj_value_,
        obj_grad_,
        constr_,
        constr_deriv_
    );
  }
  else
  {
    dens_ = optimizer_->Iterate(
        obj_value_,
        Teuchos::null,
        constr_,
        Teuchos::null
    );
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::FinishIteration(
    bool& doGradient
)
{
  return optimizer_->FinishIteration(
      obj_value_,
      constr_,
      doGradient
  );
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::Converged(
    bool& doGradient
)
{
  bool converged = false;

  // dont use gradients if not needed (safety)
  if (doGradient)
  {
    converged = optimizer_->Converged(
        obj_value_,
        obj_grad_,
        constr_,
        constr_deriv_
    );
  }
  else
  {
    converged = optimizer_->Converged(
        obj_value_,
        Teuchos::null,
        constr_,
        Teuchos::null
    );
  }

  return converged;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::DataComplete(bool checkAdjoint) const
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

  if ((num_sols!=adjointvel_->size()) and (checkAdjoint))
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> TOPOPT::Optimizer::CreateFieldTest()
{
  return Teuchos::rcp(new OPTI::OptiResultTest(*optimizer_));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const double TOPOPT::Optimizer::Iter() const {return optimizer_->Iter();}


