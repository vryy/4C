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
#include "../drt_io/io.H"
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
    const Teuchos::ParameterList& params,
    Teuchos::RCP<IO::DiscretizationWriter>& output
) :
optidis_(optidis),
fluiddis_(fluiddis),
params_(params),
gradienttype_(DRT::INPUT::IntegralValue<INPAR::TOPOPT::GradientType>(params_,"GRADIENT_TYPE"))
{
  const Teuchos::ParameterList& optimizer_params = params_.sublist("TOPOLOGY OPTIMIZER");

  // fluid fields for optimization
  fluidvel_ = Teuchos::rcp(new std::map<int,Teuchos::RCP<Epetra_Vector> >);
  adjointvel_ = Teuchos::rcp(new std::map<int,Teuchos::RCP<Epetra_Vector> >);

  // topology density fields
  dens_ = Teuchos::rcp(new Epetra_Vector(*optidis_->NodeColMap(),false));

  // value of the objective function
  obj_ = 0.0;
  // gradient of the objective function
  obj_der_ = Teuchos::rcp(new Epetra_Vector(*optidis_->NodeRowMap()));

  /// number of constraints
  switch (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(optimizer_params,"TESTCASE"))
  {
  case INPAR::TOPOPT::optitest_no:
  case INPAR::TOPOPT::optitest_workflow_without_fluiddata:
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
  default:
  {
    dserror("unknown optimization case");
    break;
  }
  }

  // value of the constraint(s);
  constr_ = Teuchos::rcp(new Epetra_SerialDenseVector(num_constr_));
  // gradient of the constraint(s)
  constr_der_ = Teuchos::rcp(new Epetra_MultiVector(*optidis_->NodeRowMap(),num_constr_,false));

  // set initial density field if present
  SetInitialDensityField(DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialDensityField>(optimizer_params,"INITIALFIELD"),
      optimizer_params.get<int>("INITFUNCNO"));

  Teuchos::RCP<Epetra_Vector> dens = Teuchos::rcp(new Epetra_Vector(*optidis_->NodeRowMap()));
  LINALG::Export(*dens_,*dens);

  if (1==1) // use this when different optimizers are present
  {
    optimizer_ = Teuchos::rcp(new OPTI::GCMMA(
        optidis_,
        params_,
        dens,
        num_constr_,
        Teuchos::null,
        Teuchos::null,
        output
    ));
  }

  // write output using the derived optimizer -> optimizer must be present here
  Output();

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeValues(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector>& constraints
)
{
  // initialize local values processor-wise
  double loc_obj = 0.0;
  Epetra_SerialDenseVector loc_constr = Epetra_SerialDenseVector(num_constr_);

  const bool doAdjoint = false;

  // check if all data is present
  CheckData(doAdjoint);

  // Transform fluid field to colmap
  TransformFlowFields(doAdjoint,true);

  Teuchos::ParameterList params;

  params.set("action","compute_values");

  params.set("objective_value",&loc_obj);
  params.set("constraint_values",&loc_constr);

  params.set("fluidvel",fluidvel_);
  params.set("fluiddis",fluiddis_);

  optidis_->ClearState();

  optidis_->SetState("density",dens_);
  optidis_->Evaluate(params,Teuchos::null,Teuchos::null);

  optidis_->ClearState();

  // communicate local values over processors
  optidis_->Comm().SumAll(&loc_obj,&objective,1);
  optidis_->Comm().SumAll(loc_constr.Values(),constraints->Values(),num_constr_);

  // Re-Transform fluid field to rowmap
  TransformFlowFields(doAdjoint,false); // TODO is this required (back-mapping + 2nd argument)
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeGradients(
    Teuchos::RCP<Epetra_Vector> obj_der,
    Teuchos::RCP<Epetra_MultiVector> constr_der
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" evaluate objective gradient");
  obj_der->PutScalar(0.0);
  constr_der->PutScalar(0.0);

  bool doAdjoint = true;
  if (gradienttype_==INPAR::TOPOPT::gradientByAdjoints) doAdjoint = true;
  else if (gradienttype_==INPAR::TOPOPT::gradientByFD1) doAdjoint = false;
  else dserror("unknown type of gradient computation");

  // check if all data is present
  CheckData(doAdjoint);

  // Transform fluid field to colmap
  TransformFlowFields(doAdjoint,true);

  Teuchos::ParameterList params;

  params.set("action","compute_gradients");

  params.set("constraints_derivations",constr_der);

  params.set("fluidvel",fluidvel_);
  params.set("adjointvel",adjointvel_);
  params.set("fluiddis",fluiddis_);

//  for (std::map<int,Teuchos::RCP<Epetra_Vector> >::iterator i=fluidvel_->begin();
//      i!=fluidvel_->end();i++)
//    cout << "in optimizer at gradients fluidvel of step " << i->first << " is " << *i->second << endl;
//  for (std::map<int,Teuchos::RCP<Epetra_Vector> >::iterator i=adjointvel_->begin();
//      i!=adjointvel_->end();i++)
//    cout << "in optimizer at gradients adjointvel of step " << i->first << " is " << *i->second << endl;

  optidis_->ClearState();

  optidis_->SetState("density",dens_);

  optidis_->Evaluate(params,Teuchos::null,obj_der);

  optidis_->ClearState();

  // Re-Transform fluid field to rowmap
  TransformFlowFields(doAdjoint,false); // TODO is this required (back-mapping + 2nd argument)
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ComputeGradientDirectionForFD(const double value, const int GID)
{
  double objective = 0.0;
  Teuchos::RCP<Epetra_SerialDenseVector> constraints = Teuchos::rcp(new Epetra_SerialDenseVector(num_constr_));

  ComputeValues(objective,constraints);

  if (dens_->Map().MyGID(GID))
  {
    obj_der_->ReplaceGlobalValue(GID,0,(objective-obj_)/value);

    for (int i=0;i<constraints->Length();i++)
      constr_der_->ReplaceGlobalValue(GID,i,((*constraints)(i)-(*constr_)(i))/value);
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::Output()
{
  optimizer_->Output();
  optimizer_->OutputWriter()->WriteVector("obj_der",obj_der_);
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
    for(int lnodeid=0;lnodeid<optidis_->NumMyColNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = optidis_->lColNode(lnodeid);

      // evaluate component k of spatial function
      double initialval = DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(0,lnode->X(),0.0,NULL); // scalar
      int err = dens_->ReplaceMyValues(1,&initialval,&lnodeid); // lnodeid = ldofid
      if (err != 0) dserror("dof not on proc");
    }
    break;
  }
  default:
  {
    dserror("unknown initial field");
    break;
  }
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
        {
          dserror("unknown material %s",imat->Name().c_str());
          break;
        }
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

  opti_ele_params.set<INPAR::TOPOPT::OptiCase>("opti_case",DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_.sublist("TOPOLOGY OPTIMIZER"),"TESTCASE"));

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
    RCP<Epetra_Vector> new_vel = Teuchos::rcp(new Epetra_Vector(*vel)); // copy
    fluidvel_->insert(std::pair<int,RCP<Epetra_Vector> >(step,new_vel));
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
    RCP<Epetra_Vector> new_vel = Teuchos::rcp(new Epetra_Vector(*vel));
    adjointvel_->insert(std::pair<int,RCP<Epetra_Vector> >(step,new_vel));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::Iterate(
    bool& doGradient
)
{
  Teuchos::RCP<Epetra_Vector> dens = Teuchos::null;

  if (doGradient)
  {
    dens = optimizer_->Iterate(
        obj_,
        obj_der_,
        constr_,
        constr_der_
    );
  }
  else
  {
    dens = optimizer_->Iterate(
        obj_,
        Teuchos::null,
        constr_,
        Teuchos::null
    );
  }

  LINALG::Export(*dens,*dens_);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::FinishIteration(
    bool& doGradient
)
{
  return optimizer_->FinishIteration(
      obj_,
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
        obj_,
        obj_der_,
        constr_,
        constr_der_
    );
  }
  else
  {
    converged = optimizer_->Converged(
        obj_,
        Teuchos::null,
        constr_,
        Teuchos::null
    );
  }

  return converged;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool TOPOPT::Optimizer::CheckData(const bool doAdjoint)
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

  if ((doAdjoint==true) and (num_sols!=adjointvel_->size()))
    dserror("adjoint field and time step numbers do not fit: n_a = %i, n_t = %i",adjointvel_->size(),num_sols);

  return true;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::TransformFlowFields(const bool doAdjoint, const bool rowToCol)
{
  const Epetra_Map* targetmap = NULL;
  if (rowToCol) targetmap= fluiddis_->DofColMap();
  else targetmap = fluiddis_->DofRowMap();

  for (std::map<int,Teuchos::RCP<Epetra_Vector> >::iterator i=fluidvel_->begin();
      i!=fluidvel_->end();i++)
  {
    // if maps have been mapped before, dont change them
    if (i->second->Map().PointSameAs(*targetmap)) continue;

    RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*targetmap,false));

    // export vector to target map
    LINALG::Export(*i->second,*vec);
    i->second = vec;
  }

  // if maps have been mapped before, dont change them
  if (doAdjoint==true)
  {
    for (std::map<int,Teuchos::RCP<Epetra_Vector> >::iterator i=adjointvel_->begin();
        i!=adjointvel_->end();i++)
    {
      if (i->second->Map().PointSameAs(*targetmap)) continue;

      RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*targetmap,false));

      // export vector to target map
      LINALG::Export(*i->second,*vec);
      i->second = vec;
    }
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
void TOPOPT::Optimizer::AdoptDensityForFD(const double value, const int GID)
{
  if (dens_->Map().MyGID(GID))
    dens_->SumIntoGlobalValue(GID,0,value);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const int TOPOPT::Optimizer::Iter() const
{
  return optimizer_->Iter();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::ReadRestart(const int step)
{
  // required here for objective derivation
  Teuchos::RCP<IO::DiscretizationReader> reader = optimizer_->ReadRestart(step);

  Teuchos::RCP<Epetra_Vector> dens = optimizer_->X();
  LINALG::Export(*dens,*dens_);
  reader->ReadVector(obj_der_,"obj_der");

  return;
}





