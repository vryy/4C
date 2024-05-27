/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop

\level 0

*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_str_timeada_joint.hpp"

#include "4C_adapter_str_timeloop.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_structure_new_solver_factory.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedataio.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"
#include "4C_structure_new_timint_factory.hpp"

#include <Teuchos_ParameterList.hpp>

#include <tuple>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureTimeAdaJoint::StructureTimeAdaJoint(Teuchos::RCP<Structure> structure)
    : StructureTimeAda(structure), sta_(Teuchos::null), sta_wrapper_(Teuchos::null)
{
  if (stm_->is_setup()) setup_auxiliar();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAdaJoint::setup_auxiliar()
{
  const Teuchos::ParameterList& sdyn = GLOBAL::Problem::Instance()->structural_dynamic_params();
  const Teuchos::ParameterList& jep = sdyn.sublist("TIMEADAPTIVITY").sublist("JOINT EXPLICIT");

  // get the parameters of the auxiliary integrator
  Teuchos::ParameterList adyn(sdyn);
  adyn.remove("TIMEADAPTIVITY");
  for (auto i = jep.begin(); i != jep.end(); ++i) adyn.setEntry(jep.name(i), jep.entry(i));

  // construct the auxiliary time integrator
  sta_ = STR::TIMINT::BuildStrategy(adyn);

  ///// setup dataio
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  //
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  ioflags->set("STDOUTEVRY", 0);
  //
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::ParameterList& nox = xparams->sublist("NOX");
  nox = problem->StructuralNoxParams();
  //
  Teuchos::RCP<IO::DiscretizationWriter> output = stm_->discretization()->Writer();
  //
  Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio = Teuchos::rcp(new STR::TIMINT::BaseDataIO());
  dataio->Init(*ioflags, adyn, *xparams, output);
  dataio->Setup();

  ///// setup datasdyn
  Teuchos::RCP<std::set<enum INPAR::STR::ModelType>> modeltypes =
      Teuchos::rcp(new std::set<enum INPAR::STR::ModelType>());
  modeltypes->insert(INPAR::STR::model_structure);
  //
  Teuchos::RCP<std::set<enum INPAR::STR::EleTech>> eletechs =
      Teuchos::rcp(new std::set<enum INPAR::STR::EleTech>());
  //
  Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<CORE::LINALG::Solver>>>
      linsolvers = STR::SOLVER::BuildLinSolvers(*modeltypes, adyn, *stm_->discretization());
  //
  Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn = STR::TIMINT::BuildDataSDyn(adyn);
  datasdyn->Init(stm_->discretization(), adyn, *xparams, modeltypes, eletechs, linsolvers);
  datasdyn->Setup();

  // setup global state
  Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate =
      STR::TIMINT::build_data_global_state();
  dataglobalstate->Init(stm_->discretization(), adyn, datasdyn);
  dataglobalstate->Setup();

  // setup auxiliary integrator
  sta_->Init(dataio, datasdyn, dataglobalstate);
  sta_->Setup();

  // setup wrapper
  sta_wrapper_ = Teuchos::rcp(new ADAPTER::StructureTimeLoop(sta_));

  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    const STR::TIMINT::Base& sti = *stm_;
    const auto& gstate = sti.data_global_state();
    dataglobalstate->GetDisN()->Update(1.0, *(gstate.GetDisN()), 0.0);
    dataglobalstate->GetVelN()->Update(1.0, *(gstate.GetVelN()), 0.0);
    dataglobalstate->GetAccN()->Update(1.0, *(gstate.GetAccN()), 0.0);
  }

  // check explicitness
  if (sta_->IsImplicit())
  {
    FOUR_C_THROW("Implicit might work, but please check carefully");
  }

  // check order
  if (sta_->method_order_of_accuracy_dis() > stm_->method_order_of_accuracy_dis())
  {
    ada_ = ada_upward;
  }
  else if (sta_->method_order_of_accuracy_dis() < stm_->method_order_of_accuracy_dis())
  {
    ada_ = ada_downward;
  }
  else if (sta_->MethodName() == stm_->MethodName())
  {
    ada_ = ada_ident;
  }
  else
  {
    ada_ = ada_orderequal;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string ADAPTER::StructureTimeAdaJoint::MethodTitle() const
{
  return "JointExplicit_" + sta_->MethodTitle();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimeAdaJoint::method_order_of_accuracy_dis() const
{
  return sta_->method_order_of_accuracy_dis();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimeAdaJoint::method_order_of_accuracy_vel() const
{
  return sta_->method_order_of_accuracy_vel();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::StructureTimeAdaJoint::method_lin_err_coeff_dis() const
{
  return sta_->method_lin_err_coeff_dis();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::StructureTimeAdaJoint::method_lin_err_coeff_vel() const
{
  return sta_->method_lin_err_coeff_vel();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
enum ADAPTER::StructureTimeAda::AdaEnum ADAPTER::StructureTimeAdaJoint::MethodAdaptDis() const
{
  return ada_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAdaJoint::integrate_step_auxiliar()
{
  // set current step size
  sta_->SetDeltaTime(stepsize_);
  sta_->SetTimeNp(time_ + stepsize_);

  // integrate the auxiliary time integrator one step in time
  // buih: another solution is to use the wrapper, but it will do more than necessary
  const STR::TIMINT::Base& sta = *sta_;
  const auto& gstate = sta.data_global_state();

  sta_->IntegrateStep();

  // copy onto target
  locerrdisn_->Update(1.0, *(gstate.GetDisNp()), 0.0);

  // reset
  sta_->reset_step();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAdaJoint::update_auxiliar()
{
  // copy the data from main integrator to the auxiliary one
  // for reference: the vector map of the global state vectors may need to be checked to ensure they
  // are the same
  const STR::TIMINT::Base& stm = *stm_;
  const STR::TIMINT::BaseDataGlobalState& gstate_i = stm.data_global_state();

  const STR::TIMINT::Base& sta = *sta_;
  const STR::TIMINT::BaseDataGlobalState& gstate_a_const = sta.data_global_state();
  STR::TIMINT::BaseDataGlobalState& gstate_a =
      const_cast<STR::TIMINT::BaseDataGlobalState&>(gstate_a_const);

  gstate_a.GetDisNp()->Update(1.0, (*gstate_i.GetDisN()), 0.0);
  gstate_a.GetVelNp()->Update(1.0, (*gstate_i.GetVelN()), 0.0);
  gstate_a.GetAccNp()->Update(1.0, (*gstate_i.GetAccN()), 0.0);
  gstate_a.GetMultiDis()->UpdateSteps((*gstate_i.GetDisN()));
  gstate_a.GetMultiVel()->UpdateSteps((*gstate_i.GetVelN()));
  gstate_a.GetMultiAcc()->UpdateSteps((*gstate_i.GetAccN()));

  gstate_a.GetTimeNp() = gstate_i.GetTimeNp();
  gstate_a.GetDeltaTime()->UpdateSteps((*gstate_i.GetDeltaTime())[0]);
  gstate_a.GetFviscoNp()->Update(1.0, (*gstate_i.GetFviscoN()), 0.0);
  gstate_a.GetFviscoN()->Update(1.0, (*gstate_i.GetFviscoN()), 0.0);
  gstate_a.GetFinertialNp()->Update(1.0, (*gstate_i.GetFinertialN()), 0.0);
  gstate_a.GetFinertialN()->Update(1.0, (*gstate_i.GetFinertialN()), 0.0);
  gstate_a.GetFintNp()->Update(1.0, (*gstate_i.GetFintN()), 0.0);
  gstate_a.GetFintN()->Update(1.0, (*gstate_i.GetFintN()), 0.0);
  gstate_a.GetFextNp()->Update(1.0, (*gstate_i.GetFextN()), 0.0);
  gstate_a.GetFextN()->Update(1.0, (*gstate_i.GetFextN()), 0.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAdaJoint::reset_step()
{
  // base reset
  ADAPTER::StructureTimeAda::reset_step();
  // set current step size
  sta_->SetDeltaTime(stepsize_);
  sta_->SetTimeNp(time_ + stepsize_);
  // reset the integrator
  sta_->reset_step();
}

FOUR_C_NAMESPACE_CLOSE
