/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop

\level 0

*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_str_timeada_joint.hpp"

#include "4C_adapter_str_timeloop.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
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
Adapter::StructureTimeAdaJoint::StructureTimeAdaJoint(Teuchos::RCP<Structure> structure)
    : StructureTimeAda(structure), sta_(Teuchos::null), sta_wrapper_(Teuchos::null)
{
  if (stm_->is_setup()) setup_auxiliar();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaJoint::setup_auxiliar()
{
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();
  const Teuchos::ParameterList& jep = sdyn.sublist("TIMEADAPTIVITY").sublist("JOINT EXPLICIT");

  // get the parameters of the auxiliary integrator
  Teuchos::ParameterList adyn(sdyn);
  adyn.remove("TIMEADAPTIVITY");
  for (auto i = jep.begin(); i != jep.end(); ++i) adyn.setEntry(jep.name(i), jep.entry(i));

  // construct the auxiliary time integrator
  sta_ = STR::TimeInt::build_strategy(adyn);

  ///// setup dataio
  Global::Problem* problem = Global::Problem::Instance();
  //
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  ioflags->set("STDOUTEVRY", 0);
  //
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::ParameterList& nox = xparams->sublist("NOX");
  nox = problem->StructuralNoxParams();
  //
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = stm_->discretization()->Writer();
  //
  Teuchos::RCP<STR::TimeInt::BaseDataIO> dataio = Teuchos::rcp(new STR::TimeInt::BaseDataIO());
  dataio->init(*ioflags, adyn, *xparams, output);
  dataio->setup();

  ///// setup datasdyn
  Teuchos::RCP<std::set<enum Inpar::STR::ModelType>> modeltypes =
      Teuchos::rcp(new std::set<enum Inpar::STR::ModelType>());
  modeltypes->insert(Inpar::STR::model_structure);
  //
  Teuchos::RCP<std::set<enum Inpar::STR::EleTech>> eletechs =
      Teuchos::rcp(new std::set<enum Inpar::STR::EleTech>());
  //
  Teuchos::RCP<std::map<enum Inpar::STR::ModelType, Teuchos::RCP<Core::LinAlg::Solver>>>
      linsolvers = STR::SOLVER::build_lin_solvers(*modeltypes, adyn, *stm_->discretization());
  //
  Teuchos::RCP<STR::TimeInt::BaseDataSDyn> datasdyn = STR::TimeInt::build_data_sdyn(adyn);
  datasdyn->init(stm_->discretization(), adyn, *xparams, modeltypes, eletechs, linsolvers);
  datasdyn->setup();

  // setup global state
  Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> dataglobalstate =
      STR::TimeInt::build_data_global_state();
  dataglobalstate->init(stm_->discretization(), adyn, datasdyn);
  dataglobalstate->setup();

  // setup auxiliary integrator
  sta_->init(dataio, datasdyn, dataglobalstate);
  sta_->setup();

  // setup wrapper
  sta_wrapper_ = Teuchos::rcp(new Adapter::StructureTimeLoop(sta_));

  const int restart = Global::Problem::Instance()->restart();
  if (restart)
  {
    const STR::TimeInt::Base& sti = *stm_;
    const auto& gstate = sti.data_global_state();
    dataglobalstate->get_dis_n()->Update(1.0, *(gstate.get_dis_n()), 0.0);
    dataglobalstate->get_vel_n()->Update(1.0, *(gstate.get_vel_n()), 0.0);
    dataglobalstate->get_acc_n()->Update(1.0, *(gstate.get_acc_n()), 0.0);
  }

  // check explicitness
  if (sta_->is_implicit())
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
  else if (sta_->method_name() == stm_->method_name())
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
std::string Adapter::StructureTimeAdaJoint::MethodTitle() const
{
  return "JointExplicit_" + sta_->method_title();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::StructureTimeAdaJoint::method_order_of_accuracy_dis() const
{
  return sta_->method_order_of_accuracy_dis();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::StructureTimeAdaJoint::method_order_of_accuracy_vel() const
{
  return sta_->method_order_of_accuracy_vel();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::StructureTimeAdaJoint::method_lin_err_coeff_dis() const
{
  return sta_->method_lin_err_coeff_dis();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::StructureTimeAdaJoint::method_lin_err_coeff_vel() const
{
  return sta_->method_lin_err_coeff_vel();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
enum Adapter::StructureTimeAda::AdaEnum Adapter::StructureTimeAdaJoint::MethodAdaptDis() const
{
  return ada_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaJoint::integrate_step_auxiliar()
{
  // set current step size
  sta_->SetDeltaTime(stepsize_);
  sta_->SetTimeNp(time_ + stepsize_);

  // integrate the auxiliary time integrator one step in time
  // buih: another solution is to use the wrapper, but it will do more than necessary
  const STR::TimeInt::Base& sta = *sta_;
  const auto& gstate = sta.data_global_state();

  sta_->IntegrateStep();

  // copy onto target
  locerrdisn_->Update(1.0, *(gstate.get_dis_np()), 0.0);

  // reset
  sta_->reset_step();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaJoint::update_auxiliar()
{
  // copy the data from main integrator to the auxiliary one
  // for reference: the vector map of the global state vectors may need to be checked to ensure they
  // are the same
  const STR::TimeInt::Base& stm = *stm_;
  const STR::TimeInt::BaseDataGlobalState& gstate_i = stm.data_global_state();

  const STR::TimeInt::Base& sta = *sta_;
  const STR::TimeInt::BaseDataGlobalState& gstate_a_const = sta.data_global_state();
  STR::TimeInt::BaseDataGlobalState& gstate_a =
      const_cast<STR::TimeInt::BaseDataGlobalState&>(gstate_a_const);

  gstate_a.get_dis_np()->Update(1.0, (*gstate_i.get_dis_n()), 0.0);
  gstate_a.get_vel_np()->Update(1.0, (*gstate_i.get_vel_n()), 0.0);
  gstate_a.get_acc_np()->Update(1.0, (*gstate_i.get_acc_n()), 0.0);
  gstate_a.get_multi_dis()->UpdateSteps((*gstate_i.get_dis_n()));
  gstate_a.get_multi_vel()->UpdateSteps((*gstate_i.get_vel_n()));
  gstate_a.get_multi_acc()->UpdateSteps((*gstate_i.get_acc_n()));

  gstate_a.get_time_np() = gstate_i.get_time_np();
  gstate_a.get_delta_time()->UpdateSteps((*gstate_i.get_delta_time())[0]);
  gstate_a.get_fvisco_np()->Update(1.0, (*gstate_i.get_fvisco_n()), 0.0);
  gstate_a.get_fvisco_n()->Update(1.0, (*gstate_i.get_fvisco_n()), 0.0);
  gstate_a.get_finertial_np()->Update(1.0, (*gstate_i.get_finertial_n()), 0.0);
  gstate_a.get_finertial_n()->Update(1.0, (*gstate_i.get_finertial_n()), 0.0);
  gstate_a.get_fint_np()->Update(1.0, (*gstate_i.get_fint_n()), 0.0);
  gstate_a.get_fint_n()->Update(1.0, (*gstate_i.get_fint_n()), 0.0);
  gstate_a.get_fext_np()->Update(1.0, (*gstate_i.get_fext_n()), 0.0);
  gstate_a.get_fext_n()->Update(1.0, (*gstate_i.get_fext_n()), 0.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaJoint::reset_step()
{
  // base reset
  Adapter::StructureTimeAda::reset_step();
  // set current step size
  sta_->SetDeltaTime(stepsize_);
  sta_->SetTimeNp(time_ + stepsize_);
  // reset the integrator
  sta_->reset_step();
}

FOUR_C_NAMESPACE_CLOSE
