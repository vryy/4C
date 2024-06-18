/*----------------------------------------------------------------------*/
/*! \file

\brief is the base for the different types of mesh and level-set based coupling conditions and
thereby builds the bridge between the xfluid class and the cut-library

\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_xfem_coupling_base.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fluid_ele_parameter_xfem.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_utils_function.hpp"
#include "4C_xfem_interface_utils.hpp"
#include "4C_xfem_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

Inpar::XFEM::EleCouplingCondType XFEM::CondType_stringToEnum(const std::string& condname)
{
  if (condname == "XFEMSurfFSIPart")
    return Inpar::XFEM::CouplingCond_SURF_FSI_PART;
  else if (condname == "XFEMSurfFSIMono")
    return Inpar::XFEM::CouplingCond_SURF_FSI_MONO;
  else if (condname == "XFEMSurfFPIMono" || condname == "XFEMSurfFPIMono_ps_ps" ||
           condname == "XFEMSurfFPIMono_ps_pf" || condname == "XFEMSurfFPIMono_pf_ps" ||
           condname == "XFEMSurfFPIMono_pf_pf")
    return Inpar::XFEM::CouplingCond_SURF_FPI_MONO;
  else if (condname == "XFEMSurfFluidFluid")
    return Inpar::XFEM::CouplingCond_SURF_FLUIDFLUID;
  else if (condname == "XFEMLevelsetWeakDirichlet")
    return Inpar::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET;
  else if (condname == "XFEMLevelsetNeumann")
    return Inpar::XFEM::CouplingCond_LEVELSET_NEUMANN;
  else if (condname == "XFEMLevelsetNavierSlip")
    return Inpar::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP;
  else if (condname == "XFEMLevelsetTwophase")
    return Inpar::XFEM::CouplingCond_LEVELSET_TWOPHASE;
  else if (condname == "XFEMLevelsetCombustion")
    return Inpar::XFEM::CouplingCond_LEVELSET_COMBUSTION;
  else if (condname == "XFEMSurfWeakDirichlet")
    return Inpar::XFEM::CouplingCond_SURF_WEAK_DIRICHLET;
  else if (condname == "XFEMSurfNeumann")
    return Inpar::XFEM::CouplingCond_SURF_NEUMANN;
  else if (condname == "XFEMSurfNavierSlip")
    return Inpar::XFEM::CouplingCond_SURF_NAVIER_SLIP;
  else if (condname == "XFEMSurfNavierSlipTwoPhase")
    return Inpar::XFEM::CouplingCond_SURF_NAVIER_SLIP_TWOPHASE;
  // else FOUR_C_THROW("condition type not supported: %s", condname.c_str());

  return Inpar::XFEM::CouplingCond_NONE;
}

/*--------------------------------------------------------------------------*
 * constructor
 *--------------------------------------------------------------------------*/
XFEM::CouplingBase::CouplingBase(
    Teuchos::RCP<Core::FE::Discretization>& bg_dis,  ///< background discretization
    const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                   ///< discretization is identified
    Teuchos::RCP<Core::FE::Discretization>&
        cond_dis,           ///< full discretization from which the cutter discretization is derived
    const int coupling_id,  ///< id of composite of coupling conditions
    const double time,      ///< time
    const int step          ///< time step
    )
    : nsd_(Global::Problem::Instance()->NDim()),
      bg_dis_(bg_dis),
      cond_name_(cond_name),
      cond_dis_(cond_dis),
      coupling_id_(coupling_id),
      cutter_dis_(Teuchos::null),
      coupl_dis_(Teuchos::null),
      coupl_name_(""),
      averaging_strategy_(Inpar::XFEM::invalid),
      myrank_(bg_dis_->Comm().MyPID()),
      dt_(-1.0),
      time_(time),
      step_(step),
      issetup_(false),
      isinit_(false)
{
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::CouplingBase::init()
{
  // TODO: correct handling of init and setup flags for derived classes

  // ---------------------------------------------------------------------------
  // We need to call setup() after init()
  // ---------------------------------------------------------------------------
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // do Init
  // ---------------------------------------------------------------------------

  if (dofset_coupling_map_.empty()) FOUR_C_THROW("Call set_dof_set_coupling_map() first!");

  SetCouplingDofsets();

  // set the name of the coupling object to allow access from outside via the name
  set_coupling_name();

  // set list of conditions that will be copied to the new cutter discretization
  set_conditions_to_copy();

  // create a cutter discretization from conditioned nodes of the given coupling discretization or
  // simply clone the discretization
  set_cutter_discretization();

  // set unique element conditions
  set_element_conditions();

  // set condition specific parameters
  set_condition_specific_parameters();

  // set the averaging strategy
  set_averaging_strategy();

  // set coupling discretization
  set_coupling_discretization();

  // initialize element level configuration map (no evaluation)
  init_configuration_map();

  // ---------------------------------------------------------------------------
  // set isInit flag
  // ---------------------------------------------------------------------------
  isinit_ = true;

  // good bye
  return;
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::CouplingBase::setup()
{
  check_init();

  // ---------------------------------------------------------------------------
  // do setup
  // ---------------------------------------------------------------------------

  // initialize state vectors according to cutter discretization
  init_state_vectors();

  // prepare the output writer for the cutter discretization
  prepare_cutter_output();

  // do condition specific setup
  do_condition_specific_setup();

  // initialize the configuration map
  setup_configuration_map();

  // ---------------------------------------------------------------------------
  // set isSetup flag
  // ---------------------------------------------------------------------------

  issetup_ = true;
}


/*--------------------------------------------------------------------------*
 * Initialize Configuration Map --> No Terms are evaluated at the interface
 *--------------------------------------------------------------------------*/
void XFEM::CouplingBase::init_configuration_map()
{
  // Configuration of Consistency Terms
  // all components:
  configuration_map_[Inpar::XFEM::F_Con_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Con_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Con_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XF_Con_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XS_Con_Col] = std::pair<bool, double>(false, 0.0);
  // normal terms:
  configuration_map_[Inpar::XFEM::F_Con_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Con_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Con_n_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XF_Con_n_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XS_Con_n_Col] = std::pair<bool, double>(false, 0.0);
  // tangential terms:
  configuration_map_[Inpar::XFEM::F_Con_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Con_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Con_t_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XF_Con_t_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XS_Con_t_Col] = std::pair<bool, double>(false, 0.0);

  // Configuration of Adjoint Consistency Terms
  // all components:
  configuration_map_[Inpar::XFEM::F_Adj_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XF_Adj_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XS_Adj_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Adj_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Adj_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::FStr_Adj_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XStr_Adj_Col] = std::pair<bool, double>(false, 0.0);
  // normal terms:
  configuration_map_[Inpar::XFEM::F_Adj_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XF_Adj_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XS_Adj_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Adj_n_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Adj_n_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::FStr_Adj_n_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XStr_Adj_n_Col] = std::pair<bool, double>(false, 0.0);
  // tangential terms:
  configuration_map_[Inpar::XFEM::F_Adj_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XF_Adj_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XS_Adj_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Adj_t_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Adj_t_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::FStr_Adj_t_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::XStr_Adj_t_Col] = std::pair<bool, double>(false, 0.0);

  // Configuration of Penalty Terms
  // all components:
  configuration_map_[Inpar::XFEM::F_Pen_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Pen_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Pen_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Pen_Col] = std::pair<bool, double>(false, 0.0);
  // linearization of penalty terms: at the moment exclusively used for inflow stab
  configuration_map_[Inpar::XFEM::F_Pen_Row_linF1] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Pen_Row_linF2] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Pen_Row_linF3] = std::pair<bool, double>(false, 0.0);
  // normal terms:
  configuration_map_[Inpar::XFEM::F_Pen_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Pen_n_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Pen_n_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Pen_n_Col] = std::pair<bool, double>(false, 0.0);
  // tangential terms:
  configuration_map_[Inpar::XFEM::F_Pen_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Pen_t_Row] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_Pen_t_Col] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_Pen_t_Col] = std::pair<bool, double>(false, 0.0);

  // Starting from here are some special Terms
  configuration_map_[Inpar::XFEM::F_LB_Rhs] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_LB_Rhs] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::F_TJ_Rhs] = std::pair<bool, double>(false, 0.0);
  configuration_map_[Inpar::XFEM::X_TJ_Rhs] = std::pair<bool, double>(false, 0.0);
  return;
}


void XFEM::CouplingBase::set_element_conditions()
{
  // number of column cutter boundary elements
  int nummycolele = cutter_dis_->NumMyColElements();

  cutterele_conds_.clear();
  cutterele_conds_.reserve(nummycolele);

  // initialize the vector invalid coupling-condition type "NONE"
  EleCoupCond init_pair = EleCoupCond(Inpar::XFEM::CouplingCond_NONE, nullptr);
  for (int lid = 0; lid < nummycolele; lid++) cutterele_conds_.push_back(init_pair);

  //-----------------------------------------------------------------------------------
  // loop all column cutting elements on this processor
  for (int lid = 0; lid < nummycolele; lid++)
  {
    Core::Elements::Element* cutele = cutter_dis_->lColElement(lid);

    // loop all possible XFEM-coupling conditions
    for (size_t cond = 0; cond < conditions_to_copy_.size(); cond++)
    {
      Inpar::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(conditions_to_copy_[cond]);

      // non-coupling condition found (e.g. FSI coupling)
      if (cond_type == Inpar::XFEM::CouplingCond_NONE) continue;

      // get all conditions with given condition name
      std::vector<Core::Conditions::Condition*> mycond;
      Core::Conditions::FindElementConditions(cutele, conditions_to_copy_[cond], mycond);

      std::vector<Core::Conditions::Condition*> mynewcond;
      get_condition_by_coupling_id(mycond, coupling_id_, mynewcond);

      Core::Conditions::Condition* cond_unique = nullptr;

      // safety checks
      if (mynewcond.size() == 0)
      {
        continue;  // try the next condition type
      }
      else if (mynewcond.size() == 1)  // unique condition found
      {
        cond_unique = mynewcond[0];
      }
      else if (mynewcond.size() > 1)
      {
        // get the right condition
        FOUR_C_THROW(
            "%i conditions of the same name with coupling id %i, for element %i! %s "
            "coupling-condition not unique!",
            mynewcond.size(), coupling_id_, cutele->Id(), conditions_to_copy_[cond].c_str());
      }

      // non-unique conditions for one cutter element
      if (cutterele_conds_[lid].first != Inpar::XFEM::CouplingCond_NONE)
      {
        FOUR_C_THROW(
            "There are two different condition types for the same cutter dis element with id %i: "
            "1st %i, 2nd %i. Make the XFEM coupling conditions unique!",
            cutele->Id(), cutterele_conds_[lid].first, cond_type);
      }

      // store the unique condition pointer to the cutting element
      cutterele_conds_[lid] = EleCoupCond(cond_type, cond_unique);
    }
  }

  //-----------------------------------------------------------------------------------
  // check if all column cutter elements have a valid condition type
  // loop all column cutting elements on this processor
  for (int lid = 0; lid < nummycolele; lid++)
  {
    if (cutterele_conds_[lid].first == Inpar::XFEM::CouplingCond_NONE)
      FOUR_C_THROW("cutter element with local id %i has no valid coupling-condition", lid);
  }
}

void XFEM::CouplingBase::get_condition_by_coupling_id(
    const std::vector<Core::Conditions::Condition*>& mycond, const int coupling_id,
    std::vector<Core::Conditions::Condition*>& mynewcond)
{
  mynewcond.clear();

  // select the conditions with specified "couplingID"
  for (auto* cond : mycond)
  {
    const int id = cond->parameters().Get<int>("label");

    if (id == coupling_id) mynewcond.push_back(cond);
  }
}

void XFEM::CouplingBase::Status(const int coupling_idx, const int side_start_gid)
{
  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank_ == 0)
  {
    printf(
        "   "
        "+----------+-----------+-----------------------------+---------+--------------------------"
        "---+-----------------------------+-----------------------------+--------------------------"
        "---+\n");
    printf("   | %8i | %9i | %27s | %7i | %27s | %27s | %27s | %27s |\n", coupling_idx,
        side_start_gid, type_to_string_for_print(CondType_stringToEnum(cond_name_)).c_str(),
        coupling_id_, DisNameToString(cutter_dis_).c_str(), DisNameToString(cond_dis_).c_str(),
        DisNameToString(coupl_dis_).c_str(),
        averaging_to_string_for_print(averaging_strategy_).c_str());
  }
}



void XFEM::CouplingBase::set_averaging_strategy()
{
  const Inpar::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch (cond_type)
  {
    case Inpar::XFEM::CouplingCond_SURF_FSI_MONO:
    {
      // ask the first cutter element
      const int lid = 0;
      const int val = cutterele_conds_[lid].second->parameters().Get<int>("COUPSTRATEGY");
      averaging_strategy_ = static_cast<Inpar::XFEM::AveragingStrategy>(val);
      // check unhandled cased
      if (averaging_strategy_ == Inpar::XFEM::Mean || averaging_strategy_ == Inpar::XFEM::Harmonic)
        FOUR_C_THROW(
            "XFEM::CouplingBase::set_averaging_strategy(): Strategy Mean/Harmoninc not available "
            "for "
            "FSI monolithic, ... coming soon!");
      break;
    }
    case Inpar::XFEM::CouplingCond_SURF_FPI_MONO:
    {
      averaging_strategy_ = Inpar::XFEM::Xfluid_Sided;
      break;
    }
    case Inpar::XFEM::CouplingCond_SURF_FLUIDFLUID:
    {
      // ask the first cutter element
      const int lid = 0;
      const int val = cutterele_conds_[lid].second->parameters().Get<int>("COUPSTRATEGY");
      averaging_strategy_ = static_cast<Inpar::XFEM::AveragingStrategy>(val);
      break;
    }
    case Inpar::XFEM::CouplingCond_LEVELSET_TWOPHASE:
    case Inpar::XFEM::CouplingCond_LEVELSET_COMBUSTION:
    {
      averaging_strategy_ = Inpar::XFEM::Harmonic;
      break;
    }
    case Inpar::XFEM::CouplingCond_SURF_FSI_PART:
    case Inpar::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
    case Inpar::XFEM::CouplingCond_SURF_NEUMANN:
    case Inpar::XFEM::CouplingCond_SURF_NAVIER_SLIP:
    case Inpar::XFEM::CouplingCond_SURF_NAVIER_SLIP_TWOPHASE:
    case Inpar::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
    case Inpar::XFEM::CouplingCond_LEVELSET_NEUMANN:
    case Inpar::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
    {
      averaging_strategy_ = Inpar::XFEM::Xfluid_Sided;
      break;
    }
    default:
      FOUR_C_THROW("which is the averaging strategy for this type of coupling %i?", cond_type);
      break;
  }
}


void XFEM::CouplingBase::set_coupling_discretization()
{
  const Inpar::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch (cond_type)
  {
    case Inpar::XFEM::CouplingCond_SURF_FPI_MONO:
    {
      coupl_dis_ = cutter_dis_;
      break;
    }
    case Inpar::XFEM::CouplingCond_SURF_FSI_MONO:
    case Inpar::XFEM::CouplingCond_SURF_FLUIDFLUID:
    {
      // depending on the weighting strategy
      if (averaging_strategy_ == Inpar::XFEM::Xfluid_Sided)
      {
        coupl_dis_ = cutter_dis_;
      }
      else if (averaging_strategy_ == Inpar::XFEM::Embedded_Sided or
               averaging_strategy_ == Inpar::XFEM::Mean)
      {
        coupl_dis_ = cond_dis_;
      }
      else
        FOUR_C_THROW("Invalid coupling strategy for XFF or XFSI application");
      break;
    }
    case Inpar::XFEM::CouplingCond_LEVELSET_TWOPHASE:
    case Inpar::XFEM::CouplingCond_LEVELSET_COMBUSTION:
    {
      coupl_dis_ = bg_dis_;
      break;
    }
    case Inpar::XFEM::CouplingCond_SURF_FSI_PART:
    case Inpar::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:  // set this to Teuchos::null when the
                                                         // values are read from the function
                                                         // instead of the ivelnp vector
    case Inpar::XFEM::CouplingCond_SURF_NEUMANN:
    case Inpar::XFEM::CouplingCond_SURF_NAVIER_SLIP:
    case Inpar::XFEM::CouplingCond_SURF_NAVIER_SLIP_TWOPHASE:
    {
      coupl_dis_ = cutter_dis_;
      break;
    }
    case Inpar::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
    case Inpar::XFEM::CouplingCond_LEVELSET_NEUMANN:
    case Inpar::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
    {
      coupl_dis_ = Teuchos::null;
      break;
    }
    default:
      FOUR_C_THROW("which is the coupling discretization for this type of coupling %i?", cond_type);
      break;
  }
}

void XFEM::CouplingBase::evaluate_dirichlet_function(Core::LinAlg::Matrix<3, 1>& ivel,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond, double time)
{
  std::vector<double> final_values(3, 0.0);

  evaluate_function(final_values, x.A(), cond, time);

  ivel(0, 0) = final_values[0];
  ivel(1, 0) = final_values[1];
  ivel(2, 0) = final_values[2];
}

void XFEM::CouplingBase::evaluate_neumann_function(Core::LinAlg::Matrix<3, 1>& itraction,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond, double time)
{
  std::vector<double> final_values(3, 0.0);

  //---------------------------------------
  const auto condtype = cond->parameters().Get<std::string>("type");

  // get usual body force
  if (!(condtype == "neum_dead" or condtype == "neum_live"))
    FOUR_C_THROW("Unknown Neumann condition");
  //---------------------------------------

  evaluate_function(final_values, x.A(), cond, time);

  itraction(0, 0) = final_values[0];
  itraction(1, 0) = final_values[1];
  itraction(2, 0) = final_values[2];
}

void XFEM::CouplingBase::evaluate_neumann_function(Core::LinAlg::Matrix<6, 1>& itraction,
    const Core::LinAlg::Matrix<3, 1>& x, const Core::Conditions::Condition* cond, double time)
{
  std::vector<double> final_values(6, 0.0);

  //---------------------------------------
  const auto condtype = cond->parameters().Get<std::string>("type");

  // get usual body force
  if (!(condtype == "neum_dead" or condtype == "neum_live"))
    FOUR_C_THROW("Unknown Neumann condition");
  //---------------------------------------

  evaluate_function(final_values, x.A(), cond, time);

  for (unsigned i = 0; i < 6; ++i) itraction(i, 0) = final_values[i];
}

void XFEM::CouplingBase::evaluate_function(std::vector<double>& final_values, const double* x,
    const Core::Conditions::Condition* cond, const double time)
{
  if (cond == nullptr) FOUR_C_THROW("invalid condition");

  const int numdof = cond->parameters().Get<int>("numdof");

  if (numdof != (int)final_values.size())
    FOUR_C_THROW("you specified NUMDOF %i in the input file, however, only %i dofs allowed!",
        numdof, (int)final_values.size());

  //---------------------------------------
  // get values and switches from the condition
  const auto* onoff = &cond->parameters().Get<std::vector<int>>("onoff");
  const auto* val = &cond->parameters().Get<std::vector<double>>("val");
  const auto* functions = cond->parameters().GetIf<std::vector<int>>("funct");

  // uniformly distributed random noise

  auto& secondary = const_cast<Core::Conditions::Condition&>(*cond);
  const auto* percentage = secondary.parameters().GetIf<double>("randnoise");

  if (time < -1e-14) FOUR_C_THROW("Negative time in curve/function evaluation: time = %f", time);

  //---------------------------------------
  // set this condition
  //---------------------------------------
  for (int dof = 0; dof < numdof; ++dof)
  {
    // get factor given by spatial function
    int functnum = -1;
    if (functions) functnum = (*functions)[dof];

    // initialization of time-curve factor and function factor
    double functionfac = 1.0;

    double num = (*onoff)[dof] * (*val)[dof];

    if (functnum > 0)
    {
      functionfac = Global::Problem::Instance()
                        ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                        .evaluate(x, time, dof % numdof);
    }

    // uniformly distributed noise
    double noise = 0.0;
    if (percentage != nullptr)
    {
      const double perc = *percentage;

      if (fabs(perc) > 1e-14)
      {
        const double randomnumber = Global::Problem::Instance()
                                        ->Random()
                                        ->Uni();  // uniformly distributed between -1.0, 1.0
        noise = perc * randomnumber;
      }
    }

    final_values[dof] = num * (functionfac + noise);
  }  // loop dofs
}

void XFEM::CouplingBase::evaluate_scalar_function(double& final_values, const double* x,
    const double& val, const Core::Conditions::Condition* cond, const double time)
{
  if (cond == nullptr) FOUR_C_THROW("invalid condition");

  const int numdof = 1;

  //---------------------------------------
  // get values and switches from the condition
  const auto* function = cond->parameters().GetIf<int>("funct");

  // uniformly distributed random noise
  auto& secondary = const_cast<Core::Conditions::Condition&>(*cond);
  const auto* percentage = secondary.parameters().GetIf<double>("randnoise");

  if (time < -1e-14) FOUR_C_THROW("Negative time in curve/function evaluation: time = %f", time);

  //---------------------------------------
  // set this condition
  //---------------------------------------
  for (int dof = 0; dof < numdof; ++dof)
  {
    // get factor given by spatial function
    int functnum = -1;
    if (function) functnum = *function;

    // initialization of time-curve factor and function factor
    double functionfac = 1.0;

    double num = val;

    if (functnum > 0)
    {
      functionfac = Global::Problem::Instance()
                        ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                        .evaluate(x, time, dof % numdof);
    }

    // uniformly distributed noise
    double noise = 0.0;
    if (percentage != nullptr)
    {
      const double perc = *percentage;

      if (fabs(perc) > 1e-14)
      {
        const double randomnumber = Global::Problem::Instance()
                                        ->Random()
                                        ->Uni();  // uniformly distributed between -1.0, 1.0
        noise = perc * randomnumber;
      }
    }

    final_values = num * (functionfac + noise);
  }  // loop dofs
}

/*--------------------------------------------------------------------------*
 * get viscosity of the master fluid
 *--------------------------------------------------------------------------*/
void XFEM::CouplingBase::GetViscosityMaster(Core::Elements::Element* xfele,  ///< xfluid ele
    double& visc_m)  ///< viscosity mastersided
{
  // Get Materials of master
  Teuchos::RCP<Core::Mat::Material> mat_m;

  // Todo: As soon as the master side may not be position = outside anymore we need to take that
  // into account
  // by an additional input parameter here (e.g. XFSI with TwoPhase)
  XFEM::UTILS::get_volume_cell_material(xfele, mat_m, Core::Geo::Cut::Point::outside);
  if (mat_m->MaterialType() == Core::Materials::m_fluid)
    visc_m = Teuchos::rcp_dynamic_cast<Mat::NewtonianFluid>(mat_m)->Viscosity();
  else
    FOUR_C_THROW("get_coupling_specific_average_weights: Master Material not a fluid material?");
  return;
}

/*--------------------------------------------------------------------------*
 * get weighting paramters
 *--------------------------------------------------------------------------*/
void XFEM::CouplingBase::GetAverageWeights(Core::Elements::Element* xfele,  ///< xfluid ele
    Core::Elements::Element* coup_ele,                                      ///< coup_ele ele
    double& kappa_m,  ///< Weight parameter (parameter +/master side)
    double& kappa_s,  ///< Weight parameter (parameter -/slave  side)
    bool& non_xfluid_coupling)
{
  non_xfluid_coupling = (get_averaging_strategy() != Inpar::XFEM::Xfluid_Sided);

  if (get_averaging_strategy() != Inpar::XFEM::Harmonic)
    XFEM::UTILS::GetStdAverageWeights(get_averaging_strategy(), kappa_m);
  else
    get_coupling_specific_average_weights(xfele, coup_ele, kappa_m);

  kappa_s = 1.0 - kappa_m;
  return;
}

/*--------------------------------------------------------------------------------
 * compute viscous part of Nitsche's penalty term scaling for Nitsche's method
 *--------------------------------------------------------------------------------*/
void XFEM::CouplingBase::get_visc_penalty_stabfac(Core::Elements::Element* xfele,  ///< xfluid ele
    Core::Elements::Element* coup_ele,                                             ///< coup_ele ele
    const double& kappa_m,  ///< Weight parameter (parameter +/master side)
    const double& kappa_s,  ///< Weight parameter (parameter -/slave  side)
    const double& inv_h_k,  ///< the inverse characteristic element length h_k
    const Discret::ELEMENTS::FluidEleParameterXFEM*
        params,                     ///< parameterlist which specifies interface configuration
    double& NIT_visc_stab_fac,      ///< viscous part of Nitsche's penalty term
    double& NIT_visc_stab_fac_tang  ///< viscous part of Nitsche's penalty term in tang direction
)
{
  get_visc_penalty_stabfac(xfele, coup_ele, kappa_m, kappa_s, inv_h_k, NIT_visc_stab_fac,
      NIT_visc_stab_fac_tang, params->NITStabScaling(), params->NITStabScalingTang(),
      params->IsPseudo2D(), params->visc_stab_trac_estimate());
}

/*--------------------------------------------------------------------------------
 * compute viscous part of Nitsche's penalty term scaling for Nitsche's method
 *--------------------------------------------------------------------------------*/
void XFEM::CouplingBase::get_visc_penalty_stabfac(Core::Elements::Element* xfele,  ///< xfluid ele
    Core::Elements::Element* coup_ele,                                             ///< coup_ele ele
    const double& kappa_m,           ///< Weight parameter (parameter +/master side)
    const double& kappa_s,           ///< Weight parameter (parameter -/slave  side)
    const double& inv_h_k,           ///< the inverse characteristic element length h_k
    double& NIT_visc_stab_fac,       ///< viscous part of Nitsche's penalty term
    double& NIT_visc_stab_fac_tang,  ///< viscous part of Nitsche's penalty term in tang direction
    const double& NITStabScaling, const double& NITStabScalingTang, const bool& IsPseudo2D,
    const Inpar::XFEM::ViscStabTraceEstimate ViscStab_TraceEstimate)
{
  double penscaling = 0.0;
  if (get_averaging_strategy() != Inpar::XFEM::Embedded_Sided)
  {
    double visc_m = 0.0;
    GetViscosityMaster(
        xfele, visc_m);  // As long as mastersided we just have a fluid, directly use this ...
    penscaling = visc_m * kappa_m * inv_h_k;
  }

  if (get_averaging_strategy() != Inpar::XFEM::Xfluid_Sided)
  {
    double penscaling_s = 0.0;
    get_penalty_scaling_slave(coup_ele, penscaling_s);
    penscaling += penscaling_s * kappa_s * inv_h_k;
  }

  XFEM::UTILS::nit_compute_visc_penalty_stabfac(xfele->Shape(), penscaling, NITStabScaling,
      IsPseudo2D, ViscStab_TraceEstimate, NIT_visc_stab_fac);

  XFEM::UTILS::nit_compute_visc_penalty_stabfac(xfele->Shape(), penscaling, NITStabScalingTang,
      IsPseudo2D, ViscStab_TraceEstimate, NIT_visc_stab_fac_tang);
  return;
}

FOUR_C_NAMESPACE_CLOSE
