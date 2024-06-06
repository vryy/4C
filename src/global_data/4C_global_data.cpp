/*----------------------------------------------------------------------*/
/*! \file

\brief global list of problems

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_global_data.hpp"

#include "4C_comm_utils.hpp"
#include "4C_contact_constitutivelaw_bundle.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_inpar_problemtype.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_discret_faces.hpp"
#include "4C_lib_discret_hdg.hpp"
#include "4C_mat_materialdefinition.hpp"
#include "4C_particle_engine_particlereader.hpp"
#include "4C_rebalance_graph_based.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <chrono>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<Global::Problem*> Global::Problem::instances_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Global::Problem* Global::Problem::Instance(int num)
{
  if (num > static_cast<int>(instances_.size()) - 1)
  {
    instances_.resize(num + 1);
    instances_[num] = new Problem();
  }
  return instances_[num];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::Done()
{
  // destroy singleton objects when the problem object is still alive
  for (auto* instance : instances_)
  {
    // skip null pointers arising from non-consecutive numbering of problem instances
    if (!instance) continue;
  }

  // This is called at the very end of a 4C run.
  //
  // It removes all global problem objects. Therefore all
  // discretizations as well and everything inside those.
  //
  // There is a whole lot going on here...
  for (auto& instance : instances_)
  {
    delete instance;
    instance = nullptr;
  }
  instances_.clear();

  // close the parallel output environment to make sure all files are properly closed
  Core::IO::cout.close();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Global::Problem::Problem()
    : probtype_(Core::ProblemType::none), restartstep_(0), communicators_(Teuchos::null)
{
  materials_ = Teuchos::rcp(new Mat::PAR::Bundle());
  contactconstitutivelaws_ = Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Bundle());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::ProblemType Global::Problem::GetProblemType() const { return probtype_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string Global::Problem::ProblemName() const
{
  std::map<std::string, Core::ProblemType> map = Inpar::PROBLEMTYPE::StringToProblemTypeMap();
  std::map<std::string, Core::ProblemType>::const_iterator i;

  for (i = map.begin(); i != map.end(); ++i)
  {
    if (i->second == probtype_) return i->first;
  }
  FOUR_C_THROW("Could not determine valid problem name");
  return "Undefined";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Global::Problem::Restart() const { return restartstep_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Global::Problem::NDim() const
{
  const Teuchos::ParameterList& sizeparams = ProblemSizeParams();
  return sizeparams.get<int>("DIM");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& Global::Problem::SolverParams(int solverNr) const
{
  std::stringstream ss;
  ss << "SOLVER " << solverNr;
  return parameters_->sublist(ss.str());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& Global::Problem::UMFPACKSolverParams()
{
  Teuchos::ParameterList& subParams = parameters_->sublist("UMFPACK SOLVER");
  subParams.set("SOLVER", "UMFPACK");
  subParams.set("NAME", "temporary UMFPACK solver");

  return parameters_->sublist("UMFPACK SOLVER");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::SetCommunicators(
    Teuchos::RCP<Core::Communication::Communicators> communicators)
{
  if (communicators_ != Teuchos::null) FOUR_C_THROW("Communicators were already set.");
  communicators_ = communicators;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Communication::Communicators> Global::Problem::GetCommunicators() const
{
  if (communicators_ == Teuchos::null) FOUR_C_THROW("No communicators allocated yet.");
  return communicators_;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::OpenControlFile(const Epetra_Comm& comm, const std::string& inputfile,
    std::string prefix, const std::string& restartkenner)
{
  if (Restart())
  {
    inputcontrol_ = Teuchos::rcp(new Core::IO::InputControl(restartkenner, comm));

    if (restartstep_ < 0)
    {
      int r = Core::IO::GetLastPossibleRestartStep(*inputcontrol_);
      SetRestartStep(r);
    }
  }

  outputcontrol_ = Teuchos::rcp(
      new Core::IO::OutputControl(comm, ProblemName(), spatial_approximation_type(), inputfile,
          restartkenner, std::move(prefix), NDim(), Restart(), IOParams().get<int>("FILESTEPS"),
          Core::UTILS::IntegralValue<bool>(IOParams(), "OUTPUT_BIN"), true));

  if (!Core::UTILS::IntegralValue<int>(IOParams(), "OUTPUT_BIN") && comm.MyPID() == 0)
  {
    Core::IO::cout << "==================================================\n"
                   << "=== WARNING: No binary output will be written. ===\n"
                   << "==================================================\n"
                   << Core::IO::endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::write_input_parameters()
{
  std::string s = OutputControlFile()->FileName();
  s.append(".parameter");
  std::ofstream stream(s.c_str());
  Input::PrintDatHeader(stream, *parameters_, "", false);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& parameter_list)
{
  try
  {
    // Test parameter list against valid parameters, set default values
    // and set validator objects to extract numerical values for string
    // parameters.
    parameter_list->validateParametersAndSetDefaults(*getValidParameters());
  }
  catch (Teuchos::Exceptions::InvalidParameter& err)
  {
    std::cerr << "\n\n" << err.what();
    FOUR_C_THROW("Input parameter validation failed. Fix your input file.");
  }

  parameters_ = parameter_list;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> Global::Problem::getValidParameters() const
{
  // call the external method to get the valid parameters
  // this way the parameter configuration is separate from the source
  return Input::ValidParameters();
}


Teuchos::RCP<const Teuchos::ParameterList> Global::Problem::getParameterList() const
{
  return parameters_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::AddDis(const std::string& name, Teuchos::RCP<Discret::Discretization> dis)
{
  // safety checks
  if (dis == Teuchos::null) FOUR_C_THROW("Received Teuchos::null.");
  if (dis->Name().empty()) FOUR_C_THROW("discretization has empty name string.");

  if (!discretizationmap_.insert(std::make_pair(name, dis)).second)
  {
    // if the same key already exists we have to inform the user since
    // the insert statement did not work in this case
    FOUR_C_THROW("Could not insert discretization '%s' under (duplicate) key '%s'.",
        dis->Name().c_str(), name.c_str());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Discret::Discretization> Global::Problem::GetDis(const std::string& name) const
{
  auto iter = discretizationmap_.find(name);

  if (iter != discretizationmap_.end())
  {
    return iter->second;
  }
  else
  {
    FOUR_C_THROW("Could not find discretization '%s'.", name.c_str());
    return Teuchos::null;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<std::string> Global::Problem::GetDisNames() const
{
  unsigned mysize = NumFields();
  std::vector<std::string> vec;
  vec.reserve(mysize);

  std::transform(discretizationmap_.begin(), discretizationmap_.end(), std::back_inserter(vec),
      [](const auto& key_value) { return key_value.first; });

  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Global::Problem::DoesExistDis(const std::string& name) const
{
  auto iter = discretizationmap_.find(name);
  return iter != discretizationmap_.end();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::SetRestartStep(int r) { restartstep_ = r; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::Problem::SetProblemType(Core::ProblemType targettype) { probtype_ = targettype; }


void Global::Problem::SetFunctionManager(Core::UTILS::FunctionManager&& function_manager_in)
{
  functionmanager_ = std::move(function_manager_in);
}

void Global::Problem::set_spatial_approximation_type(
    Core::FE::ShapeFunctionType shape_function_type)
{
  shapefuntype_ = shape_function_type;
}

FOUR_C_NAMESPACE_CLOSE
