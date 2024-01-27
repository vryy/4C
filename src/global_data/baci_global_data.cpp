/*----------------------------------------------------------------------*/
/*! \file

\brief global list of problems

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_global_data.H"

#include "baci_comm_utils.H"
#include "baci_contact_constitutivelaw_bundle.H"
#include "baci_global_legacy_module.H"
#include "baci_inpar_problemtype.H"
#include "baci_inpar_validparameters.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_materialdefinition.H"
#include "baci_lib_discret.H"
#include "baci_lib_discret_faces.H"
#include "baci_lib_discret_hdg.H"
#include "baci_lib_utils_createdis.H"
#include "baci_particle_engine_particlereader.H"
#include "baci_rebalance.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <chrono>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<DRT::Problem*> DRT::Problem::instances_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Problem* DRT::Problem::Instance(int num)
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
void DRT::Problem::Done()
{
  // destroy singleton objects when the problem object is still alive
  for (auto* instance : instances_)
  {
    // skip null pointers arising from non-consecutive numbering of problem instances
    if (!instance) continue;
  }

  // This is called at the very end of a baci run.
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
  IO::cout.close();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Problem::Problem()
    : probtype_(ProblemType::none),
      restartstep_(0),
      restarttime_(0.0),
      communicators_(Teuchos::null)
{
  materials_ = Teuchos::rcp(new MAT::PAR::Bundle());
  contactconstitutivelaws_ = Teuchos::rcp(new CONTACT::CONSTITUTIVELAW::Bundle());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ProblemType DRT::Problem::GetProblemType() const { return probtype_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string DRT::Problem::ProblemName() const
{
  std::map<std::string, ProblemType> map = INPAR::PROBLEMTYPE::StringToProblemTypeMap();
  std::map<std::string, ProblemType>::const_iterator i;

  for (i = map.begin(); i != map.end(); ++i)
  {
    if (i->second == probtype_) return i->first;
  }
  dserror("Could not determine valid problem name");
  return "Undefined";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::Problem::Restart() const { return restartstep_; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::Problem::RestartTime() const { return restarttime_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::Problem::NDim() const
{
  const Teuchos::ParameterList& sizeparams = ProblemSizeParams();
  return sizeparams.get<int>("DIM");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::Problem::Walltime()
{
  const std::chrono::time_point<std::chrono::high_resolution_clock> now =
      std::chrono::high_resolution_clock::now();

  return std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() *
         1.0e-3;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& DRT::Problem::SolverParams(int solverNr) const
{
  std::stringstream ss;
  ss << "SOLVER " << solverNr;
  return parameters_->sublist(ss.str());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& DRT::Problem::UMFPACKSolverParams()
{
  Teuchos::ParameterList& subParams = parameters_->sublist("UMFPACK SOLVER");
  subParams.set("SOLVER", "UMFPACK");
  subParams.set("NAME", "temporary UMFPACK solver");

  return parameters_->sublist("UMFPACK SOLVER");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetCommunicators(Teuchos::RCP<CORE::COMM::Communicators> communicators)
{
  if (communicators_ != Teuchos::null) dserror("Communicators were already set.");
  communicators_ = communicators;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CORE::COMM::Communicators> DRT::Problem::GetCommunicators() const
{
  if (communicators_ == Teuchos::null) dserror("No communicators allocated yet.");
  return communicators_;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenControlFile(const Epetra_Comm& comm, const std::string& inputfile,
    std::string prefix, const std::string& restartkenner)
{
  if (Restart())
  {
    inputcontrol_ = Teuchos::rcp(new IO::InputControl(restartkenner, comm));

    if (restartstep_ < 0)
    {
      int r = IO::GetLastPossibleRestartStep(*inputcontrol_);
      SetRestartStep(r);
    }
  }

  outputcontrol_ = Teuchos::rcp(new IO::OutputControl(comm, ProblemName(),
      SpatialApproximationType(), inputfile, restartkenner, std::move(prefix), NDim(), Restart(),
      IOParams().get<int>("FILESTEPS"), INPUT::IntegralValue<int>(IOParams(), "OUTPUT_BIN"), true));

  if (!INPUT::IntegralValue<int>(IOParams(), "OUTPUT_BIN") && comm.MyPID() == 0)
  {
    IO::cout << "==================================================\n"
             << "=== WARNING: No binary output will be written. ===\n"
             << "==================================================\n"
             << IO::endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::WriteInputParameters()
{
  std::string s = OutputControlFile()->FileName();
  s.append(".parameter");
  std::ofstream stream(s.c_str());
  INPUT::PrintDatHeader(stream, *parameters_, "", false);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& parameter_list)
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
    dserror("Input parameter validation failed. Fix your input file.");
  }

  parameters_ = parameter_list;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> DRT::Problem::getValidParameters() const
{
  // call the external method to get the valid parameters
  // this way the parameter configuration is separate from the source
  return INPUT::ValidParameters();
}


Teuchos::RCP<const Teuchos::ParameterList> DRT::Problem::getParameterList() const
{
  return parameters_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::AddDis(const std::string& name, Teuchos::RCP<Discretization> dis)
{
  // safety checks
  if (dis == Teuchos::null) dserror("Received Teuchos::null.");
  if (dis->Name().empty()) dserror("Discretization has empty name string.");

  if (!discretizationmap_.insert(std::make_pair(name, dis)).second)
  {
    // if the same key already exists we have to inform the user since
    // the insert statement did not work in this case
    dserror("Could not insert discretization '%s' under (duplicate) key '%s'.", dis->Name().c_str(),
        name.c_str());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> DRT::Problem::GetDis(const std::string& name) const
{
  auto iter = discretizationmap_.find(name);

  if (iter != discretizationmap_.end())
  {
    return iter->second;
  }
  else
  {
    dserror("Could not find discretization '%s'.", name.c_str());
    return Teuchos::null;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<std::string> DRT::Problem::GetDisNames() const
{
  unsigned mysize = NumFields();
  std::vector<std::string> vec;
  vec.reserve(mysize);

  std::map<std::string, Teuchos::RCP<Discretization>>::const_iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
    vec.push_back(iter->first);
  }

  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::Problem::DoesExistDis(const std::string& name) const
{
  auto iter = discretizationmap_.find(name);
  return iter != discretizationmap_.end();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetRestartStep(int r) { restartstep_ = r; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetProblemType(ProblemType targettype) { probtype_ = targettype; }


void DRT::Problem::SetFunctionManager(CORE::UTILS::FunctionManager&& function_manager_in)
{
  functionmanager_ = std::move(function_manager_in);
}

void DRT::Problem::SetSpatialApproximationType(CORE::FE::ShapeFunctionType shape_function_type)
{
  shapefuntype_ = shape_function_type;
}

BACI_NAMESPACE_CLOSE
