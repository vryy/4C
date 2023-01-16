/*----------------------------------------------------------------------*/
/*! \file

\brief global list of problems

\level 1


*/
/*----------------------------------------------------------------------*/

#include <sys/time.h>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_ParameterListExceptions.hpp>

#include <Epetra_Time.h>
#include <Epetra_Comm.h>

#include "lib_conditiondefinition.H"
#include "lib_dofset_independent.H"
#include "lib_materialdefinition.H"
#include "lib_function.H"
#include "lib_globalproblem.H"
#include "lib_inputreader.H"
#include "lib_elementreader.H"
#include "lib_meshreader.H"
#include "lib_nodereader.H"
#include "lib_particlereader.H"
#include "lib_utils_createdis.H"
#include "lib_discret.H"
#include "lib_discret_faces.H"
#include "lib_discret_xwall.H"
#include "lib_discret_hdg.H"
#include "lib_discret_xfem.H"
#include "lib_linedefinition.H"
#include "contact_constitutivelaw_constitutivelaw_definition.H"
#include "contact_constitutivelaw_bundle.H"
#include "inpar_validconditions.H"
#include "inpar_validparameters.H"
#include "inpar_validmaterials.H"
#include "mat_elchmat.H"
#include "mat_elchphase.H"
#include "mat_micromaterial.H"
#include "mat_newman_multiscale.H"
#include "mat_scatra_mat_multiscale.H"
#include "rebalance_utils.H"
#include "comm_utils.H"
#include "inpar_validcontactconstitutivelaw.H"
#include "inpar_problemtype.H"
#include "io.H"
#include "io_control.H"

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
  struct timeval tp;
  gettimeofday(&tp, nullptr);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string DRT::Problem::SpatialApproximation() const
{
  return INPAR::PROBLEMTYPE::ShapeFunctionTypeToString(shapefuntype_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadParameter(DRT::INPUT::DatFileReader& reader)
{
  Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList("DAT FILE"));

  reader.ReadGidSection("--DISCRETISATION", *list);
  reader.ReadGidSection("--PROBLEM SIZE", *list);
  reader.ReadGidSection("--PROBLEM TYP", *list);
  reader.ReadGidSection("--BINNING STRATEGY", *list);
  reader.ReadGidSection("--BOUNDINGVOLUME STRATEGY", *list);
  reader.ReadGidSection("--IO", *list);
  reader.ReadGidSection("--IO/EVERY ITERATION", *list);
  reader.ReadGidSection("--IO/MONITOR STRUCTURE DBC", *list);
  reader.ReadGidSection("--IO/RUNTIME VTK OUTPUT", *list);
  reader.ReadGidSection("--IO/RUNTIME VTK OUTPUT/FLUID", *list);
  reader.ReadGidSection("--IO/RUNTIME VTK OUTPUT/STRUCTURE", *list);
  reader.ReadGidSection("--IO/RUNTIME VTK OUTPUT/BEAMS", *list);
  reader.ReadGidSection("--IO/RUNTIME VTP OUTPUT STRUCTURE", *list);
  reader.ReadGidSection("--DESIGN DESCRIPTION", *list);
  reader.ReadGidSection("--PATIENT SPECIFIC", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GEMM", *list);
  reader.ReadGidSection("--STAT INVERSE ANALYSIS", *list);
  reader.ReadGidSection("--MORTAR COUPLING", *list);
  reader.ReadGidSection("--MORTAR COUPLING/PARALLEL REDISTRIBUTION", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED/COMBO", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED/STEEPESTASCENT", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED/LAGRANGE_MULTIPLIER_FUNCTION", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED/PLOT", *list);
  reader.ReadGidSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING", *list);
  reader.ReadGidSection(
      "--CARDIOVASCULAR 0D-STRUCTURE COUPLING/SYS-PUL CIRCULATION PARAMETERS", *list);
  reader.ReadGidSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING/RESPIRATORY PARAMETERS", *list);
  reader.ReadGidSection("--BROWNIAN DYNAMICS", *list);
  reader.ReadGidSection("--BEAM INTERACTION", *list);
  reader.ReadGidSection("--BEAM INTERACTION/SPHERE BEAM LINK", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO BEAM CONTACT", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO SPHERE CONTACT", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE CONTACT", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE/RUNTIME VTK OUTPUT", *list);
  reader.ReadGidSection("--BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING", *list);
  reader.ReadGidSection(
      "--BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING/RUNTIME VTK OUTPUT", *list);
  reader.ReadGidSection("--BEAM INTERACTION/CROSSLINKING", *list);
  reader.ReadGidSection("--THERMAL DYNAMIC", *list);
  reader.ReadGidSection("--THERMAL DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--THERMAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--TSI DYNAMIC", *list);
  reader.ReadGidSection("--TSI DYNAMIC/MONOLITHIC", *list);
  reader.ReadGidSection("--TSI DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--TSI CONTACT", *list);
  reader.ReadGidSection("--POROELASTICITY DYNAMIC", *list);
  reader.ReadGidSection("--POROSCATRA CONTROL", *list);
  reader.ReadGidSection("--POROFLUIDMULTIPHASE DYNAMIC", *list);
  reader.ReadGidSection("--POROFLUIDMULTIPHASE DYNAMIC/ARTERY COUPLING", *list);
  reader.ReadGidSection("--POROMULTIPHASE DYNAMIC", *list);
  reader.ReadGidSection("--POROMULTIPHASE DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--POROMULTIPHASE DYNAMIC/MONOLITHIC", *list);
  reader.ReadGidSection("--POROMULTIPHASESCATRA DYNAMIC", *list);
  reader.ReadGidSection("--POROMULTIPHASESCATRA DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--POROMULTIPHASESCATRA DYNAMIC/MONOLITHIC", *list);
  reader.ReadGidSection("--ELASTO HYDRO DYNAMIC", *list);
  reader.ReadGidSection("--ELASTO HYDRO DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--ELASTO HYDRO DYNAMIC/MONOLITHIC", *list);
  reader.ReadGidSection("--SSI CONTROL", *list);
  reader.ReadGidSection("--SSI CONTROL/ELCH", *list);
  reader.ReadGidSection("--SSI CONTROL/MANIFOLD", *list);
  reader.ReadGidSection("--SSI CONTROL/MONOLITHIC", *list);
  reader.ReadGidSection("--SSI CONTROL/PARTITIONED", *list);
  reader.ReadGidSection("--SSTI CONTROL", *list);
  reader.ReadGidSection("--SSTI CONTROL/MONOLITHIC", *list);
  reader.ReadGidSection("--SSTI CONTROL/THERMO", *list);
  reader.ReadGidSection("--FLUID DYNAMIC", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/RESIDUAL-BASED STABILIZATION", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/EDGE-BASED STABILIZATION", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/POROUS-FLOW STABILIZATION", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TURBULENCE MODEL", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/SUBGRID VISCOSITY", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/WALL MODEL", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/MULTIFRACTAL SUBGRID SCALES", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/TURBULENT INFLOW", *list);
  reader.ReadGidSection("--FLUID DYNAMIC/NONLINEAR SOLVER TOLERANCES", *list);
  reader.ReadGidSection("--TWO PHASE FLOW", *list);
  reader.ReadGidSection("--TWO PHASE FLOW/SURFACE TENSION", *list);
  reader.ReadGidSection("--TWO PHASE FLOW/SMEARED", *list);
  reader.ReadGidSection("--LUBRICATION DYNAMIC", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/NONLINEAR", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/S2I COUPLING", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/ARTERY COUPLING", *list);
  reader.ReadGidSection("--STI DYNAMIC", *list);
  reader.ReadGidSection("--STI DYNAMIC/MONOLITHIC", *list);
  reader.ReadGidSection("--STI DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--FS3I DYNAMIC", *list);
  reader.ReadGidSection("--FS3I DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION", *list);
  reader.ReadGidSection("--FS3I DYNAMIC/AC", *list);
  reader.ReadGidSection("--ALE DYNAMIC", *list);
  reader.ReadGidSection("--FSI DYNAMIC", *list);
  reader.ReadGidSection("--FSI DYNAMIC/CONSTRAINT", *list);
  reader.ReadGidSection("--FSI DYNAMIC/MONOLITHIC SOLVER", *list);
  reader.ReadGidSection("--FSI DYNAMIC/PARTITIONED SOLVER", *list);
  reader.ReadGidSection("--FSI DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--FLUID BEAM INTERACTION", *list);
  reader.ReadGidSection("--FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING", *list);
  reader.ReadGidSection(
      "--FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING/RUNTIME VTK OUTPUT", *list);
  reader.ReadGidSection("--IMMERSED METHOD", *list);
  reader.ReadGidSection("--IMMERSED METHOD/PARTITIONED SOLVER", *list);
  reader.ReadGidSection("--FPSI DYNAMIC", *list);
  reader.ReadGidSection("--ARTERIAL DYNAMIC", *list);
  reader.ReadGidSection("--REDUCED DIMENSIONAL AIRWAYS DYNAMIC", *list);
  reader.ReadGidSection("--COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC", *list);
  reader.ReadGidSection("--SEARCH TREE", *list);
  reader.ReadGidSection("--XFEM GENERAL", *list);
  reader.ReadGidSection("--CUT GENERAL", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC/GENERAL", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC/XFPSI MONOLITHIC", *list);
  reader.ReadGidSection("--LOMA CONTROL", *list);
  reader.ReadGidSection("--ELCH CONTROL", *list);
  reader.ReadGidSection("--ELCH CONTROL/DIFFCOND", *list);
  reader.ReadGidSection("--BIOFILM CONTROL", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL/TOPOLOGY OPTIMIZER", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL/TOPOLOGY ADJOINT FLUID", *list);
  reader.ReadGidSection("--PARTICLE DYNAMIC", *list);
  reader.ReadGidSection("--PARTICLE DYNAMIC/INITIAL AND BOUNDARY CONDITIONS", *list);
  reader.ReadGidSection("--PARTICLE DYNAMIC/SPH", *list);
  reader.ReadGidSection("--PARTICLE DYNAMIC/DEM", *list);
  reader.ReadGidSection("--PASI DYNAMIC", *list);
  reader.ReadGidSection("--LEVEL-SET CONTROL", *list);
  reader.ReadGidSection("--LEVEL-SET CONTROL/REINITIALIZATION", *list);
  reader.ReadGidSection("--WEAR", *list);
  reader.ReadGidSection("--BEAM CONTACT", *list);
  reader.ReadGidSection("--BEAM CONTACT/RUNTIME VTK OUTPUT", *list);
  reader.ReadGidSection("--BEAM POTENTIAL", *list);
  reader.ReadGidSection("--BEAM POTENTIAL/RUNTIME VTK OUTPUT", *list);
  reader.ReadGidSection("--SEMI-SMOOTH PLASTICITY", *list);
  reader.ReadGidSection("--ELECTROMAGNETIC DYNAMIC", *list);
  reader.ReadGidSection("--VOLMORTAR COUPLING", *list);
  reader.ReadGidSection("--TUTORIAL DYNAMIC", *list);
  reader.ReadGidSection("--TUTORIAL DYNAMIC/NONLINEAR TRUSS", *list);
  reader.ReadGidSection("--TUTORIAL DYNAMIC/FIXED POINT SCHEME", *list);
  reader.ReadGidSection("--CARDIAC MONODOMAIN CONTROL", *list);
  reader.ReadGidSection("--MOR", *list);
  reader.ReadGidSection("--MESH PARTITIONING", *list);

  reader.ReadSection("--STRUCT NOX", *list);
  reader.ReadSection("--STRUCT NOX/Direction", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton/Modified", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton/Linear Solver", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Steepest Descent", *list);
  reader.ReadSection("--STRUCT NOX/Line Search", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Full Step", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Backtrack", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Polynomial", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/More'-Thuente", *list);
  reader.ReadSection("--STRUCT NOX/Pseudo Transient", *list);
  reader.ReadSection("--STRUCT NOX/Trust Region", *list);
  reader.ReadSection("--STRUCT NOX/Printing", *list);
  reader.ReadSection("--STRUCT NOX/Status Test", *list);
  reader.ReadSection("--STRUCT NOX/Solver Options", *list);

  reader.ReadSection("--LOCA", *list);
  reader.ReadSection("--LOCA/Stepper", *list);
  reader.ReadSection("--LOCA/Stepper/Eigensolver", *list);
  reader.ReadSection("--LOCA/Bifurcation", *list);
  reader.ReadSection("--LOCA/Predictor", *list);
  reader.ReadSection("--LOCA/Predictor/First Step Predictor", *list);
  reader.ReadSection("--LOCA/Predictor/Last Step Predictor", *list);
  reader.ReadSection("--LOCA/Step Size", *list);
  reader.ReadSection("--LOCA/Constraints", *list);

  // read in solver sections
  // Note: the maximum number of solver blocks in dat files is hardwired here.
  // If you change this do not forget to edit the corresponding parts in
  // validparameters.cpp, too!
  for (int i = 1; i < 10; i++)
  {
    std::stringstream ss;
    ss << "--SOLVER " << i;
    reader.ReadGidSection(ss.str(), *list);

    // adapt path of XML file if necessary
    Teuchos::ParameterList& sublist = list->sublist(ss.str().substr(2));
    std::vector<std::string> listOfFileNameParameters = {"AMGNXN_XML_FILE", "MUELU_XML_FILE"};

    for (auto& filenameParameter : listOfFileNameParameters)
    {
      auto* xml_filename = sublist.getPtr<std::string>(filenameParameter);
      if (xml_filename != nullptr and *xml_filename != "none")
      {
        // make path relative to input file path if it is not an absolute path
        if ((*xml_filename)[0] != '/')
        {
          std::string filename = reader.MyInputfileName();
          std::string::size_type pos = filename.rfind('/');
          if (pos != std::string::npos)
          {
            std::string tmp = filename.substr(0, pos + 1);
            xml_filename->insert(xml_filename->begin(), tmp.begin(), tmp.end());
          }
        }
      }
    }
  }

  reader.ReadGidSection("--UMFPACK SOLVER", *list);


  // read in STRUCT NOX/Status Test and modify the xml file name, if there
  // is one.
  if (list->sublist("STRUCT NOX").sublist("Status Test").isParameter("XML File"))
  {
    // adapt path of XML file if necessary
    Teuchos::ParameterList& sublist = list->sublist("STRUCT NOX").sublist("Status Test");
    auto* statustest_xmlfile = sublist.getPtr<std::string>("XML File");
    // make path relative to input file path if it is not an absolute path
    if (((*statustest_xmlfile)[0] != '/') and ((*statustest_xmlfile) != "none"))
    {
      std::string filename = reader.MyInputfileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string tmp = filename.substr(0, pos + 1);
        statustest_xmlfile->insert(statustest_xmlfile->begin(), tmp.begin(), tmp.end());
      }
    }
  }  // STRUCT NOX/Status Test

  // check for invalid parameters
  setParameterList(list);

  //---------------------------------------------------------------------
  // Now we have successfully read the whole input file. It's time to access some data

  // 1) get the problem type
  const Teuchos::ParameterList& type = ProblemTypeParams();
  probtype_ = DRT::INPUT::IntegralValue<ProblemType>(type, "PROBLEMTYP");

  // 2) get the spatial approximation type
  shapefuntype_ = DRT::INPUT::IntegralValue<ShapeFunctionType>(type, "SHAPEFCT");

  // 3) do the restart business with the four options we support (partially)
  if (restartstep_ == 0)
  {
    // no restart flag on the command line, so check the restart flag from the input file
    restartstep_ = type.get<int>("RESTART");
  }
  else  // SetRestartStep() has been called before!
  {
    // There is a non-zero restart flag on the command line, so we ignore the input file.
    // The RESTART flag in the input file should be zero or have the same value!
    const int restartflaginfile = type.get<int>("RESTART");
    if ((restartflaginfile > 0) and (restartflaginfile != restartstep_))
      dserror("Restart flags in input file and command line are non-zero and different!");
  }

  // Set restart time
  restarttime_ = type.get<double>("RESTARTTIME");
  if (restarttime_ > 0.0)
  {
    // Currently this option is implemented only for scalar structure interaction problems
    if (GetProblemType() != ProblemType::ssi)
      dserror("Restart with time option currently only implemented for SSI problems");
    // The value restartstep_ is used very deep down in Baci. Therefore we demand the user
    // to set this value, if one wants to use restart time option.
    // If this feature should be expanded for all problemtyps another handling of the
    // restartstep_ has to be considered
    if (restartstep_ == 0)
      dserror("Restart with time option needs a RESTART flag different from 0");
  }

  // Set restart time based on walltime
  const double restartinterval = IOParams().get<double>("RESTARTWALLTIMEINTERVAL");
  const int restartevry = IOParams().get<int>("RESTARTEVRY");
  RestartManager()->SetupRestartManager(restartinterval, restartevry);

  // 4) set random seed
  // time is in seconds, therefore we add the global processor id to obtain a unique seed on each
  // proc
  {
    int rs = type.get<int>("RANDSEED");
    if (rs < 0)
      rs = static_cast<int>(time(nullptr)) +
           42 * DRT::Problem::Instance(0)->GetCommunicators()->GlobalComm()->MyPID();

    srand((unsigned int)rs);  // Set random seed for stdlibrary. This is deprecated, as it does not
                              // produce random numbers on some platforms!
    Random()->SetRandSeed((unsigned int)rs);  // Use this instead.
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& DRT::Problem::SolverParams(int solverNr) const
{
  std::stringstream ss;
  ss << "SOLVER " << solverNr;
  return getParameterList()->sublist(ss.str());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& DRT::Problem::UMFPACKSolverParams()
{
  Teuchos::RCP<Teuchos::ParameterList> params = getNonconstParameterList();

  Teuchos::ParameterList& subParams = params->sublist("UMFPACK SOLVER");
  subParams.set("SOLVER", "UMFPACK");
  subParams.set("NAME", "temporary UMFPACK solver");

  return getParameterList()->sublist("UMFPACK SOLVER");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetCommunicators(Teuchos::RCP<COMM_UTILS::Communicators> communicators)
{
  if (communicators_ != Teuchos::null) dserror("Communicators were already set.");
  communicators_ = communicators;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<COMM_UTILS::Communicators> DRT::Problem::GetCommunicators() const
{
  if (communicators_ == Teuchos::null) dserror("No communicators allocated yet.");
  return communicators_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMaterials(DRT::INPUT::DatFileReader& reader)
{
  // create list of known materials
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>> vm =
      DRT::INPUT::ValidMaterials();
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& matlist = *vm;

  // test for each material definition (input file --MATERIALS section)
  // and store in #matmap_
  for (auto& mat : matlist)
  {
    // read material from DAT file of type #matlist[m]
    mat->Read(*this, reader, materials_);
  }

  // check if every material was identified
  const std::string material_section = "--MATERIALS";
  std::vector<const char*> section = reader.Section(material_section);
  int nummat = 0;
  if (section.size() > 0)
  {
    for (auto& section_i : section)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(section_i));

      std::string mat;
      std::string number;
      std::string name;
      (*condline) >> mat >> number >> name;
      if ((not(*condline)) or (mat != "MAT"))
        dserror("invalid material line in '%s'", name.c_str());

      // extract material ID
      int matid = -1;
      {
        char* ptr;
        matid = static_cast<int>(strtol(number.c_str(), &ptr, 10));
        if (ptr == number.c_str())
          dserror("failed to read material object number '%s'", number.c_str());
      }

      // processed?
      if (materials_->Find(matid) == -1)
        dserror("Material 'MAT %d' with name '%s' could not be identified", matid, name.c_str());

      // count number of materials provided in file
      nummat += 1;
    }
  }

  // make fast access parameters
  materials_->MakeParameters();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadContactConstitutiveLaws(DRT::INPUT::DatFileReader& reader)
{
  // create list of known contact constitutive laws
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>> vm =
      DRT::INPUT::ValidContactConstitutiveLaws();
  std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist = *vm;

  // test for each contact constitutive law definition (input file --CONTACT CONSTITUTIVE LAWS
  // section) and store it
  for (auto& m : coconstlawlist)
  {
    // read contact constitutive law from DAT file of type
    m->Read(*this, reader, contactconstitutivelaws_);
  }

  // check if every contact constitutive law was identified
  const std::string contact_const_laws = "--CONTACT CONSTITUTIVE LAWS";
  std::vector<const char*> section = reader.Section(contact_const_laws);
  int numlaws = 0;
  if (section.size() > 0)
  {
    for (auto& section_i : section)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(section_i));

      std::string coconstlaw;
      std::string number;
      std::string name;
      (*condline) >> coconstlaw >> number >> name;
      if ((not(*condline)) or (coconstlaw != "LAW"))
        dserror("invalid contact constitutive law line in '%s'", name.c_str());

      // extract contact constitutive law ID
      int coconstlawid = -1;
      {
        char* ptr;
        coconstlawid = static_cast<int>(strtol(number.c_str(), &ptr, 10));
        if (ptr == number.c_str())
          dserror("failed to read contact constitutive law object number '%s'", number.c_str());
      }

      // processed?
      if (contactconstitutivelaws_->Find(coconstlawid) == -1)
        dserror("Contact constitutive law 'LAW %d' with name '%s' could not be identified",
            coconstlawid, name.c_str());

      // count number of contact constitutive laws provided in file
      numlaws += 1;
    }
  }

  // make fast access parameters
  contactconstitutivelaws_->MakeParameters();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadCloningMaterialMap(DRT::INPUT::DatFileReader& reader)
{
  Teuchos::RCP<DRT::INPUT::Lines> lines = DRT::UTILS::ValidCloningMaterialMapLines();

  // perform the actual reading and extract the input parameters
  std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> input_line_vec = lines->Read(reader);
  for (auto& input_line : input_line_vec)
  {
    // extract what was read from the input file
    std::string src_field;
    input_line->ExtractString("SRC_FIELD", src_field);
    int src_matid(-1);
    input_line->ExtractInt("SRC_MAT", src_matid);
    std::string tar_field;
    input_line->ExtractString("TAR_FIELD", tar_field);
    int tar_matid(-1);
    input_line->ExtractInt("TAR_MAT", tar_matid);

    // create the key pair
    std::pair<std::string, std::string> fields(src_field, tar_field);

    // enter the material pairing into the map
    std::pair<int, int> matmap(src_matid, tar_matid);
    clonefieldmatmap_[fields].insert(matmap);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadTimeFunctionResult(DRT::INPUT::DatFileReader& reader)
{
  //---------------------------------------- input of spatial functions
  functionmanager_.ReadInput(reader);
  //-------------------------------------- input of result descriptions
  resulttest_.ReadInput(reader);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadConditions(DRT::INPUT::DatFileReader& reader)
{
  Epetra_Time time(*reader.Comm());
  if (reader.Comm()->MyPID() == 0)
  {
    IO::cout << "Read/generate conditions                          in....";
    IO::cout.flush();
  }

  //------------------------------- read number of design objects we have
  // this currently serves to determine how many node sets we might have
  const Teuchos::ParameterList& design = DesignDescriptionParams();
  int ndnode = design.get<int>("NDPOINT");
  int ndline = design.get<int>("NDLINE");
  int ndsurf = design.get<int>("NDSURF");
  int ndvol = design.get<int>("NDVOL");

  //--------------------------------------------- read generic node sets
  // read design nodes <-> nodes
  std::vector<std::vector<int>> dnode_fenode(ndnode);
  reader.ReadDesign("DNODE", dnode_fenode);

  // read design lines <-> nodes
  std::vector<std::vector<int>> dline_fenode(ndline);
  reader.ReadDesign("DLINE", dline_fenode);

  // read design surfaces <-> nodes
  std::vector<std::vector<int>> dsurf_fenode(ndsurf);
  reader.ReadDesign("DSURF", dsurf_fenode);

  // read design volumes <-> nodes
  std::vector<std::vector<int>> dvol_fenode(ndvol);
  reader.ReadDesign("DVOL", dvol_fenode);

  // check for meshfree discretisation to add node set topologies
  std::vector<std::vector<std::vector<int>>*> nodeset(4);
  nodeset[0] = &dnode_fenode;
  nodeset[1] = &dline_fenode;
  nodeset[2] = &dsurf_fenode;
  nodeset[3] = &dvol_fenode;

  // create list of known conditions
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>> vc =
      DRT::INPUT::ValidConditions();
  std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist = *vc;

  // test for each condition definition (input file condition section)
  // - read all conditions that match the definition
  // - add the nodal clouds to the conditions
  // - add the conditions to the appropriate discretizations
  //
  // Note that this will reset (un-FillComplete) the discretizations.
  for (auto& condition : condlist)
  {
    std::multimap<int, Teuchos::RCP<DRT::Condition>> cond;

    // read conditions from dat file
    condition->Read(*this, reader, cond);

    // add nodes to conditions
    std::multimap<int, Teuchos::RCP<DRT::Condition>>::const_iterator curr;
    for (curr = cond.begin(); curr != cond.end(); ++curr)
    {
      switch (curr->second->GType())
      {
        case Condition::Point:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dnode_fenode.size())
          {
            dserror(
                "DPoint %d not in range [0:%d[\n"
                "DPoint condition on non existent DPoint?",
                curr->first, dnode_fenode.size());
          }
          curr->second->Add("Node Ids", dnode_fenode[curr->first]);
          break;
        case Condition::Line:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dline_fenode.size())
          {
            dserror(
                "DLine %d not in range [0:%d[\n"
                "DLine condition on non existent DLine?",
                curr->first, dline_fenode.size());
          }
          curr->second->Add("Node Ids", dline_fenode[curr->first]);
          break;
        case Condition::Surface:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dsurf_fenode.size())
          {
            dserror(
                "DSurface %d not in range [0:%d[\n"
                "DSurface condition on non existent DSurface?",
                curr->first, dsurf_fenode.size());
          }
          curr->second->Add("Node Ids", dsurf_fenode[curr->first]);
          break;
        case Condition::Volume:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dvol_fenode.size())
          {
            dserror(
                "DVolume %d not in range [0:%d[\n"
                "DVolume condition on non existent DVolume?",
                curr->first, dvol_fenode.size());
          }
          curr->second->Add("Node Ids", dvol_fenode[curr->first]);
          break;
        default:
          dserror("geometry type unspecified");
          break;
      }

      // Iterate through all discretizations and sort the appropriate condition
      // into the correct discretization it applies to

      std::map<std::string, Teuchos::RCP<Discretization>>::iterator iter;
      for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
      {
        Teuchos::RCP<DRT::Discretization> actdis = iter->second;

        const std::vector<int>* nodes = curr->second->Nodes();
        if (nodes->size() == 0)
          dserror("%s condition %d has no nodal cloud", condition->Description().c_str(),
              curr->second->Id());

        int foundit = 0;
        for (int node : *nodes)
        {
          foundit = actdis->HaveGlobalNode(node);
          if (foundit) break;
        }
        int found = 0;
        actdis->Comm().SumAll(&foundit, &found, 1);
        if (found)
        {
          // Insert a copy since we might insert the same condition in many discretizations.
          actdis->SetCondition(condition->Name(), Teuchos::rcp(new Condition(*curr->second)));
        }
      }
    }
  }

  if (reader.Comm()->MyPID() == 0)
  {
    std::cout << time.ElapsedTime() << " secs\n";
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadKnots(DRT::INPUT::DatFileReader& reader)
{
  // get information on the spatial approximation --- we only read knots
  // in the nurbs case
  ShapeFunctionType distype = SpatialApproximationType();

  // get problem dimension
  int dim = NDim();

  // Iterate through all discretizations and sort the appropriate condition
  // into the correct discretization it applies to

  std::map<std::string, Teuchos::RCP<Discretization>>::iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
    Teuchos::RCP<DRT::Discretization> actdis = iter->second;

    if (distype == ShapeFunctionType::shapefunction_nurbs)
    {
      // cast discretisation to nurbs variant to be able
      // to add the knotvector
      auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*actdis));

      if (nurbsdis == nullptr)
        dserror("Discretization %s is not a NurbsDiscretization! Panic.", actdis->Name().c_str());

      // define an empty knot vector object
      Teuchos::RCP<DRT::NURBS::Knotvector> disknots = Teuchos::null;

      // read the knotvector data from the input
      reader.ReadKnots(dim, actdis->Name(), disknots);

      if (disknots == Teuchos::null)
      {
        dserror("Knotvector read failed in Nurbs discretisation\n");
      }

      // make sure atdis is fillcompleted, to be able to call
      // ElementRowMap() on it
      // do not initialize elements, since this would require knot
      // vector values
      if (!actdis->Filled())
      {
        actdis->FillComplete(false, false, false);
      }

      // the smallest gid in the discretisation determines the access
      // pattern via the element offset
      int smallest_gid_in_dis = actdis->ElementRowMap()->MinAllGID();

      // consistency checks
      disknots->FinishKnots(smallest_gid_in_dis);

      // add knots to discretisation
      nurbsdis->SetKnotVector(disknots);
    }
  }  // loop fields
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadParticles(DRT::INPUT::DatFileReader& reader)
{
  // no need to read in particles in case of restart
  if (Restart()) return;

  // the basic particle reader
  DRT::INPUT::ParticleReader particlereader(reader, "--PARTICLES");

  // do the actual reading of particles
  particlereader.Read(particles_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenControlFile(const Epetra_Comm& comm, const std::string& inputfile,
    std::string prefix, const std::string& restartkenner)
{
  if (Restart()) inputcontrol_ = Teuchos::rcp(new IO::InputControl(restartkenner, comm));

  outputcontrol_ =
      Teuchos::rcp(new IO::OutputControl(comm, ProblemName(), SpatialApproximationType(), inputfile,
          restartkenner, std::move(prefix), NDim(), Restart(), IOParams().get<int>("FILESTEPS"),
          DRT::INPUT::IntegralValue<int>(IOParams(), "OUTPUT_BIN"), true));

  if (!DRT::INPUT::IntegralValue<int>(IOParams(), "OUTPUT_BIN") && comm.MyPID() == 0)
  {
    IO::cout << "==================================================\n"
             << "=== WARNING: No binary output will be written. ===\n"
             << "==================================================\n"
             << IO::endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenErrorFile(
    const Epetra_Comm& comm, std::string prefix, const bool enforceopening)
{
  bool openfile = enforceopening;
  if (!enforceopening)
  {
    // what's given in the input file?
    openfile = DRT::INPUT::IntegralValue<int>(IOParams(), "OUTPUT_BIN");
  }
  errorfilecontrol_ =
      Teuchos::rcp(new IO::ErrorFileControl(comm, std::move(prefix), Restart(), openfile));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::WriteInputParameters()
{
  std::string s = OutputControlFile()->FileName();
  s.append(".parameter");
  std::ofstream stream(s.c_str());
  DRT::INPUT::PrintDatHeader(stream, *getParameterList(), "", false, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadFields(DRT::INPUT::DatFileReader& reader, const bool readmesh)
{
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> xfluiddis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> structaledis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> thermdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> lubricationdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> cellscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> fluidscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> structscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> artscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> arterydis = Teuchos::null;  //_1D_ARTERY_
  Teuchos::RCP<DRT::Discretization> airwaydis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> optidis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> porofluiddis = Teuchos::null;  // fpsi, poroelast
  Teuchos::RCP<DRT::Discretization> elemagdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> celldis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> pboxdis = Teuchos::null;

  // decide which kind of spatial representation is required
  const ShapeFunctionType distype = SpatialApproximationType();

  // the basic mesh reader. now add desired node and element readers to it!
  DRT::INPUT::MeshReader meshreader(reader, "--NODE COORDS");
  meshreader.AddNodeReader(Teuchos::rcp(new DRT::INPUT::NodeReader(reader, "--NODE COORDS")));

  switch (GetProblemType())
  {
    case ProblemType::fsi:
    case ProblemType::fsi_redmodels:
    case ProblemType::fsi_lung:
    {
      if (distype == ShapeFunctionType::shapefunction_nurbs)
      {
        structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
      }
      else if (DRT::INPUT::IntegralValue<int>(FluidDynamicParams().sublist("WALL MODEL"), "X_WALL"))
      {
        structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXWall("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }
      else
      {
        structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
        if (DRT::INPUT::IntegralValue<bool>(
                XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
          xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      if (xfluiddis != Teuchos::null)
        xfluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);
      if (xfluiddis != Teuchos::null) AddDis("xfluid", xfluiddis);
      AddDis("ale", aledis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

      if (xfluiddis != Teuchos::null)
      {
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS")));
      }
      else
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

      break;
    }
    case ProblemType::gas_fsi:
    case ProblemType::ac_fsi:
    case ProblemType::thermo_fsi:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          dserror("Nurbs discretization not possible for fs3i!");
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1", reader.Comm()));
          structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));
      fluidscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis)));
      structscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structscatradis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);
      AddDis("ale", aledis);
      AddDis("scatra1", fluidscatradis);
      AddDis("scatra2", structscatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(fluidscatradis, reader, "--TRANSPORT ELEMENTS")));
      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(structscatradis, reader, "--TRANSPORT2 ELEMENTS")));

#ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false, false, false);
#endif

      break;
    }
    case ProblemType::biofilm_fsi:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          dserror("Nurbs discretization not possible for biofilm problems!");
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          structaledis = Teuchos::rcp(new DRT::Discretization("structale", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));
      structaledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structaledis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);
      AddDis("ale", aledis);
      AddDis("structale", structaledis);


      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

#ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false, false, false);
#endif

      // fluid scatra field
      fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis)));
      AddDis("scatra1", fluidscatradis);

      // structure scatra field
      structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structscatradis)));
      AddDis("scatra2", structscatradis);

      break;
    }
    case ProblemType::fsi_xfem:
    case ProblemType::fluid_xfem:
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      AddDis("structure", structdis);
      meshreader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(), "GEOMETRY"),
          nullptr);

      if (DRT::INPUT::IntegralValue<int>(XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
        fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
        AddDis("fluid", fluiddis);

        xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
        xfluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis)));
        AddDis("xfluid", xfluiddis);

        meshreader.AddElementReader(Teuchos::rcp(
            new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID")));
      }
      else
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid", reader.Comm()));
        fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
        AddDis("fluid", fluiddis);

        meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
            DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(), "GEOMETRY"),
            nullptr);
      }

      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));
      AddDis("ale", aledis);
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));
      break;
    }
    case ProblemType::fpsi_xfem:
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("porofluid", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("structure", structdis);
      AddDis("porofluid", porofluiddis);
      AddDis("fluid", fluiddis);
      AddDis("ale", aledis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(), "GEOMETRY"),
          nullptr);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

      break;
    }
    case ProblemType::ale:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
          break;
        }
        default:
        {
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("ale", aledis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

      break;
    }
    case ProblemType::fluid:
    case ProblemType::fluid_redmodels:
    {
      if (distype == ShapeFunctionType::shapefunction_hdg)
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationHDG("fluid", reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      }
      else if (distype == ShapeFunctionType::shapefunction_nurbs)
      {
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));

        // create discretization writer - in constructor set ingto and owned by corresponding
        // discret
        fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      }
      else if (DRT::INPUT::IntegralValue<int>(FluidDynamicParams().sublist("WALL MODEL"), "X_WALL"))
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXWall("fluid", reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      }
      else
      {
        // fluiddis  = Teuchos::rcp(new DRT::Discretization("fluid",reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      }

      AddDis("fluid", fluiddis);

      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }
    case ProblemType::lubrication:
    {
      // create empty discretizations
      lubricationdis = Teuchos::rcp(new DRT::Discretization("lubrication", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      lubricationdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(lubricationdis)));

      AddDis("lubrication", lubricationdis);

      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS")));

      break;
    }
    case ProblemType::cardiac_monodomain:
    case ProblemType::scatra:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra", reader.Comm()));
          break;
        }
        case ShapeFunctionType::shapefunction_hdg:
        {
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::DiscretizationHDG("scatra", reader.Comm()));
          break;
        }
        default:
        {
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));


      AddDis("fluid", fluiddis);
      AddDis("scatra", scatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

      break;
    }
    case ProblemType::sti:
    {
      // safety checks
      if (distype == ShapeFunctionType::shapefunction_nurbs)
        dserror("Scatra-thermo interaction does not work for nurbs discretizations yet!");

      // create empty discretizations for scalar and thermo fields
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
      thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));

      // create discretization writers
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
      thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));

      // add empty discretizations to global problem
      AddDis("scatra", scatradis);
      AddDis("thermo", thermdis);

      // add element reader to node reader
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

      break;
    }
    case ProblemType::fluid_ale:
    case ProblemType::freesurf:
    {
      if (distype == ShapeFunctionType::shapefunction_hdg)
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationHDG("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }
      else if (distype == ShapeFunctionType::shapefunction_nurbs)
      {
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
      }
      else if (DRT::INPUT::IntegralValue<int>(FluidDynamicParams().sublist("WALL MODEL"), "X_WALL"))
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXWall("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }
      else
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
        if (DRT::INPUT::IntegralValue<bool>(
                XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
          xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }


      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      if (xfluiddis != Teuchos::null)
      {
        xfluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis)));
      }
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));


      AddDis("fluid", fluiddis);
      if (xfluiddis != Teuchos::null)
      {
        AddDis("xfluid", xfluiddis);  // xfem discretization on slot 1
      }
      AddDis("ale", aledis);

      if (xfluiddis != Teuchos::null)
      {
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS")));
      }
      else
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

      break;
    }
    case ProblemType::tsi:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          thermdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("thermo", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));

      AddDis("structure", structdis);
      AddDis("thermo", thermdis);

      meshreader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(), "GEOMETRY"),
          nullptr);
      meshreader.AddAdvancedReader(thermdis, reader, "THERMO",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(ThermalDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }
    case ProblemType::thermo:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          thermdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("thermo", reader.Comm()));
          break;
        }
        default:
        {
          thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));

      AddDis("thermo", thermdis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(thermdis, reader, "--THERMO ELEMENTS")));

      break;
    }

    case ProblemType::structure:
    case ProblemType::invana:
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));

      AddDis("structure", structdis);

      meshreader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }

    case ProblemType::polymernetwork:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      pboxdis = Teuchos::rcp(new DRT::Discretization("boundingbox", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      pboxdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(pboxdis)));

      AddDis("structure", structdis);
      AddDis("boundingbox", pboxdis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(pboxdis, reader, "--PERIODIC BOUNDINGBOX ELEMENTS")));

      break;
    }

    case ProblemType::loma:
    {
      // create empty discretizations
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      AddDis("fluid", fluiddis);
      AddDis("scatra", scatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

      break;
    }

    case ProblemType::two_phase_flow:
    case ProblemType::fluid_xfem_ls:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      if (GetProblemType() == ProblemType::fluid_xfem_ls)
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid", reader.Comm()));
      else
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);
      AddDis("scatra", scatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(), "GEOMETRY"),
          nullptr);
      // meshreader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader,
      // "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
      break;
    }

    case ProblemType::elch:
    {
      // create empty discretizations
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
          break;
        }
        default:
        {
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("fluid", fluiddis);
      AddDis("scatra", scatradis);
      AddDis("ale", aledis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

      break;
    }
    case ProblemType::art_net:  // _1D_ARTERY_
    {
      // create empty discretizations
      arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));

      // create empty discretizations
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          dserror("Nurbs Discretization not possible for artery");
          break;
        }
        default:
        {
          scatradis = Teuchos::rcp(new DRT::Discretization("artery_scatra", reader.Comm()));
          break;
        }
      }

      AddDis("artery", arterydis);
      AddDis("artery_scatra", scatradis);

      // create discretization writer - in constructor set into and owned by corresponding discret
      arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

      break;
    }
    case ProblemType::red_airways:  // _reduced D airways
    {
      // create empty discretizations
      airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));

      AddDis("red_airway", airwaydis);

      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));

      break;
    }
    case ProblemType::struct_ale:  // structure with ale
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("structure", structdis);
      AddDis("ale", aledis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

      break;
    }
    case ProblemType::fluid_topopt:
    {
      // create empty discretizations
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      optidis = Teuchos::rcp(new DRT::Discretization("opti", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      optidis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(optidis)));

      AddDis("fluid", fluiddis);
      AddDis("opti", optidis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

      break;
    }
    case ProblemType::poroelast:
    case ProblemType::poromultiphase:
    {
      // create empty discretizations
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          porofluiddis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));

      AddDis("structure", structdis);
      AddDis("porofluid", porofluiddis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));

      if (DRT::INPUT::IntegralValue<bool>(PoroMultiPhaseDynamicParams(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));
        AddDis("artery", arterydis);
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));
      }

      break;
    }
    case ProblemType::poromultiphasescatra:
    {
      // create empty discretizations
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          porofluiddis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      AddDis("structure", structdis);
      AddDis("porofluid", porofluiddis);
      AddDis("scatra", scatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

      if (DRT::INPUT::IntegralValue<bool>(PoroMultiPhaseScatraDynamicParams(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));
        AddDis("artery", arterydis);
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));

        artscatradis = Teuchos::rcp(new DRT::Discretization("artery_scatra", reader.Comm()));
        artscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(artscatradis)));
        AddDis("artery_scatra", artscatradis);
        meshreader.AddElementReader(Teuchos::rcp(
            new DRT::INPUT::ElementReader(artscatradis, reader, "--TRANSPORT ELEMENTS")));
      }

      break;
    }
    case ProblemType::porofluidmultiphase:
    {
      // create empty discretizations
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          porofluiddis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid", reader.Comm()));
          break;
        }
        default:
        {
          porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));

      AddDis("porofluid", porofluiddis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));

      if (DRT::INPUT::IntegralValue<bool>(PoroFluidMultiPhaseDynamicParams(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));
        AddDis("artery", arterydis);
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));
      }
      break;
    }
    case ProblemType::fpsi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("structure", structdis);
      AddDis("porofluid", porofluiddis);
      AddDis("fluid", fluiddis);
      AddDis("ale", aledis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

      break;
    }
    case ProblemType::fbi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }
    case ProblemType::immersed_fsi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

      break;
    }
    case ProblemType::fps3i:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

      AddDis("structure", structdis);
      AddDis("porofluid", porofluiddis);
      AddDis("fluid", fluiddis);
      AddDis("ale", aledis);


      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

      // fluid scatra field
      fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis)));
      AddDis("scatra1", fluidscatradis);

      // poro structure scatra field
      structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structscatradis)));
      AddDis("scatra2", structscatradis);

      break;
    }
    case ProblemType::poroscatra:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      AddDis("structure", structdis);
      AddDis("porofluid", porofluiddis);
      AddDis("scatra", scatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
      break;
    }
    case ProblemType::ehl:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      lubricationdis = Teuchos::rcp(new DRT::Discretization("lubrication", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      lubricationdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(lubricationdis)));

      AddDis("structure", structdis);
      AddDis("lubrication", lubricationdis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS")));

      break;
    }
    case ProblemType::ssi:
    case ProblemType::ssti:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      AddDis("structure", structdis);
      AddDis("scatra", scatradis);

      // consider case of additional scatra manifold
      if (DRT::INPUT::IntegralValue<bool>(SSIControlParams().sublist("MANIFOLD"), "ADD_MANIFOLD"))
      {
        auto scatra_manifold_dis =
            Teuchos::rcp(new DRT::Discretization("scatra_manifold", reader.Comm()));
        scatra_manifold_dis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(scatra_manifold_dis)));
        AddDis("scatra_manifold", scatra_manifold_dis);
      }

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

      if (GetProblemType() == ProblemType::ssti)
      {
        thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));
        thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));
        AddDis("thermo", thermdis);
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(thermdis, reader, "--TRANSPORT ELEMENTS")));
      }

      break;
    }
    case ProblemType::particle:
    case ProblemType::pasi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));

      AddDis("structure", structdis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

      break;
    }
    case ProblemType::level_set:
    {
      // create empty discretizations
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

      AddDis("scatra", scatradis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
      break;
    }
    case ProblemType::np_support:
    {
      // no discretizations and nodes needed for supporting procs
      break;
    }
    case ProblemType::elemag:
    {
      // create empty discretizations
      elemagdis = Teuchos::rcp(new DRT::DiscretizationHDG("elemag", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      elemagdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(elemagdis)));

      AddDis("elemag", elemagdis);

      std::set<std::string> elemagelementtypes;
      elemagelementtypes.insert("ELECTROMAGNETIC");
      elemagelementtypes.insert("ELECTROMAGNETICDIFF");

      meshreader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(
          elemagdis, reader, "--ELECTROMAGNETIC ELEMENTS", elemagelementtypes)));

      break;
    }
    case ProblemType::redairways_tissue:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));

      AddDis("structure", structdis);
      AddDis("red_airway", airwaydis);

      meshreader.AddElementReader(
          Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      meshreader.AddElementReader(Teuchos::rcp(
          new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));
    }
    break;
    case ProblemType::tutorial:
    {
    }
    break;
    default:
      dserror("Unknown problem type: %d", GetProblemType());
      break;
  }

  // add artery or airways discretizations only for the following problem types
  switch (GetProblemType())
  {
    case ProblemType::fsi_redmodels:
    case ProblemType::fsi_lung:
    case ProblemType::fluid_ale:
    case ProblemType::fluid_redmodels:
    {
      if (distype == ShapeFunctionType::shapefunction_polynomial)
      {
        // create empty discretizations
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        // create discretization writer - in constructor set into and owned by corresponding discret
        arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));
        AddDis("artery", arterydis);
        meshreader.AddElementReader(
            Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));

        airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway", reader.Comm()));
        // create discretization writer - in constructor set into and owned by corresponding discret
        airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));
        AddDis("red_airway", airwaydis);
        meshreader.AddElementReader(Teuchos::rcp(
            new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));
      }
    }
    break;
    default:
      break;
  }

  if (readmesh)  // now read and allocate!
  {
    // we read nodes and elements for the desired fields as specified above
    meshreader.ReadAndPartition();

    NestedParallelismType npType = DRT::Problem::Instance()->GetCommunicators()->NpType();
    // care for special applications
    switch (GetProblemType())
    {
      case ProblemType::elch:
      case ProblemType::fsi:
      case ProblemType::fsi_redmodels:
      case ProblemType::fsi_lung:
      case ProblemType::invana:
      case ProblemType::scatra:
      case ProblemType::structure:
      {
        // read microscale fields from second, third, ... input file if necessary
        // (in case of multi-scale material models)
        if (npType != copy_dat_file) ReadMicroFields(reader);
        break;
      }
      case ProblemType::np_support:
      {
        // read microscale fields from second, third, ... inputfile for supporting processors
        ReadMicrofieldsNPsupport();
        break;
      }
      case ProblemType::tutorial:
      {
        break;
      }
      default:
        break;
    }
  }  // if(readmesh)
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMicroFields(DRT::INPUT::DatFileReader& reader)
{
  // check whether micro material is specified
  const int id_struct =
      DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_struct_multiscale);
  const int id_scatra =
      DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_scatra_multiscale);
  const int id_elch =
      DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_newman_multiscale);

  // return if no multiscale material is used
  if (id_struct == -1 and id_scatra == -1 and id_elch == -1) return;

  // safety check
  if ((id_struct != -1 and id_scatra != -1) or (id_struct != -1 and id_elch != -1) or
      (id_scatra != -1 and id_elch != -1))
    dserror("Cannot have more than one multi-scale material!");

  // store name of macro-scale discretization in string
  std::string macro_dis_name("");
  if (id_struct != -1)
    macro_dis_name = "structure";
  else
    macro_dis_name = "scatra";

  // fetch communicators
  Teuchos::RCP<Epetra_Comm> lcomm = communicators_->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = communicators_->GlobalComm();

  DRT::Problem* macro_problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> macro_dis = macro_problem->GetDis(macro_dis_name);

  // repartition macro problem for a good distribution of elements with micro material
  if (macro_dis_name == "structure")
    DRT::UTILS::REBALANCING::RedistributeAndFillCompleteDiscretizationUsingWeights(
        macro_dis, true, true, true);

  // make sure that we read the micro discretizations only on the processors on
  // which elements with the corresponding micro material are evaluated

  std::set<int> my_multimat_IDs;

  // take care also of ghosted elements! -> ElementColMap!
  for (int i = 0; i < macro_dis->ElementColMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = macro_dis->lColElement(i);
    Teuchos::RCP<MAT::Material> actmat = actele->Material();

    if (id_elch != -1 and actmat->MaterialType() == INPAR::MAT::m_elchmat)
    {
      // extract wrapped material
      auto elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(actmat);
      auto elchphase =
          Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
      actmat = elchphase->MatById(elchphase->MatID(0));
    }

    if ((actmat->MaterialType() == INPAR::MAT::m_struct_multiscale and
            macro_dis_name == "structure") or
        (actmat->MaterialType() == INPAR::MAT::m_scatra_multiscale and
            macro_dis_name == "scatra") or
        (actmat->MaterialType() == INPAR::MAT::m_newman_multiscale and macro_dis_name == "scatra"))
    {
      MAT::PAR::Parameter* actparams = actmat->Parameter();
      my_multimat_IDs.insert(actparams->Id());
    }
  }

  // check which macro procs have an element with micro material
  int foundmicromat = 0;
  int foundmicromatmyrank = -1;
  if (my_multimat_IDs.size() != 0)
  {
    foundmicromat = 1;
    foundmicromatmyrank = lcomm->MyPID();
  }

  // find out how many procs have micro material
  int nummicromat = 0;
  lcomm->SumAll(&foundmicromat, &nummicromat, 1);
  // broadcast number of procs that have micro material
  gcomm->Broadcast(&nummicromat, 1, 0);

  // every proc needs to know which procs have micro material in order to distribute colors
  // array is filled with either its local proc id or -1 when no micro mat was found
  std::vector<int> foundmyranks;
  foundmyranks.resize(lcomm->NumProc(), -1);
  lcomm->GatherAll(&foundmicromatmyrank, &foundmyranks[0], 1);

  // determine color of macro procs with any contribution to micro material, only important for
  // procs with micro material color starts with 0 and is incremented for each group
  int color = -1;
  if (foundmicromat == 1)
  {
    for (int foundmyrank : foundmyranks)
    {
      if (foundmyrank != -1) ++color;
      if (foundmyrank == foundmicromatmyrank) break;
    }
  }
  else
  {
    color = MPI_UNDEFINED;
  }

  // do the splitting of the communicator (macro proc must always be proc in subcomm with lowest key
  // --> 0 is inserted here)
  MPI_Comm mpi_local_comm;
  MPI_Comm_split((Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm, true)->GetMpiComm()), color,
      0 /*important here*/, &mpi_local_comm);

  // sort out macro procs that do not have micro material
  if (foundmicromat == 1)
  {
    // create the sub communicator that includes one macro proc and some supporting procs
    Teuchos::RCP<Epetra_Comm> subgroupcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
    communicators_->SetSubComm(subgroupcomm);

    // find out how many micro problems have to be solved on this macro proc
    int microcount = 0;
    for (const auto& material_map : *materials_->Map())
    {
      int matid = material_map.first;
      if (my_multimat_IDs.find(matid) != my_multimat_IDs.end()) microcount++;
    }
    // and broadcast it to the corresponding group of procs
    subgroupcomm->Broadcast(&microcount, 1, 0);

    for (const auto& material_map : *materials_->Map())
    {
      int matid = material_map.first;

      if (my_multimat_IDs.find(matid) != my_multimat_IDs.end())
      {
        Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);

        // initialize variables storing micro-scale information
        int microdisnum(-1);
        std::string micro_dis_name = "";
        std::string micro_inputfile_name("");
        DRT::Problem* micro_problem(nullptr);

        // structure case
        if (macro_dis_name == "structure")
        {
          // access multi-scale structure material
          auto* micromat = static_cast<MAT::MicroMaterial*>(mat.get());

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->MicroDisNum();
          subgroupcomm->Broadcast(&microdisnum, 1, 0);

          // set name of micro-scale discretization
          micro_dis_name = "structure";

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->MicroInputFileName();

          // instantiate micro-scale problem
          micro_problem = DRT::Problem::Instance(microdisnum);
        }

        // scalar transport case
        else
        {
          // access multi-scale scalar transport material
          MAT::ScatraMultiScale* micromat = nullptr;
          if (id_scatra != -1)
            micromat = dynamic_cast<MAT::ScatraMatMultiScale*>(mat.get());
          else if (id_elch != -1)
            micromat = dynamic_cast<MAT::NewmanMultiScale*>(mat.get());
          else
            dserror("How the heck did you get here?!");

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->MicroDisNum();
          subgroupcomm->Broadcast(&microdisnum, 1, 0);

          // set unique name of micro-scale discretization
          std::stringstream name;
          name << "scatra_multiscale_" << microdisnum;
          micro_dis_name = name.str();

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->MicroInputFileName();

          // instantiate micro-scale problem
          micro_problem = DRT::Problem::Instance(microdisnum);
        }

        if (micro_inputfile_name[0] != '/')
        {
          std::string filename = reader.MyInputfileName();
          std::string::size_type pos = filename.rfind('/');
          if (pos != std::string::npos)
          {
            std::string path = filename.substr(0, pos + 1);
            micro_inputfile_name.insert(micro_inputfile_name.begin(), path.begin(), path.end());
          }
        }

        // broadcast micro input file name
        int length = static_cast<int>(micro_inputfile_name.length());
        subgroupcomm->Broadcast(&length, 1, 0);
        subgroupcomm->Broadcast((const_cast<char*>(micro_inputfile_name.c_str())), length, 0);

        // start with actual reading
        DRT::INPUT::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

        Teuchos::RCP<DRT::Discretization> dis_micro =
            Teuchos::rcp(new DRT::Discretization(micro_dis_name, micro_reader.Comm()));

        // replace standard dofset inside micro discretization by independent dofset
        // to avoid inconsistent dof numbering in non-nested parallel settings with more than one
        // micro discretization
        if (communicators_->NpType() == no_nested_parallelism)
          dis_micro->ReplaceDofSet(Teuchos::rcp(new DRT::IndependentDofSet()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        dis_micro->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(dis_micro)));

        micro_problem->AddDis(micro_dis_name, dis_micro);

        micro_problem->ReadParameter(micro_reader);

        // read materials of microscale
        // CAUTION: materials for microscale cannot be read until
        // micro_reader is activated, since else materials will again be
        // read from macroscale inputfile. Besides, materials MUST be read
        // before elements are read since elements establish a connection
        // to the corresponding material! Thus do not change position of
        // function calls!
        materials_->SetReadFromProblem(microdisnum);

        micro_problem->ReadMaterials(micro_reader);

        DRT::INPUT::MeshReader micromeshreader(micro_reader, "--NODE COORDS");

        if (micro_dis_name == "structure")
        {
          micromeshreader.AddElementReader(Teuchos::rcp(
              new DRT::INPUT::ElementReader(dis_micro, micro_reader, "--STRUCTURE ELEMENTS")));
        }
        else
          micromeshreader.AddElementReader(Teuchos::rcp(
              new DRT::INPUT::ElementReader(dis_micro, micro_reader, "--TRANSPORT ELEMENTS")));

        micromeshreader.ReadAndPartition();

        // read conditions of microscale
        // -> note that no time curves and spatial functions can be read!
        micro_problem->ReadTimeFunctionResult(micro_reader);
        micro_problem->ReadConditions(micro_reader);

        // At this point, everything for the microscale is read,
        // subsequent reading is only for macroscale
        dis_micro->FillComplete();

        // broadcast restart information
        subgroupcomm->Broadcast(&restartstep_, 1, 0);

        // set the problem number from which to call materials again to zero
        // (i.e. macro problem), cf. MAT::Material::Factory!
        materials_->ResetReadFromProblem();
      }
    }
    materials_->ResetReadFromProblem();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMicrofieldsNPsupport()
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetCommunicators()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetCommunicators()->GlobalComm();

  // receive number of procs that have micro material
  int nummicromat = 0;
  gcomm->Broadcast(&nummicromat, 1, 0);

  // prepare the supporting procs for a splitting of gcomm

  // groups should be equally sized
  // in a first step every macro proc that needs support gets procpergroup supporting procs
  int procpergroup = int(floor((lcomm->NumProc()) / nummicromat));
  std::vector<int> supgrouplayout(nummicromat, procpergroup);
  // remaining procs are added to the groups in the beginning
  int remainingProcs = lcomm->NumProc() - procpergroup * nummicromat;
  for (int k = 0; k < remainingProcs; ++k)
  {
    supgrouplayout[k]++;
  }

  // secondly: colors are distributed
  // color starts with 0 and is incremented for each group
  int color = -1;
  int gsum = 0;
  do
  {
    color++;
    gsum += supgrouplayout[color];
  } while (gsum <= lcomm->MyPID());

  // do the splitting of the communicator
  MPI_Comm mpi_local_comm;
  MPI_Comm_split((Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm, true)->GetMpiComm()), color,
      gcomm->MyPID(), &mpi_local_comm);

  // create the sub communicator that includes one macro proc and some supporting procs
  Teuchos::RCP<Epetra_Comm> subgroupcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
  communicators_->SetSubComm(subgroupcomm);

  // number of micro problems for this sub group
  int microcount = 0;
  subgroupcomm->Broadcast(&microcount, 1, 0);

  for (int n = 0; n < microcount; n++)
  {
    // broadcast microdis number
    int microdisnum = -1;
    subgroupcomm->Broadcast(&microdisnum, 1, 0);

    DRT::Problem* micro_problem = DRT::Problem::Instance(microdisnum);

    // broadcast micro input file name
    int length = -1;
    std::string micro_inputfile_name;
    subgroupcomm->Broadcast(&length, 1, 0);
    micro_inputfile_name.resize(length);
    subgroupcomm->Broadcast((const_cast<char*>(micro_inputfile_name.c_str())), length, 0);

    // start with actual reading
    DRT::INPUT::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

    Teuchos::RCP<DRT::Discretization> structdis_micro =
        Teuchos::rcp(new DRT::Discretization("structure", micro_reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis_micro->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis_micro)));

    micro_problem->AddDis("structure", structdis_micro);

    micro_problem->ReadParameter(micro_reader);

    // read materials of microscale
    // CAUTION: materials for microscale cannot be read until
    // micro_reader is activated, since else materials will again be
    // read from macroscale inputfile. Besides, materials MUST be read
    // before elements are read since elements establish a connection
    // to the corresponding material! Thus do not change position of
    // function calls!
    materials_->SetReadFromProblem(microdisnum);

    micro_problem->ReadMaterials(micro_reader);

    DRT::INPUT::MeshReader micromeshreader(micro_reader, "--NODE COORDS");
    micromeshreader.AddElementReader(Teuchos::rcp(
        new DRT::INPUT::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS")));
    micromeshreader.ReadAndPartition();

    // read conditions of microscale
    // -> note that no time curves and spatial functions can be read!

    micro_problem->ReadConditions(micro_reader);

    // At this point, everything for the microscale is read,
    // subsequent reading is only for macroscale
    structdis_micro->FillComplete();

    // broadcast restart information
    subgroupcomm->Broadcast(&restartstep_, 1, 0);

    // set the problem number from which to call materials again to zero
    // (i.e. macro problem), cf. MAT::Material::Factory!
    materials_->ResetReadFromProblem();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  try
  {
    // Test parameter list against valid parameters, set default values
    // and set validator objects to extract numerical values for string
    // parameters.
    paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  }
  catch (Teuchos::Exceptions::InvalidParameter& err)
  {
    std::cerr << "\n\n" << err.what();
    dserror("Input parameter validation failed. Fix your input file.");
  }

  // yes, it is my list
  setMyParamList(paramList);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> DRT::Problem::getValidParameters() const
{
  // call the external method to get the valid parameters
  // this way the parameter configuration is separate from the source
  return DRT::INPUT::ValidParameters();
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
  // For debug: what's currently in the map:
  /*
  std::map<std::string,Teuchos::RCP<Discretization> >::iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
    std::cout << "key : " << iter->first << "    " << "discret.name = " << iter->second->Name() <<
  std::endl << std::endl;
  }
  */
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
void DRT::Problem::SetRestartStep(int r)
{
  if (r < 0) dserror("Negative restart step not allowed");

  restartstep_ = r;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetProblemType(ProblemType targettype) { probtype_ = targettype; }
