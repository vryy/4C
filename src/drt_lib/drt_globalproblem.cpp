 /*----------------------------------------------------------------------*/
/*!
\file drt_globalproblem.cpp

\brief global list of problems

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*/
/*----------------------------------------------------------------------*/

#include <sys/time.h>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_ParameterListExceptions.hpp>

#include <Epetra_Time.h>
#include <Epetra_Comm.h>

#include "drt_conditiondefinition.H"
#include "drt_materialdefinition.H"
#include "drt_function.H"
#include "drt_singletondestruction.H"
#include "drt_globalproblem.H"
#include "drt_inputreader.H"
#include "drt_elementreader.H"
#include "drt_nodereader.H"
#include "drt_timecurve.H"
#include "drt_utils_parallel.H"
#include "drt_utils_createdis.H"
#include "drt_discret.H"
#include "drt_discret_faces.H"
#include "drt_discret_xwall.H"
#include "drt_discret_hdg.H"
#include "drt_discret_combust.H"
#include "drt_discret_xfem.H"
#include "drt_linedefinition.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/drt_validmaterials.H"
#include "../drt_inpar/inpar.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/scatra_mat_multiscale.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_meshfree_discret/drt_meshfree_discret.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_problemtype.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_mlmc.H"

/*----------------------------------------------------------------------*/
// the instances
/*----------------------------------------------------------------------*/
std::vector<DRT::Problem* > DRT::Problem::instances_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Problem* DRT::Problem::Instance(int num)
{
  if (num > static_cast<int>(instances_.size())-1)
  {
    instances_.resize(num+1);
    instances_[num] = new Problem();
  }
  return instances_[num];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::Done()
{
  // destroy singleton objects when the problem object is still alive
  for ( std::vector<Problem* >::iterator i=instances_.begin(); i!=instances_.end(); ++i )
  {
    Problem * p = *i;
    for (std::vector<DRT::SingletonDestruction *>::iterator j=p->sds_.begin(); j!=p->sds_.end(); ++j)
    {
      DRT::SingletonDestruction * sd = *j;
      sd->Done();
    }
    p->sds_.clear();
  }

  // This is called at the very end of a baci run.
  //
  // It removes all global problem objects. Therefore all
  // discretizations as well and everything inside those.
  //
  // There is a whole lot going on here...
  for ( std::vector<Problem* >::iterator i=instances_.begin(); i!=instances_.end(); ++i )
  {
    delete *i;
    *i = 0;
  }
  instances_.clear();

  // close the parallel output environment to make sure all files are properly closed
  IO::cout.close();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Problem::Problem() :
  probtype_(prb_none),
  restartstep_(0),
  restarttime_(0.0),
  npgroup_(Teuchos::null)
{
  materials_ = Teuchos::rcp(new MAT::PAR::Bundle());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::Problem::~Problem()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PROBLEM_TYP DRT::Problem::ProblemType() const
{
  return probtype_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string DRT::Problem::ProblemName() const
{
  std::map<std::string,PROBLEM_TYP> map = DRT::StringToProblemTypeMap();
  std::map<std::string,PROBLEM_TYP>::const_iterator i;

  for (i = map.begin(); i != map.end();++i)
  {
    if (i->second == probtype_)
      return i->first;
  }
  dserror("Could not determine valid problem name");
  return "Undefined";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::Problem::Restart() const
{
  return restartstep_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::Problem::RestartTime() const
{
  return restarttime_;
}


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
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6);
}


#if 0 // Currently unused, might come back into usage though
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::Problem::BandWidthOpt() const
{
  // in case of filters (e.g. post_drt_ensight) we don't have the list
  Teuchos::RCP<const Teuchos::ParameterList> list = getParameterList();
  if (list==Teuchos::null) return false;

  const Teuchos::ParameterList& typeparams = ProblemTypeParams();
  bool yesno = Teuchos::getIntegralValue<int>(typeparams,"BANDWIDTHOPT");
  return yesno;
}
#endif

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string DRT::Problem::SpatialApproximation() const
{
  // TODO: fix downstream use in fluid_ele.cpp
  // in case we do not have the global list available (post-processing etc)
  if (getParameterList() == Teuchos::null ||
      (not getParameterList()->isSublist("PROBLEM TYP")))
    return "Polynomial";

  // decide which kind of spatial representation is required
  const Teuchos::ParameterList& ptype = ProblemTypeParams();

  std::string basis_fct_type = ptype.get<std::string>("SHAPEFCT");

  return basis_fct_type;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadParameter(DRT::INPUT::DatFileReader& reader)
{
  Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList("DAT FILE"));

  reader.ReadGidSection("--DISCRETISATION", *list);
  reader.ReadGidSection("--PROBLEM SIZE", *list);
  reader.ReadGidSection("--PROBLEM TYP", *list);
  reader.ReadGidSection("--MESHFREE", *list);
  reader.ReadGidSection("--IO", *list);
  reader.ReadGidSection("--DESIGN DESCRIPTION", *list);
  reader.ReadGidSection("--PATIENT SPECIFIC", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--STRUCTURAL DYNAMIC/GEMM", *list);
  reader.ReadGidSection("--INVERSE ANALYSIS", *list);
  reader.ReadGidSection("--STAT INVERSE ANALYSIS", *list);
  reader.ReadGidSection("--MULTI LEVEL MONTE CARLO", *list);
  reader.ReadGidSection("--MORTAR COUPLING", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/AUGMENTED/STEEPESTASCENT", *list);
  reader.ReadGidSection("--CONTACT DYNAMIC/XCONTACT", *list);
  reader.ReadGidSection("--XCONTACT DYNAMIC", *list);
  reader.ReadGidSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING", *list);
  reader.ReadGidSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING/SYS-PUL CIRCULATION PARAMETERS", *list);
  reader.ReadGidSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING/RESPIRATORY PARAMETERS", *list);
  reader.ReadGidSection("--FLUCTUATING HYDRODYNAMICS", *list);
  reader.ReadGidSection("--BROWNIAN DYNAMICS", *list);
  reader.ReadGidSection("--BEAM INTERACTION", *list);
  reader.ReadGidSection("--BEAM INTERACTION/CONTRACTILE CELLS", *list);
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
  reader.ReadGidSection("--POROMULTIPHASE DYNAMIC", *list);
  reader.ReadGidSection("--POROMULTIPHASESCATRA DYNAMIC", *list);
  reader.ReadGidSection("--ELASTO HYDRO DYNAMIC", *list);
  reader.ReadGidSection("--ELASTO HYDRO DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--ELASTO HYDRO DYNAMIC/MONOLITHIC", *list);
  reader.ReadGidSection("--SSI CONTROL", *list);
  reader.ReadGidSection("--SSI CONTROL/PARTITIONED", *list);
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
  reader.ReadGidSection("--TWO PHASE FLOW", *list);
  reader.ReadGidSection("--TWO PHASE FLOW/SURFACE TENSION", *list);
  reader.ReadGidSection("--TWO PHASE FLOW/SMEARED", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION FLUID", *list);
  reader.ReadGidSection("--COMBUSTION CONTROL/COMBUSTION GFUNCTION", *list);
  reader.ReadGidSection("--LUBRICATION DYNAMIC", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/NONLINEAR", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--SCALAR TRANSPORT DYNAMIC/S2I COUPLING", *list);
  reader.ReadGidSection("--STI DYNAMIC", *list);
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
  reader.ReadGidSection("--IMMERSED METHOD", *list);
  reader.ReadGidSection("--IMMERSED METHOD/PARTITIONED SOLVER", *list);
  reader.ReadGidSection("--CELL DYNAMIC", *list);
  reader.ReadGidSection("--CELL DYNAMIC/SCALAR TRANSPORT DOF IDS", *list);
  reader.ReadGidSection("--CELL DYNAMIC/STRUCTURAL DYNAMIC", *list);
  reader.ReadGidSection("--CELL DYNAMIC/STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadGidSection("--CELL DYNAMIC/STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadGidSection("--CELL DYNAMIC/SCALAR TRANSPORT", *list);
  reader.ReadGidSection("--CELL DYNAMIC/SCALAR TRANSPORT/NONLINEAR", *list);
  reader.ReadGidSection("--CELL DYNAMIC/SCALAR TRANSPORT/STABILIZATION", *list);
  reader.ReadGidSection("--CELL DYNAMIC/PROTEOLYSIS MODULE", *list);
  reader.ReadGidSection("--CELL DYNAMIC/CONFINEMENT MODULE", *list);
  reader.ReadGidSection("--CELL DYNAMIC/ADHESION MODULE", *list);
  reader.ReadGidSection("--CELL DYNAMIC/PROTRUSION MODULE", *list);
  reader.ReadGidSection("--CELL DYNAMIC/FLOW INTERACTION MODULE", *list);
  reader.ReadGidSection("--CELL DYNAMIC/ENDOEXOCYTOSIS MODULE", *list);
  reader.ReadGidSection("--FPSI DYNAMIC", *list);
  reader.ReadGidSection("--ARTERIAL DYNAMIC", *list);
  reader.ReadGidSection("--REDUCED DIMENSIONAL AIRWAYS DYNAMIC", *list);
  reader.ReadGidSection("--COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC", *list);
  reader.ReadGidSection("--SEARCH TREE", *list);
  reader.ReadGidSection("--XFEM GENERAL", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC/GENERAL", *list);
  reader.ReadGidSection("--XFLUID DYNAMIC/STABILIZATION", *list);
  reader.ReadGidSection("--LOMA CONTROL", *list);
  reader.ReadGidSection("--ELCH CONTROL", *list);
  reader.ReadGidSection("--ELCH CONTROL/DIFFCOND", *list);
  reader.ReadGidSection("--BIOFILM CONTROL", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL/TOPOLOGY OPTIMIZER", *list);
  reader.ReadGidSection("--TOPOLOGY OPTIMIZATION CONTROL/TOPOLOGY ADJOINT FLUID", *list);
  reader.ReadGidSection("--CAVITATION DYNAMIC", *list);
  reader.ReadGidSection("--PARTICLE DYNAMIC", *list);
  reader.ReadGidSection("--PARTICLE DYNAMIC/GENALPHA", *list);
  reader.ReadGidSection("--PASI DYNAMIC", *list);
  reader.ReadGidSection("--PASI DYNAMIC/PARTITIONED", *list);
  reader.ReadGidSection("--LEVEL-SET CONTROL", *list);
  reader.ReadGidSection("--LEVEL-SET CONTROL/PARTICLE", *list);
  reader.ReadGidSection("--LEVEL-SET CONTROL/REINITIALIZATION", *list);
  reader.ReadGidSection("--WEAR", *list);
  reader.ReadGidSection("--BEAM CONTACT", *list);
  reader.ReadGidSection("--BEAM POTENTIAL", *list);
  reader.ReadGidSection("--SEMI-SMOOTH PLASTICITY", *list);
  reader.ReadGidSection("--ACOUSTIC DYNAMIC", *list);
  reader.ReadGidSection("--ACOUSTIC DYNAMIC/PA IMAGE RECONSTRUCTION", *list);
  reader.ReadGidSection("--VOLMORTAR COUPLING", *list);
  reader.ReadGidSection("--NONLINEAR SOLVER", *list);
  reader.ReadGidSection("--TUTORIAL DYNAMIC", *list);
  reader.ReadGidSection("--TUTORIAL DYNAMIC/NONLINEAR TRUSS", *list);
  reader.ReadGidSection("--TUTORIAL DYNAMIC/FIXED POINT SCHEME", *list);
  reader.ReadGidSection("--CARDIAC MONODOMAIN CONTROL", *list);

  reader.ReadSection("--STRUCT NOX", *list);
  reader.ReadSection("--STRUCT NOX/Direction", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton", *list);
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
  // drt_validparameters.cpp, too!
  for (int i = 1; i<10; i++) {
    std::stringstream ss;
    ss << "--SOLVER " << i;
    reader.ReadGidSection(ss.str(), *list);

    // adapt path of XML file if necessary
    Teuchos::ParameterList& sublist = list->sublist(ss.str().substr(2));

    std::string* stratimikos_xmlfile = sublist.getPtr<std::string>("STRATIMIKOS_XMLFILE");
    if (stratimikos_xmlfile != NULL and *stratimikos_xmlfile != "none")
    {
      // make path relative to input file path if it is not an absolute path
      if ((*stratimikos_xmlfile)[0]!='/')
      {
        std::string filename = reader.MyInputfileName();
        std::string::size_type pos = filename.rfind('/');
        if (pos!=std::string::npos)
        {
          std::string tmp = filename.substr(0,pos+1);
          stratimikos_xmlfile->insert(stratimikos_xmlfile->begin(), tmp.begin(), tmp.end());
        }
      }
    }

    std::string* amgnxn_xmlfile = sublist.getPtr<std::string>("AMGNXN_XML_FILE");
    if (amgnxn_xmlfile != NULL and *amgnxn_xmlfile != "none")
    {
      // make path relative to input file path if it is not an absolute path
      if ((*amgnxn_xmlfile)[0]!='/')
      {
        std::string filename = reader.MyInputfileName();
        std::string::size_type pos = filename.rfind('/');
        if (pos!=std::string::npos)
        {
          std::string tmp = filename.substr(0,pos+1);
          amgnxn_xmlfile->insert(amgnxn_xmlfile->begin(), tmp.begin(), tmp.end());
        }
      }
    }
  }

  reader.ReadGidSection("--UMFPACK SOLVER",*list);

  // read in random field sections
  // Note: the maximum number of random fields in dat files is hardwired here.
  // If you change this do not forget to edit the corresponding parts in
  // drt_validparameters.cpp, too!
  for (int i = 1; i<4; i++) {
    std::stringstream ss;
    ss << "--RANDOM FIELD " << i;
    reader.ReadGidSection(ss.str(), *list);
  }

  // read in random variable sections
  // Note: the maximum number of random variable in dat files is hardwired here.
  // If you change this do not forget to edit the corresponding parts in
  // drt_validparameters.cpp, too!
  for (int i = 1; i<4; i++) {
    std::stringstream ss;
    ss << "--RANDOM VARIABLE " << i;
    reader.ReadGidSection(ss.str(), *list);
  }

  // read in STRUCT NOX/Status Test and modify the xml file name, if there
  // is one.
  if (list->sublist("STRUCT NOX").sublist("Status Test").isParameter("XML File"))
  {
    // adapt path of XML file if necessary
    Teuchos::ParameterList& sublist = list->sublist("STRUCT NOX").sublist("Status Test");
    std::string* statustest_xmlfile = sublist.getPtr<std::string>("XML File");
    // make path relative to input file path if it is not an absolute path
    if (((*statustest_xmlfile)[0]!='/') and ((*statustest_xmlfile) != "none"))
    {
      std::string filename = reader.MyInputfileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos!=std::string::npos)
      {
        std::string tmp = filename.substr(0,pos+1);
        statustest_xmlfile->insert(statustest_xmlfile->begin(), tmp.begin(), tmp.end());
      }
    }
  } // STRUCT NOX/Status Test

  // check for invalid parameters
  setParameterList(list);

  //---------------------------------------------------------------------
  // Now we have successfully read the whole input file. It's time to access some data

  // 1) get the problem type
  const Teuchos::ParameterList& type = ProblemTypeParams();
  probtype_ = DRT::INPUT::IntegralValue<PROBLEM_TYP>(type,"PROBLEMTYP");

  // 2) do the restart business with the four options we support (partially)
  if (restartstep_==0 )
  {
    // no restart flag on the command line, so check the restart flag from the input file
    restartstep_ = type.get<int>("RESTART");
  }
  else // SetRestartStep() has been called before!
  {
    // There is a non-zero restart flag on the command line, so we ignore the input file.
    // The RESTART flag in the input file should be zero or have the same value!
    const int restartflaginfile = type.get<int>("RESTART");
    if ((restartflaginfile > 0) and (restartflaginfile != restartstep_))
      dserror("Restart flags in input file and command line are non-zero and different!");
  }

  // Set restart time
  restarttime_ = type.get<double>("RESTARTTIME");
  if (restarttime_>0.0)
  {
    // Currently this option is implemented only for scalar structure interaction problems
    if (ProblemType() != prb_ssi)
      dserror("Restart with time option currently only implemented for SSI problems");
    // The value restartstep_ is used very deep down in Baci. Therefore we demand the user
    // to set this value, if one wants to use restart time option.
    // If this feature should be expanded for all problemtyps another handling of the
    // restartstep_ has to be considered
    if (restartstep_==0)
      dserror("Restart with time option needs a RESTART flag different from 0");
  }

  // Set restart time based on walltime
  const double restartinterval = IOParams().get<double>("RESTARTWALLTIMEINTERVAL");
  const int restartevry = IOParams().get<int>("RESTARTEVRY");
  RestartManager()->SetupRestartManager(restartinterval, restartevry);

  // 3) set random seed
  // time is in seconds, therefore we add the global processor id to obtain a unique seed on each proc
  {
    int rs = type.get<int>("RANDSEED");
    if (rs < 0)
      rs = (int)time(NULL) + 42*DRT::Problem::Instance(0)->GetNPGroup()->GlobalComm()->MyPID();

    srand( (unsigned int)rs ); // Set random seed for stdlibrary. This is deprecated, as it does not produce random numbers on some platforms!
    Random()->SetRandSeed( (unsigned int)rs ); // Use this instead.
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
const Teuchos::ParameterList& DRT::Problem::RandomFieldParams(
    int randomfieldNr) const
{
  std::stringstream ss;
  ss << "RANDOM FIELD " << randomfieldNr;
  return getParameterList()->sublist(ss.str());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& DRT::Problem::RandomVariableParams(
    int randomvariableNr) const
{
  std::stringstream ss;
  ss << "RANDOM VARIABLE " << randomvariableNr;
  return getParameterList()->sublist(ss.str());
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::ParameterList& DRT::Problem::UMFPACKSolverParams()
{
  Teuchos::RCP<Teuchos::ParameterList> params = getNonconstParameterList();

  Teuchos::ParameterList& subParams = params->sublist("UMFPACK SOLVER");
  subParams.set("SOLVER", "UMFPACK");
  subParams.set("NAME"  , "temporary UMFPACK solver");

  return getParameterList()->sublist("UMFPACK SOLVER");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::NPGroup(
  int groupId,
  int ngroup,
  std::map<int, int> lpidgpid,
  Teuchos::RCP<Epetra_Comm> lcomm,
  Teuchos::RCP<Epetra_Comm> gcomm,
  NP_TYPE npType
  )
{
  if (npgroup_ != Teuchos::null) dserror("NPGroup was already set.");
  npgroup_ = Teuchos::rcp(new COMM_UTILS::NestedParGroup(groupId,
                                                         ngroup,
                                                         lpidgpid,
                                                         lcomm,
                                                         gcomm,
                                                         npType));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::NPGroup(
    COMM_UTILS::NestedParGroup& npgroup
  )
{
  if (npgroup_ != Teuchos::null) dserror("NPGroup was already set.");
  npgroup_ = Teuchos::rcp(new COMM_UTILS::NestedParGroup(npgroup));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<COMM_UTILS::NestedParGroup> DRT::Problem::GetNPGroup()
{
  if (npgroup_ == Teuchos::null) dserror("No NPGroup allocated yet.");
  return npgroup_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMaterials(DRT::INPUT::DatFileReader& reader)
{
  // create list of known materials
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> > > vm = DRT::INPUT::ValidMaterials();
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >& matlist = *vm;

  // test for each material definition (input file --MATERIALS section)
  // and store in #matmap_
  for (unsigned m=0; m<matlist.size(); ++m)
  {
    // read material from DAT file of type #matlist[m]
    matlist[m]->Read(*this,reader,materials_);
  }

  // check if every material was identified
  const std::string name = "--MATERIALS";
  std::vector<const char*> section = reader.Section(name);
  int nummat = 0;
  if (section.size() > 0)
  {
    for (std::vector<const char*>::iterator i=section.begin();
         i!=section.end();
         ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      std::string mat;
      std::string number;
      std::string name;
      (*condline) >> mat >> number >> name;
      if ( (not (*condline)) or (mat != "MAT") )
        dserror("invalid material line in '%s'",name.c_str());

      // extract material ID
      int matid = -1;
      {
        char* ptr;
        matid = strtol(number.c_str(),&ptr,10);
        if (ptr == number.c_str())
          dserror("failed to read material object number '%s'",
                  number.c_str());
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


  // inform user
  //std::cout << "Number of successfully read materials is " << nummat << std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadCloningMaterialMap(DRT::INPUT::DatFileReader& reader)
{
  Teuchos::RCP<DRT::INPUT::Lines> lines = DRT::UTILS::ValidCloningMaterialMapLines();

  // perform the actual reading and extract the input parameters
  std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> > input = lines->Read(reader);
  for (size_t i =0 ; i < input.size(); i++)
  {
    // extract what was read from the input file
    std::string src_field;
    (input[i])->ExtractString("SRC_FIELD",src_field);
    int src_matid(-1);
    (input[i])->ExtractInt("SRC_MAT",src_matid);
    std::string tar_field;
    (input[i])->ExtractString("TAR_FIELD",tar_field);
    int tar_matid(-1);
    (input[i])->ExtractInt("TAR_MAT",tar_matid);

    // create the key pair
    std::pair<std::string,std::string> fields(src_field,tar_field);

    // enter the material pairing into the map
    std::pair<int,int> matmap(src_matid,tar_matid);
    clonefieldmatmap_[fields].insert(matmap);
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadTimeFunctionResult(DRT::INPUT::DatFileReader& reader)
{
  //---------------------------------------------- input of time curves
  timecurvemanager_.ReadInput(reader);
  //---------------------------------------- input of spatial functions
  functionmanager_.ReadInput(reader);
  //-------------------------------------- input of result descriptions
  resulttest_.ReadInput(reader);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadConditions(DRT::INPUT::DatFileReader& reader)
{
  Epetra_Time time(*reader.Comm());
  if (reader.Comm()->MyPID()==0)
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
  int ndvol  = design.get<int>("NDVOL");
  int ndparticle  = design.get<int>("NDPARTICLE");

  //--------------------------------------------- read generic node sets
  // read design nodes <-> nodes
  std::vector<std::vector<int> > dnode_fenode(ndnode);
  reader.ReadDesign("DNODE",dnode_fenode);

  // read design lines <-> nodes
  std::vector<std::vector<int> > dline_fenode(ndline);
  reader.ReadDesign("DLINE",dline_fenode);

  // read design surfaces <-> nodes
  std::vector<std::vector<int> > dsurf_fenode(ndsurf);
  reader.ReadDesign("DSURF",dsurf_fenode);

  // read design volumes <-> nodes
  std::vector<std::vector<int> > dvol_fenode(ndvol);
  reader.ReadDesign("DVOL",dvol_fenode);

  // read design particles
  std::vector<std::vector<int> > dparticle(ndparticle);
  reader.ReadDesign("DPARTICLE",dparticle);

  // check for meshfree discretisation to add node set topologies
  std::vector<std::vector<std::vector<int> >* > nodeset(4);
  nodeset[0] = &dnode_fenode;
  nodeset[1] = &dline_fenode;
  nodeset[2] = &dsurf_fenode;
  nodeset[3] = &dvol_fenode;
  std::map<std::string,Teuchos::RCP<Discretization> >::iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
    Teuchos::RCP<DRT::MESHFREE::MeshfreeDiscretization> actdis =
      Teuchos::rcp_dynamic_cast<DRT::MESHFREE::MeshfreeDiscretization>(iter->second);
    if (actdis!=Teuchos::null)
    {
      int foundit = 0;
      for(unsigned n=0; n<4; ++n)
      {
        for(size_t m=0; m<nodeset[n]->size(); ++m)
        {
          std::vector<int> nodes = (*nodeset[n])[m];
          for (unsigned i=0; i<nodes.size(); ++i)
          {
            const int node = nodes[i];
            foundit = actdis->HaveGlobalNode(node);
            if (foundit)
              break;
          }
        }
      }

      int found=0;
      actdis->Comm().SumAll(&foundit,&found,1);
      if (found)
      {
        actdis->AddNodeSetTopology(nodeset);
      }
    }
  }

  // create list of known conditions
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > vc = DRT::INPUT::ValidConditions();
  std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist = *vc;

  // test for each condition definition (input file condition section)
  // - read all conditions that match the definition
  // - add the nodal clouds to the conditions
  // - add the conditions to the appropriate discretizations
  //
  // Note that this will reset (un-FillComplete) the discretizations.
  for (unsigned c=0; c<condlist.size(); ++c)
  {
    std::multimap<int,Teuchos::RCP<DRT::Condition> > cond;

    // read conditions from dat file
    condlist[c]->Read(*this,reader,cond);

    // add nodes to conditions
    std::multimap<int,Teuchos::RCP<DRT::Condition> >::const_iterator curr;
    for (curr=cond.begin(); curr!=cond.end(); ++curr)
    {
      switch (curr->second->GType())
      {
      case Condition::Point:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dnode_fenode.size())
          dserror("DPoint %d not in range [0:%d[\n"
                  "DPoint condition on non existent DPoint?",
                  curr->first,dnode_fenode.size());
        curr->second->Add("Node Ids",dnode_fenode[curr->first]);
        break;
      case Condition::Line:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dline_fenode.size())
          dserror("DLine %d not in range [0:%d[\n"
                  "DLine condition on non existent DLine?",
                  curr->first,dline_fenode.size());
        curr->second->Add("Node Ids",dline_fenode[curr->first]);
        break;
      case Condition::Surface:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dsurf_fenode.size())
          dserror("DSurface %d not in range [0:%d[\n"
                  "DSurface condition on non existent DSurface?",
                  curr->first,dsurf_fenode.size());
        curr->second->Add("Node Ids",dsurf_fenode[curr->first]);
        break;
      case Condition::Volume:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dvol_fenode.size())
          dserror("DVolume %d not in range [0:%d[\n"
                  "DVolume condition on non existent DVolume?",
                  curr->first,dvol_fenode.size());
        curr->second->Add("Node Ids",dvol_fenode [curr->first]);
        break;
      case Condition::Particle:
        if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dparticle.size())
          // particle conditions are allowed for having empty nodal clouds
          break;
        else
           curr->second->Add("Node Ids",dparticle [curr->first]);
         break;
      default:
        dserror("geometry type unspecified");
        break;
      }

      // Iterate through all discretizations and sort the appropriate condition
      // into the correct discretization it applies to

      std::map<std::string,Teuchos::RCP<Discretization> >::iterator iter;
      for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
      {
        Teuchos::RCP<DRT::Discretization> actdis = iter->second;
        if(actdis->Name() != "particle") // standard case
        {
          if(curr->second->GType() == Condition::Particle)
            continue;
          const std::vector<int>* nodes = curr->second->Nodes();
          if (nodes->size()==0)
            dserror("%s condition %d has no nodal cloud",
                     condlist[c]->Description().c_str(),
                     curr->second->Id());

          int foundit = 0;
          for (unsigned i=0; i<nodes->size(); ++i)
          {
            const int node = (*nodes)[i];
            foundit = actdis->HaveGlobalNode(node);
            if (foundit)
              break;
          }
          int found=0;
          actdis->Comm().SumAll(&foundit,&found,1);
          if (found)
          {
            // Insert a copy since we might insert the same condition in many discretizations.
            actdis->SetCondition(condlist[c]->Name(),Teuchos::rcp(new Condition(*curr->second)));
          }
        }
        else
        {
          // insert any particle condition into particle discret
          if(curr->second->GType() == Condition::Particle)
          {
            actdis->SetCondition(condlist[c]->Name(),Teuchos::rcp(new Condition(*curr->second)));
          }
        }
      }
    }
  }

  // debug
#if 0
  for (unsigned i=0; i<NumFields(); ++i)
  {
    for (unsigned j=0; j<NumDis(i); ++j)
    {
      for (unsigned c=0; c<condlist.size(); ++c)
      {
        Teuchos::RCP<DRT::Discretization> actdis = Dis(i,j);
        condlist[c]->Print(cout,&*actdis,true);
      }
    }
  }
#endif

  if (reader.Comm()->MyPID()==0)
  {
    std::cout << time.ElapsedTime() << " secs\n";
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadKnots(DRT::INPUT::DatFileReader& reader)
{
  // get information on the spatial approximation --- we only read knots
  // in the nurbs case
  std::string distype = SpatialApproximation();

  // get problem dimension
  int dim = NDim();

  // Iterate through all discretizations and sort the appropriate condition
  // into the correct discretization it applies to

  std::map<std::string,Teuchos::RCP<Discretization> >::iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
      Teuchos::RCP<DRT::Discretization> actdis = iter->second;

      if(distype == "Nurbs")
      {
        // cast discretisation to nurbs variant to be able
        // to add the knotvector
        DRT::NURBS::NurbsDiscretization* nurbsdis
          =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*actdis));

        if (nurbsdis==NULL)
          dserror("Discretization %s is not a NurbsDiscretization! Panic.", actdis->Name().c_str());

        // define an empty knot vector object
        Teuchos::RCP<DRT::NURBS::Knotvector> disknots=Teuchos::null;

        // read the knotvector data from the input
        reader.ReadKnots(dim,actdis->Name(),disknots);

        if(disknots==Teuchos::null)
        {
          dserror("Knotvector read failed in Nurbs discretisation\n");
        }

        // make sure atdis is fillcompleted, to be able to call
        // ElementRowMap() on it
        // do not initialize elements, since this would require knot
        // vector values
        if(!actdis->Filled())
        {
          actdis->FillComplete(false,false,false);
        }

        // the smallest gid in the discretisation determines the access
        // pattern via the element offset
        int smallest_gid_in_dis=actdis->ElementRowMap()->MinAllGID();

        // consistency checks
        disknots->FinishKnots(smallest_gid_in_dis);

        // add knots to discretisation
        nurbsdis->SetKnotVector(disknots);
      }
  } //loop fields

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenControlFile(const Epetra_Comm& comm, std::string inputfile, std::string prefix, std::string restartkenner)
{
  if (Restart())
    inputcontrol_ = Teuchos::rcp(new IO::InputControl(restartkenner, comm));

  outputcontrol_ = Teuchos::rcp(new IO::OutputControl(comm,
                                                      ProblemName(),
                                                      SpatialApproximation(),
                                                      inputfile,
                                                      restartkenner,
                                                      prefix,
                                                      NDim(),
                                                      Restart(),
                                                      IOParams().get<int>("FILESTEPS"),
                                                      DRT::INPUT::IntegralValue<int>(IOParams(),"OUTPUT_BIN"),
                                                      true));
 if(!DRT::INPUT::IntegralValue<int>(IOParams(),"OUTPUT_BIN"))
  IO::cout<< " Warning no binary Output will be written " << IO::endl;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::OpenErrorFile(const Epetra_Comm& comm, std::string prefix, const bool enforceopening)
{
  bool openfile = enforceopening;
  if (enforceopening == false)
  {
    // what's given in the input file?
    openfile = DRT::INPUT::IntegralValue<int>(IOParams(),"OUTPUT_BIN");
  }
    errorfilecontrol_ = Teuchos::rcp(new IO::ErrorFileControl(comm, prefix, Restart(), openfile));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::WriteInputParameters()
{
  std::string s = OutputControlFile()->FileName();
  s.append(".parameter");
  std::ofstream stream(s.c_str());
  DRT::INPUT::PrintDatHeader(stream,*getParameterList(),"",false,false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadFields(DRT::INPUT::DatFileReader& reader, const bool readmesh)
{

  Teuchos::RCP<DRT::Discretization> structdis       = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> fluiddis        = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> xfluiddis       = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> aledis          = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> structaledis    = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> thermdis        = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> lubricationdis  = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> scatradis       = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> cellscatradis   = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> fluidscatradis  = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> structscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> arterydis       = Teuchos::null; //_1D_ARTERY_
  Teuchos::RCP<DRT::Discretization> airwaydis       = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> optidis         = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> particledis     = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> porofluiddis    = Teuchos::null; // fpsi, poroelast
  Teuchos::RCP<DRT::Discretization> acoudis         = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> celldis         = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> renderingdis    = Teuchos::null; // simple and general rendering discretization for particleMeshFree simulations
  Teuchos::RCP<DRT::Discretization> pboxdis         = Teuchos::null;

  // decide which kind of spatial representation is required
  std::string distype = SpatialApproximation();

  // the basic node reader. now add desired element readers to it!
  DRT::INPUT::NodeReader nodereader(reader, "--NODE COORDS");

  switch (ProblemType())
  {
  case prb_fsi:
  case prb_fsi_redmodels:
  case prb_fsi_lung:
  {
    if(distype == "Nurbs")
    {
      structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
      fluiddis  = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid"    ,reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale"      ,reader.Comm()));
    }
    else if(DRT::INPUT::IntegralValue<int>(FluidDynamicParams().sublist("WALL MODEL"),"X_WALL"))
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure" ,reader.Comm()));
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationXWall("fluid",reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::Discretization("ale",reader.Comm()));
    }
    else
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure" ,reader.Comm()));
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
      if (DRT::INPUT::IntegralValue<bool>(XFluidDynamicParams().sublist("GENERAL"),"XFLUIDFLUID"))
        xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid"   ,reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::Discretization("ale",reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    if (xfluiddis != Teuchos::null)
      xfluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("fluid", fluiddis);
    if (xfluiddis != Teuchos::null)
      AddDis("xfluid", xfluiddis);
    AddDis("ale", aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    if (xfluiddis != Teuchos::null)
    {
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS")));
    }
    else
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_gas_fsi:
  case prb_ac_fsi:
  case prb_thermo_fsi:
  {
    if(distype == "Nurbs")
    {
      dserror("Nurbs discretization not possible for fs3i!");
    }
    else
    {
      structdis       = Teuchos::rcp(new DRT::Discretization("structure" ,reader.Comm()));
      fluiddis        = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
      aledis          = Teuchos::rcp(new DRT::Discretization("ale"       ,reader.Comm()));
      fluidscatradis  = Teuchos::rcp(new DRT::Discretization("scatra1"   ,reader.Comm()));
      structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2"   ,reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(      Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    fluiddis->SetWriter(       Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    aledis->SetWriter(         Teuchos::rcp(new IO::DiscretizationWriter(aledis)));
    fluidscatradis->SetWriter( Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis)));
    structscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structscatradis)));

    AddDis("structure", structdis);
    AddDis("fluid", fluiddis);
    AddDis("ale", aledis);
    AddDis("scatra1", fluidscatradis);
    AddDis("scatra2", structscatradis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    //nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluidscatradis, reader, "--TRANSPORT ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structscatradis, reader, "--TRANSPORT2 ELEMENTS")));

#ifdef EXTENDEDPARALLELOVERLAP
    structdis->CreateExtendedOverlap(false,false,false);
#endif

    break;
  }
  case prb_biofilm_fsi:
  {
      if(distype == "Nurbs")
      {
        dserror("Nurbs discretization not possible for biofilm problems!");
      }
      else
      {
        structdis    = Teuchos::rcp(new DRT::Discretization("structure" ,reader.Comm()));
        fluiddis     = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
        aledis       = Teuchos::rcp(new DRT::Discretization("ale"       ,reader.Comm()));
        structaledis = Teuchos::rcp(new DRT::Discretization("structale" ,reader.Comm()));
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


      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      //nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

  #ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false,false,false);
  #endif

      // fluid scatra field
      fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1",reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis)));
      AddDis("scatra1", fluidscatradis);

      // structure scatra field
      structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2",reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structscatradis)));
      AddDis("scatra2", structscatradis);

      break;
    }
  case prb_fsi_xfem:
  case prb_fluid_xfem:
  {
    structdis = Teuchos::rcp(new DRT::Discretization("structure" ,reader.Comm()));
    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    AddDis("structure", structdis);
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    if(DRT::INPUT::IntegralValue<int>(XFluidDynamicParams().sublist("GENERAL"),"XFLUIDFLUID"))
    {
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      AddDis("fluid", fluiddis);

      xfluiddis  = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid",reader.Comm()));
      xfluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis)));
      AddDis("xfluid", xfluiddis);

      //Not working for XFF... Have to look into it if one wants to use the box geometry feature.
//      nodereader.AddAdvancedReader(fluiddis, reader, "FLUID",
//              DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(),"GEOMETRY"), 0);

      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID")));
    }
    else
    {
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid",reader.Comm()));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      AddDis("fluid", fluiddis);

      nodereader.AddAdvancedReader(fluiddis, reader, "FLUID",
              DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(),"GEOMETRY"), 0);

//      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    }

    aledis    = Teuchos::rcp(new DRT::Discretization("ale"       ,reader.Comm()));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));
    AddDis("ale", aledis);
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));
    break;
  }
  case prb_fpsi_xfem:
  {
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    fluiddis  = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid"    ,reader.Comm()));
    porofluiddis    = Teuchos::rcp(new DRT::Discretization("porofluid"      ,reader.Comm()));
    aledis    = Teuchos::rcp(new DRT::Discretization("ale"      ,reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);
    AddDis("fluid", fluiddis);
    AddDis("ale", aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    nodereader.AddAdvancedReader(fluiddis, reader, "FLUID",
            DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(),"GEOMETRY"), 0);

    //nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_ale:
  {
    if(distype == "Nurbs")
    {
      aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale",reader.Comm()));
    }
    else
    {
      aledis = Teuchos::rcp(new DRT::Discretization("ale",reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("ale", aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid:
  case prb_fluid_redmodels:
  {
    if(distype == "Meshfree")
    {
      fluiddis = Teuchos::rcp(new DRT::MESHFREE::MeshfreeDiscretization("fluid",reader.Comm(),MeshfreeParams()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::MeshfreeDiscretizationWriter(fluiddis)));
    }
    else if(distype == "HDG")
    {
      fluiddis = Teuchos::rcp(new DRT::DiscretizationHDG("fluid",reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    }
    else if(distype == "Nurbs")
    {
      fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));

      // create discretization writer - in constructor set ingto and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    }
    else if(DRT::INPUT::IntegralValue<int>(FluidDynamicParams().sublist("WALL MODEL"),"X_WALL"))
    {
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationXWall("fluid",reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    }
    else
    {
      //fluiddis  = Teuchos::rcp(new DRT::Discretization("fluid",reader.Comm()));
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    }

    AddDis("fluid", fluiddis);

    nodereader.AddAdvancedReader(fluiddis, reader, "FLUID",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(),"GEOMETRY"), 0);

    break;
  }
  case prb_lubrication:
  {
    // create empty discretizations
    lubricationdis = Teuchos::rcp(new DRT::Discretization("lubrication",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    lubricationdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(lubricationdis)));

    AddDis("lubrication", lubricationdis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS")));

    break;
  }
  case prb_scatra_endoexocytosis:
  case prb_cardiac_monodomain:
  case prb_scatra:
  {
    // create empty discretizations
    if(distype == "Meshfree")
    {
      fluiddis = Teuchos::rcp(new DRT::MESHFREE::MeshfreeDiscretization("fluid",reader.Comm(),MeshfreeParams()));
      scatradis = Teuchos::rcp(new DRT::MESHFREE::MeshfreeDiscretization("scatra",reader.Comm(),MeshfreeParams()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::MeshfreeDiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::MeshfreeDiscretizationWriter(scatradis)));
    }
    else
    {
      if(distype == "Nurbs")
      {
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));
        scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra",reader.Comm()));
      }
      else if(distype == "HDG")
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
        scatradis = Teuchos::rcp(new DRT::DiscretizationHDG("scatra",reader.Comm()));
      }
      else
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
        scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
    }

    AddDis("fluid", fluiddis);
    AddDis("scatra", scatradis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

    break;
  }
  case prb_sti:
  {
    // safety checks
    if(distype == "Meshfree" or distype == "Nurbs")
      dserror("Scatra-thermo interaction does not work for meshfree or nurbs discretizations yet!");

    // create empty discretizations for scalar and thermo fields
    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
    thermdis = Teuchos::rcp(new DRT::Discretization("thermo",reader.Comm()));

    // create discretization writers
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
    thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));

    // add empty discretizations to global problem
    AddDis("scatra",scatradis);
    AddDis("thermo",thermdis);

    // add element reader to node reader
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis,reader,"--TRANSPORT ELEMENTS")));

    break;
  }
  case prb_fluid_ale:
  case prb_freesurf:
  {
    if(distype == "Nurbs")
    {
      fluiddis  = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid"    ,reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale"      ,reader.Comm()));
    }
    else if(DRT::INPUT::IntegralValue<int>(FluidDynamicParams().sublist("WALL MODEL"),"X_WALL"))
    {
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationXWall("fluid",reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::Discretization("ale"        ,reader.Comm()));
    }
    else
    {
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationFaces("fluid" ,reader.Comm()));
      if (DRT::INPUT::IntegralValue<bool>(XFluidDynamicParams().sublist("GENERAL"),"XFLUIDFLUID"))
        xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::Discretization("ale"        ,reader.Comm()));
    }


    // create discretization writer - in constructor set into and owned by corresponding discret
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    if (xfluiddis!=Teuchos::null)
    {
      xfluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis)));
    }
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));


    AddDis("fluid", fluiddis);
    if (xfluiddis!=Teuchos::null)
    {
      AddDis("xfluid", xfluiddis); // xfem discretization on slot 1
    }
    AddDis("ale", aledis);

    if (xfluiddis!=Teuchos::null)
    {
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS")));
    }
    else
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_tsi:
  case prb_tfsi_aero:
  {
    if(distype == "Nurbs")
    {
      structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
      thermdis  = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("thermo"   ,reader.Comm()));
    }
    else
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
      thermdis  = Teuchos::rcp(new DRT::Discretization("thermo"   ,reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));

    AddDis("structure", structdis);
    AddDis("thermo", thermdis);

    nodereader.AddAdvancedReader(structdis, reader, "STRUCTURE",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(),"GEOMETRY"), 0);
    nodereader.AddAdvancedReader(thermdis, reader, "THERMO",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(ThermalDynamicParams(),"GEOMETRY"), 0);

    break;
  }
  case prb_thermo:
  {
    thermdis = Teuchos::rcp(new DRT::Discretization("thermo",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    thermdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(thermdis)));

    AddDis("thermo", thermdis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(thermdis, reader, "--THERMO ELEMENTS")));

    break;
  }

  case prb_structure:
  case prb_invana:
  {
    if(distype == "Meshfree")
    {
      dserror("Meshfree structure not implemented, yet.");
      structdis = Teuchos::rcp(new DRT::MESHFREE::MeshfreeDiscretization("structure",reader.Comm(),MeshfreeParams()));
    }
    else if(distype == "Nurbs")
    {
      structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
    }
    else
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));

    AddDis("structure", structdis);

    nodereader.AddAdvancedReader(structdis, reader, "STRUCTURE",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(),"GEOMETRY"), 0);

    break;
  }

  case prb_xcontact:
  {
    // create a discretization for the structure
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));

    AddDis("structure", structdis);

    nodereader.AddAdvancedReader(structdis, reader, "STRUCTURE",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(),"GEOMETRY"), 0);

    // This is moved to XCONTACT::ALGORITHM::Base::Setup()
//    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
//    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
//    AddDis("scatra", scatradis);
    // seems unnecessary, since we do not read the scatra mesh from the dat file
//    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
//    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
    break;
  }

  case prb_polymernetwork:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization( "structure", reader.Comm() ) );
    pboxdis  = Teuchos::rcp(new DRT::Discretization( "boundingbox" , reader.Comm() ) );

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter( Teuchos::rcp( new IO::DiscretizationWriter( structdis ) ) );
    pboxdis->SetWriter( Teuchos::rcp( new IO::DiscretizationWriter( pboxdis ) ) );

    AddDis( "structure", structdis );
    AddDis( "boundingbox", pboxdis );

    nodereader.AddElementReader( Teuchos::rcp( new DRT::INPUT::ElementReader( structdis, reader, "--STRUCTURE ELEMENTS") ) );
    nodereader.AddElementReader( Teuchos::rcp( new DRT::INPUT::ElementReader( pboxdis, reader, "--PERIODIC BOUNDINGBOX ELEMENTS" ) ) );

    break;
  }

  case prb_loma:
  {
    // create empty discretizations
    fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

    AddDis("fluid", fluiddis);
    AddDis("scatra", scatradis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

    break;
  }

  case prb_two_phase_flow:
  case prb_fluid_xfem_ls:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure" ,reader.Comm()));
      if(ProblemType()==prb_fluid_xfem_ls)
        fluiddis  = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid",reader.Comm()));
      else
        fluiddis  = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
      particledis = Teuchos::rcp(new DRT::Discretization("particle",reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
      fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
      scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
      particledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particledis)));

      AddDis("structure", structdis);
      AddDis("fluid", fluiddis);
      AddDis("scatra", scatradis);
      AddDis("particle", particledis);

      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
      nodereader.AddAdvancedReader(fluiddis, reader, "FLUID",
              DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(),"GEOMETRY"), 0);
      //nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
      nodereader.AddParticleReader(particledis, reader, "PARTICLE");
      break;
    }

  case prb_elch:
  {
    // create empty discretizations
    if(distype == "Nurbs")
    {
      fluiddis  = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid",reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra",reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale",reader.Comm()));
    }
    else
    {
      fluiddis  = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra",    reader.Comm()));
      aledis    = Teuchos::rcp(new DRT::Discretization("ale",       reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("fluid", fluiddis);
    AddDis("scatra", scatradis);
    AddDis("ale", aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis,reader, "--TRANSPORT ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis,   reader, "--ALE ELEMENTS")));

    break;
  }

  case prb_combust:
  {
    // create empty discretizations
    fluiddis = Teuchos::rcp(new DRT::DiscretizationCombust("fluid",reader.Comm()));
    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
    particledis = Teuchos::rcp(new DRT::Discretization("particle",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
    particledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particledis)));

    AddDis("fluid", fluiddis);
    AddDis("scatra", scatradis);
    AddDis("particle", particledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddParticleReader(particledis, reader, "PARTICLE");
    break;
  }

  case prb_art_net: // _1D_ARTERY_
  {
    // create empty discretizations
    arterydis = Teuchos::rcp(new DRT::Discretization("artery",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));

    AddDis("artery", arterydis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));

    break;
  }
  case prb_red_airways: // _reduced D airways
  {
    // create empty discretizations
    airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));

    AddDis("red_airway", airwaydis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));

    break;
  }
  case prb_struct_ale: // structure with ale
  {
    // create empty discretizations
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    aledis  = Teuchos::rcp(new DRT::Discretization("ale"   ,reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("ale", aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_fluid_topopt:
  {
    // create empty discretizations
    fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
    optidis = Teuchos::rcp(new DRT::Discretization("opti",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    optidis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(optidis)));

    AddDis("fluid", fluiddis);
    AddDis("opti", optidis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));

    break;
  }
  case prb_poroelast:
  case prb_poromultiphase:
  {
    // create empty discretizations
    if(distype == "Nurbs")
    {
      structdis  = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid",reader.Comm()));
    }
    else
    {
      structdis     = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
      porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid"   ,reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));

    break;
  }
  case prb_poromultiphasescatra:
  {
    // create empty discretizations
    if(distype == "Nurbs")
    {
      structdis  = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid",reader.Comm()));
      scatradis     = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra",reader.Comm()));
    }
    else
    {
      structdis     = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
      porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid",reader.Comm()));
      scatradis     = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);
    AddDis("scatra",    scatradis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
    break;
  }
  case prb_porofluidmultiphase:
  {
    // create empty discretizations
    if(distype == "Nurbs")
    {
      porofluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid",reader.Comm()));
    }
    else
    {
      porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid"   ,reader.Comm()));
    }

    // create discretization writer - in constructor set into and owned by corresponding discret
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));

    AddDis("porofluid", porofluiddis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));

    break;
  }
  case prb_fpsi:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
    porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
    fluiddis      = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
    aledis        = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);
    AddDis("fluid",     fluiddis);
    AddDis("ale",   aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis,  reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    break;
  }
  case prb_immersed_fsi:
  case prb_immersed_membrane_fsi:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
    fluiddis      = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    fluiddis ->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));

    AddDis("structure", structdis);
    AddDis("fluid",     fluiddis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis,  reader, "--FLUID ELEMENTS")));

    break;
  }
  case prb_immersed_ale_fsi:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
    fluiddis      = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
    aledis        = Teuchos::rcp(new DRT::DiscretizationFaces("ale",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    fluiddis ->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    aledis ->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("fluid",     fluiddis);
    AddDis("ale",     aledis);

    std::set<std::string> fluidelementtypes;
    fluidelementtypes.insert("FLUIDIMMERSED");

    std::set<std::string> structelementtypes;
    fluidelementtypes.insert("SOLIDH8");

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS",structelementtypes)));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis,  reader, "--FLUID ELEMENTS",fluidelementtypes)));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(aledis,  reader, "--ALE ELEMENTS")));

    break;
  }
  case prb_immersed_cell:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid"   ,reader.Comm()));
    celldis       = Teuchos::rcp(new DRT::Discretization("cell",reader.Comm()));
    scatradis     = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
    cellscatradis = Teuchos::rcp(new DRT::Discretization("cellscatra",reader.Comm()));
    aledis        = Teuchos::rcp(new DRT::Discretization("ale",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
    celldis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(celldis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
    cellscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(cellscatradis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);
    AddDis("cell", celldis);
    AddDis("scatra", scatradis);
    AddDis("cellscatra", cellscatradis);
    AddDis("ale", aledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(celldis, reader, "--CELL ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(cellscatradis, reader, "--CELLSCATRA ELEMENTS")));

    break;
  }
  case prb_fps3i:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization("structure",  reader.Comm()));
    porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid",  reader.Comm()));
    fluiddis      = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
    aledis        = Teuchos::rcp(new DRT::Discretization("ale",        reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);
    AddDis("fluid",     fluiddis);
    AddDis("ale",     aledis);


    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(fluiddis,  reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));

    // fluid scatra field
    fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1",reader.Comm()));
    // create discretization writer - in constructor set into and owned by corresponding discret
    fluidscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis)));
    AddDis("scatra1", fluidscatradis);

    // poro structure scatra field
    structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2",reader.Comm()));
    // create discretization writer - in constructor set into and owned by corresponding discret
    structscatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structscatradis)));
    AddDis("scatra2", structscatradis);

    break;
  }
  case prb_poroscatra:
  {
    // create empty discretizations
    structdis     = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
    porofluiddis  = Teuchos::rcp(new DRT::Discretization("porofluid" ,reader.Comm()));
    scatradis     = Teuchos::rcp(new DRT::Discretization("scatra",    reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    porofluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

    AddDis("structure", structdis);
    AddDis("porofluid", porofluiddis);
    AddDis("scatra",    scatradis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
    break;
  }
  case prb_ehl:
  {
    // create empty discretizations
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    lubricationdis = Teuchos::rcp(new DRT::Discretization("lubrication",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    lubricationdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(lubricationdis)));

    AddDis("structure", structdis);
    AddDis("lubrication", lubricationdis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS")));

    break;
  }
  case prb_ssi:
  {
    // create empty discretizations
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

    AddDis("structure", structdis);
    AddDis("scatra", scatradis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

    break;
  }
  case prb_particle:
  {
    if(distype == "Meshfree")
    {
      particledis = Teuchos::rcp(new DRT::Discretization("particle",reader.Comm()));
    }
    else
    {
      dserror("particle simulations must be distype=Meshfree");
    }
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    particledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particledis)));

    AddDis("structure", structdis);
    AddDis("particle", particledis);

    nodereader.AddAdvancedReader(structdis, reader, "STRUCTURE",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(),"GEOMETRY"), 0);
    nodereader.AddParticleReader(particledis, reader, "PARTICLE");

    break;
  }
  case prb_cavitation:
  {
    // create empty discretizations
    fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid",reader.Comm()));
    particledis = Teuchos::rcp(new DRT::Discretization("particle",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    fluiddis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(fluiddis)));
    particledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particledis)));

    AddDis("fluid", fluiddis);
    AddDis("particle", particledis);

    nodereader.AddAdvancedReader(fluiddis, reader, "FLUID",
        DRT::INPUT::IntegralValue<INPAR::GeometryType>(FluidDynamicParams(),"GEOMETRY"), 0);
    nodereader.AddParticleReader(particledis, reader, "PARTICLE");

    break;
  }
  case prb_meshfree:
  case prb_pasi:
  {
    if(distype == "Meshfree")
    {
     particledis = Teuchos::rcp(new DRT::Discretization("particle",reader.Comm()));
    }
    else
    {
     dserror("meshfree simulations must be distype=Meshfree");
    }
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));


    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    particledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particledis)));


    AddDis("structure", structdis);
    AddDis("particle", particledis);

    nodereader.AddAdvancedReader(structdis, reader, "STRUCTURE",
       DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(),"GEOMETRY"), 0);
    nodereader.AddParticleReader(particledis, reader, "PARTICLE");

    // add rendering in case is required
    if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"RENDERING"))
    {
      renderingdis = Teuchos::rcp(new DRT::Discretization("rendering",reader.Comm()));
      renderingdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(renderingdis)));
      AddDis("rendering", renderingdis);
      nodereader.AddAdvancedReader(renderingdis, reader, "MESHFREE RENDERING",
          INPAR::geometry_box, 0);
    }
    break;
  }
  case prb_level_set:
  {
    // create empty discretizations
    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));
    particledis = Teuchos::rcp(new DRT::Discretization("particle",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));
    particledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particledis)));

    AddDis("scatra", scatradis);
    AddDis("particle", particledis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));
    nodereader.AddParticleReader(particledis, reader, "PARTICLE");
    break;
  }
  case prb_np_support:
  {
    // no discretizations and nodes needed for supporting procs
    break;
  }
  case prb_acou:
  {
    // create empty discretizations
    acoudis = Teuchos::rcp(new DRT::DiscretizationHDG("acou",reader.Comm()));
    scatradis = Teuchos::rcp(new DRT::Discretization("scatra",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    acoudis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(acoudis)));
    scatradis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(scatradis)));

    AddDis("acou", acoudis);
    AddDis("scatra", scatradis);

    std::set<std::string> acouelementtypes;
    acouelementtypes.insert("ACOUSTIC");
    acouelementtypes.insert("ACOUSTICSOL");

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(acoudis, reader, "--ACOUSTIC ELEMENTS",acouelementtypes)));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS")));

    break;
  }
  case prb_redairways_tissue:
  {
    // create empty discretizations
    structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
    airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway",reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));
    airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));

    AddDis("structure", structdis);
    AddDis("red_airway", airwaydis);

    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));
  }
  break;
  case prb_uq:
  {
    // tiny bit brutal
    //get list for multi level monte carlo
    const Teuchos::ParameterList& mlmcp = this->MultiLevelMonteCarloParams();
    INPAR::MLMC::FWDProblem fwdprb = DRT::INPUT::IntegralValue<INPAR::MLMC::FWDProblem>(mlmcp,"FWDPROBLEM");

    // quick check whether to read in structure or airways
    if (fwdprb==INPAR::MLMC::structure)
    {
      if(distype == "Meshfree")
      {
        dserror("Meshfree structure not implemented, yet.");
        structdis = Teuchos::rcp(new DRT::MESHFREE::MeshfreeDiscretization("structure",reader.Comm(),MeshfreeParams()));
      }
      else if(distype == "Nurbs")
      {
        dserror("Meshfree structure not implemented, yet.");
        structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure",reader.Comm()));
      }
      else
      {
        structdis = Teuchos::rcp(new DRT::Discretization("structure",reader.Comm()));
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis)));

      AddDis("structure", structdis);

      nodereader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          DRT::INPUT::IntegralValue<INPAR::GeometryType>(StructuralDynamicParams(),"GEOMETRY"), 0);

    }
    else if (fwdprb==INPAR::MLMC::red_airways)
    {
      // create empty discretizations
      airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway",reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));

      AddDis("red_airway", airwaydis);

      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));
    }
    else
    {
      dserror("Uncertainty quantification is only implemented for structure or airways ");
    }
  }
  break;
  case prb_tutorial:
  {

  }
  break;
  default:
    dserror("Unknown problem type: %d",ProblemType());
    break;
  }

  // add artery or airways discretizations only for the following problem types
  switch (ProblemType())
  {
  case prb_fsi_redmodels:
  case prb_fsi_lung:
  case prb_fluid_ale:
  case prb_fluid_redmodels:
  {
    if(distype == "Polynomial")
    {
      // create empty discretizations
      arterydis = Teuchos::rcp(new DRT::Discretization("artery",reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      arterydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(arterydis)));
      AddDis("artery", arterydis);
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS")));

      airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway",reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      airwaydis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(airwaydis)));
      AddDis("red_airway", airwaydis);
      nodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS")));
    }
  }
  break;
  default:
    break;
  }

  if (readmesh) // now read and allocate!
  {
    // we read nodes and elements for the desired fields as specified above
    nodereader.Read();

    NP_TYPE npType = DRT::Problem::Instance()->GetNPGroup()->NpType();
    // care for special applications
    switch (ProblemType())
    {
    case prb_fsi:
    case prb_fsi_redmodels:
    case prb_fsi_lung:
    case prb_scatra:
    {
      // read microscale fields from second, third, ... inputfile if necessary
      // (in case of multi-scale material models in structure field)
      if (npType != copy_dat_file) ReadMicroFields(reader);
      break;
    }
    case prb_structure:
    case prb_invana:
    {
      // read microscale fields from second, third, ... inputfile if necessary
      // (in case of multi-scale material models)
      if (npType != copy_dat_file) ReadMicroFields(reader);
      break;
    }
    case prb_np_support:
    {
      // read microscale fields from second, third, ... inputfile for supporting processors
      ReadMicrofields_NPsupport();
      break;
    }
    case prb_tutorial:
    {
      break;
    }
    default:
      break;
    }
  } // if(readmesh)

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMicroFields(DRT::INPUT::DatFileReader& reader)
{
  // check whether micro material is specified
  const int id_struct = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_struct_multiscale);
  const int id_scatra = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_scatra_multiscale);

  // return if no multiscale material is used
  if(id_struct == -1 and id_scatra == -1)
    return;

  // safety check
  if(id_struct != -1 and id_scatra != -1)
    dserror("Cannot have multi-scale material for both structure and scalar transport!");

  // store name of macro-scale discretization in string
  std::string macro_dis_name("");
  if(id_struct != -1)
    macro_dis_name = "structure";
  else
    macro_dis_name = "scatra";

  // fetch communicators
  Teuchos::RCP<Epetra_Comm> lcomm = npgroup_->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = npgroup_->GlobalComm();

  DRT::Problem* macro_problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> macro_dis = macro_problem->GetDis(macro_dis_name);

  // repartition macro problem for a good distribution of elements with micro material
  if(macro_dis_name == "structure")
    DRT::UTILS::WeightedRepartitioning(macro_dis, true, true, true);

  // make sure that we read the micro discretizations only on the processors on
  // which elements with the corresponding micro material are evaluated

  std::set<int> my_multimat_IDs;

  // take care also of ghosted elements! -> ElementColMap!
  for (int i=0; i<macro_dis->ElementColMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = macro_dis->lColElement(i);
    Teuchos::RCP<MAT::Material> actmat = actele->Material();

    if((actmat->MaterialType() == INPAR::MAT::m_struct_multiscale and macro_dis_name == "structure")
    or (actmat->MaterialType() == INPAR::MAT::m_scatra_multiscale and macro_dis_name == "scatra"))
    {
      MAT::PAR::Parameter* actparams = actmat->Parameter();
      my_multimat_IDs.insert(actparams->Id());
    }
  }

  // check which macro procs have an element with micro material
  int foundmicromat = 0;
  int foundmicromatmyrank = -1;
  if(my_multimat_IDs.size() != 0)
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

  // determine color of macro procs with any contribution to micro material, only important for procs with micro material
  // color starts with 0 and is incremented for each group
  int color = -1;
  if(foundmicromat == 1)
  {
    for(int t=0; t<(int)foundmyranks.size(); t++)
    {
      if(foundmyranks[t] != -1)
        ++color;
      if(foundmyranks[t] == foundmicromatmyrank)
        break;
    }
  }
  else
  {
    color = MPI_UNDEFINED;
  }

  // do the splitting of the communicator (macro proc must always be proc in subcomm with lowest key --> 0 is inserted here)
  MPI_Comm  mpi_local_comm;
  MPI_Comm_split((Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm,true)->GetMpiComm()),color,0/*important here*/,&mpi_local_comm);

  // sort out macro procs that do not have micro material
  if(foundmicromat == 1)
  {
    // create the sub communicator that includes one macro proc and some supporting procs
    Teuchos::RCP<Epetra_Comm> subgroupcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
    npgroup_->SetSubComm(subgroupcomm);

    // find out how many micro problems have to be solved on this macro proc
    int microcount = 0;
    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials_->Map()->begin();
         i!=materials_->Map()->end();
         ++i)
    {
      int matid = i->first;
      if (my_multimat_IDs.find(matid)!=my_multimat_IDs.end())
        microcount++;
    }
    // and broadcast it to the corresponding group of procs
    subgroupcomm->Broadcast(&microcount, 1, 0);

    for (std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=materials_->Map()->begin();
         i!=materials_->Map()->end();
         ++i)
    {
      int matid = i->first;

      if (my_multimat_IDs.find(matid)!=my_multimat_IDs.end())
      {
        Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);

        // initialize variables storing micro-scale information
        int microdisnum(-1);
        std::string micro_dis_name = "";
        std::string micro_inputfile_name("");
        DRT::Problem* micro_problem(NULL);

        // structure case
        if(macro_dis_name == "structure")
        {
          // access multi-scale structure material
          MAT::MicroMaterial* micromat = static_cast<MAT::MicroMaterial*>(mat.get());

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->MicroDisNum();
          subgroupcomm->Broadcast(&microdisnum,1,0);

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
          MAT::ScatraMatMultiScale* micromat = static_cast<MAT::ScatraMatMultiScale*>(mat.get());

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->MicroDisNum();
          subgroupcomm->Broadcast(&microdisnum,1,0);

          // set unique name of micro-scale discretization
          std::stringstream name;
          name << "scatra_multiscale_" << microdisnum;
          micro_dis_name = name.str();

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->MicroInputFileName();

          // instantiate micro-scale problem
          micro_problem = DRT::Problem::Instance(microdisnum);
        }

        if (micro_inputfile_name[0]!='/')
        {
          std::string filename = reader.MyInputfileName();
          std::string::size_type pos = filename.rfind('/');
          if (pos!=std::string::npos)
          {
            std::string path = filename.substr(0,pos+1);
            micro_inputfile_name.insert(micro_inputfile_name.begin(), path.begin(), path.end());
          }
        }

        // broadcast micro input file name
        int length = micro_inputfile_name.length();
        subgroupcomm->Broadcast(&length, 1, 0);
        subgroupcomm->Broadcast((const_cast<char *>(micro_inputfile_name.c_str())), length, 0);

        // start with actual reading
        DRT::INPUT::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

        Teuchos::RCP<DRT::Discretization> dis_micro = Teuchos::rcp(new DRT::Discretization(micro_dis_name,micro_reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        dis_micro->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(dis_micro)));

        micro_problem->AddDis(micro_dis_name,dis_micro);

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

        DRT::INPUT::NodeReader micronodereader(micro_reader, "--NODE COORDS");

        if(micro_dis_name == "structure")
          micronodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(dis_micro,micro_reader,"--STRUCTURE ELEMENTS")));
        else
          micronodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(dis_micro,micro_reader,"--TRANSPORT ELEMENTS")));

        micronodereader.Read();

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

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::ReadMicrofields_NPsupport()
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetNPGroup()->GlobalComm();

  // receive number of procs that have micro material
  int nummicromat = 0;
  gcomm->Broadcast(&nummicromat, 1, 0);

  // prepare the supporting procs for a splitting of gcomm

  // groups should be equally sized
  // in a first step every macro proc that needs support gets procpergroup supporting procs
  int procpergroup = int(floor((lcomm->NumProc())/nummicromat));
  std::vector<int> supgrouplayout(nummicromat, procpergroup);
  // remaining procs are added to the groups in the beginning
  int remainingProcs = lcomm->NumProc() - procpergroup*nummicromat;
  for(int k=0; k<remainingProcs; ++k)
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
  }
  while( gsum <= lcomm->MyPID() );

  // do the splitting of the communicator
  MPI_Comm  mpi_local_comm;
  MPI_Comm_split((Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm,true)->GetMpiComm()),color,gcomm->MyPID(),&mpi_local_comm);

  // create the sub communicator that includes one macro proc and some supporting procs
  Teuchos::RCP<Epetra_Comm> subgroupcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
  npgroup_->SetSubComm(subgroupcomm);

  // number of micro problems for this sub group
  int microcount = 0;
  subgroupcomm->Broadcast(&microcount, 1, 0);

  for(int n=0; n<microcount; n++)
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
    subgroupcomm->Broadcast((const_cast<char *>(micro_inputfile_name.c_str())), length, 0);

    // start with actual reading
    DRT::INPUT::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

    Teuchos::RCP<DRT::Discretization> structdis_micro = Teuchos::rcp(new DRT::Discretization("structure", micro_reader.Comm()));

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

    DRT::INPUT::NodeReader micronodereader(micro_reader, "--NODE COORDS");
    micronodereader.AddElementReader(Teuchos::rcp(new DRT::INPUT::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS")));
    micronodereader.Read();

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

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::setParameterList( Teuchos::RCP< Teuchos::ParameterList > const &paramList )
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
void DRT::Problem::AddDis(const std::string name, Teuchos::RCP<Discretization> dis)
{
  // safety checks
  if (dis == Teuchos::null) dserror ("Received Teuchos::null.");
  if (dis->Name().empty())  dserror ("Discretization has empty name string.");

  if (discretizationmap_.insert(std::make_pair(name, dis)).second == false)
  {
    // if the same key already exists we have to inform the user since
    // the insert statement did not work in this case
    dserror("Could not insert discretization '%s' under (duplicate) key '%s'.",dis->Name().c_str(),name.c_str());
  }
  // For debug: what's currently in the map:
  /*
  std::map<std::string,Teuchos::RCP<Discretization> >::iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
    cout<<"key : "<<iter->first<<"    "<<"discret.name = "<<iter->second->Name()<<endl<<endl;
  }
  */
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> DRT::Problem::GetDis(const std::string name) const
{
  std::map<std::string,Teuchos::RCP<Discretization> >::const_iterator iter = discretizationmap_.find(name);

  if (iter != discretizationmap_.end())
  {
    return iter->second;
  }
  else
  {
    dserror("Could not find discretization '%s'.",name.c_str());
    return Teuchos::null;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<std::string> DRT::Problem::GetDisNames() const {

  unsigned mysize = NumFields();
  std::vector<std::string> vec;
  vec.reserve(mysize);

  std::map<std::string,Teuchos::RCP<Discretization> >::const_iterator iter;
  for (iter = discretizationmap_.begin(); iter != discretizationmap_.end(); ++iter)
  {
    vec.push_back(iter->first);
  }

  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::Problem::DoesExistDis(const std::string name) const {

  std::map<std::string,Teuchos::RCP<Discretization> >::const_iterator iter = discretizationmap_.find(name);
  if (iter != discretizationmap_.end())
  {
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::Function& DRT::Problem::Funct(int num)
{
  return functionmanager_.Funct(num);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::TimeCurve& DRT::Problem::Curve(int num)
{
  return timecurvemanager_.Curve(num);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Problem::SetRestartStep(int r)
{
  if (r<0) dserror("Negative restart step not allowed");

  restartstep_ = r;
}
