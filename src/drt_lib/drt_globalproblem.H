/*----------------------------------------------------------------------*/
/*!
\file drt_globalproblem.H

\brief internal representation of an input file

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*/
/*----------------------------------------------------------------------*/


#ifndef DRT_GLOBALPROBLEM_H
#define DRT_GLOBALPROBLEM_H

#include <vector>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include "drt_globalproblem_enums.H"
#include "drt_function.H"
#include "drt_timecurve.H"
#include "drt_resulttest.H"
#include "drt_utils.H"

namespace IO
{
  class OutputControl;
  class InputControl;
  class ErrorFileControl;
}

namespace COMM_UTILS
{
  class NestedParGroup;
}

namespace MAT
{
namespace PAR
{
class Bundle;
}
}


namespace DRT
{
class Discretization;
class SingletonDestruction;

namespace INPUT
{
class DatFileReader;
}


/// Global problem instance that keeps the discretizations
/*!

  The global problem represents the input file passed to baci. This
  class organizes the reading of a dat file (utilizing the
  DatFileReader of course). That is way, in all but the most eccentric
  cases there will be exactly one object of this class during a baci
  run. This object contains all parameters read from the input file as
  well as any material definitions and even all the discretizations.

  <h3>Input parameters</h3>

  All input parameters are known by the global problem object. These
  parameters are guaranteed to be valid (because they passed the
  validation) and are guaranteed to be there (because default values
  have been set for all parameters missing from the input file.) This
  is Teuchos::ParameterList magic, that just requires the list of
  valid parameters in the file drt_validparameters.cpp to be
  complete.

  The algorithms are meant to ask the global problem object for their
  parameters and extract them from the respective parameter list
  directly.

  <h3>Discretizations</h3>

  The global problem object knows the discretizations defined by the
  input file. In order to access a particular discretization you get
  the global problem object and ask.

  <h3>Materials</h3>

  The global problem object knows the material descriptions from the
  input file. These are not to be confused with the material classes
  the elements know and work with. The global problem object does not
  keep track of gauss point material values, all that is known here
  are the definitions from the input file.

  <h3>Purpose</h3>

  baci evolved out of ccarat in a (somewhat) smooth fashion. The old
  ccarat used global C struct variables all over the place. It is
  indeed the aim of the global problem to collect these variables in
  one place and provide a sane interface to the data. The transition
  is already pretty advanced, however there are still many places
  where global variables are used. These places will have to be
  revised some day.

  The global problem class is still evolving.

  \note It is illegal to introduce new global variables in baci. It is
  illegal to use global variables where other means to do the job are
  available. We just did not get around to remove all of them yet.

  <h3>Singleton behavior</h3>

  \warning This is a guru only section!

  The global problem behaves like a singleton, so there is always one
  instance available. But you can have more than one instance of
  Problem. In normal situations this will not be needed. So don't
  bother. Just call the static Instance() function to get the global
  instance and access your discretizations.

  In the special case that you want to read more that one input file,
  however, you will need to handle the fields from each file
  separately. This is something you cannot do with ccarat's field
  variable so in that case you will have to rely on a second Problem
  instance.

  One artefact that comes from using global variables together with
  multiple Problem objects is the notion of activating of problem
  object. This translates to setting the global variable pointers to
  this Problem object's internal variables. Normally, if there is just
  one Problem object, this is done by default. If you need more that
  one, however, you will have to activate the global problems
  yourself.

  \author u.kue
  \date 06/07
 */
class Problem : public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /// @name Instances

  /// return an instance of this class
  static Problem* Instance(int num=0);

  /// return number of problem instances
  static unsigned NumInstances() { return instances_.size(); }

  /// calculation done, clean up
  /*!
    There can be a variety of objects to a problem. Some of them might
    require proper cleanup. Make sure we always do it.
   */
  static void Done();

  /// register singleton objects that need explicit destruction
  void Register(DRT::SingletonDestruction * sd) { sds_.push_back(sd); }

  //@}

  ~Problem();

  /// @name Input

  /// read global parameters
  void ReadParameter(DRT::INPUT::DatFileReader& reader);

  /// setup the discretizations
  void ReadFields(DRT::INPUT::DatFileReader& reader, const bool readmesh = true);

  /// set up micro-scale discretizations in case of structure multi-scale simulations
  void ReadMicroFields(DRT::INPUT::DatFileReader& reader);

  /// set up supporting processors for micro-scale discretizations
  void ReadMicrofields_NPsupport();

  /// input of materials
  void ReadMaterials(DRT::INPUT::DatFileReader& reader);

  /// setup map between materials of original and cloned elements
  void ReadCloningMaterialMap(DRT::INPUT::DatFileReader& reader);

  /// input of conditions
  void ReadConditions(DRT::INPUT::DatFileReader& reader);

  /// input of time curves, functions and result tests
  void ReadTimeFunctionResult(DRT::INPUT::DatFileReader& reader);

  /// input of knots for isogeometric analysis
  void ReadKnots(DRT::INPUT::DatFileReader& reader);

  /// set restart step which was read from the command line
  void SetRestartStep(int r);

  void SetInputControlFile(Teuchos::RCP<IO::InputControl>& input){inputcontrol_ = input;}

  //@}

  /// @name General query methods
  /// Once and for all definitions

  /// give enum of my problem type
  PROBLEM_TYP ProblemType() const;

  /// give string name of my problem type
  std::string ProblemName() const;

  /// return restart step
  int Restart() const;

  /// return restart time
  double RestartTime() const;

  /// number of space dimensions (as specified in the input file)
  int NDim() const;

  /// current wall time in seconds with micro second precision
  static double Walltime();

  // do bandwith optimization of dofsets (as specified in the input file)
  //bool BandWidthOpt() const;

  // give string name of my basis function type (polynomial/nurbs)
  std::string SpatialApproximation() const;

  //@}

  /// @name Control file

  /// create control file for output and read restart data if required
  void OpenControlFile(const Epetra_Comm& comm, std::string inputfile, std::string prefix, std::string restartkenner);

  /// error files
  void OpenErrorFile(const Epetra_Comm& comm, std::string prefix, const bool enforceopening = false);

  /// control file for restart read
  Teuchos::RCP<IO::InputControl> InputControlFile() { return inputcontrol_; }

  /// control file for normal output
  Teuchos::RCP<IO::OutputControl> OutputControlFile() { return outputcontrol_; }

  /// control error files
  Teuchos::RCP<IO::ErrorFileControl> ErrorFile() { return errorfilecontrol_; }

  /// write parameters read from input file for documentation
  void WriteInputParameters();

  //@}

  /// @name Parameters read from file

  /// Set parameters from a parameter list and return with default values.
  void setParameterList(Teuchos::RCP< Teuchos::ParameterList > const &paramList);

  /// Return a const parameter list of all of the valid parameters that this->setParameterList(...) will accept.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// @name Nested parallelism group

  /// nested parallelism group is initialized
  void NPGroup(
    int groupId,  //!< group id
    int ngroup, //!< number of groups in nested parallelism
    std::map<int, int> lpidgpid,  //!< map of local and global proc ids
    Teuchos::RCP<Epetra_Comm> lcomm,  //!< local communicator
    Teuchos::RCP<Epetra_Comm> gcomm,  //!< global communicator
    NP_TYPE npType  //!< nested parallelism type
    );

  /// nested parallelism group initialized from an existing NPgroup
  void NPGroup(COMM_UTILS::NestedParGroup& npgroup //!< copy from npgroup
      );

  /// return nested parallelism group
  Teuchos::RCP<COMM_UTILS::NestedParGroup> GetNPGroup();

  //@}

  //@}

  /// @name Input parameter sections
  /// direct access to parameters from input file sections

  const Teuchos::ParameterList& MeshfreeParams() const { return getParameterList()->sublist("MESHFREE"); }
  const Teuchos::ParameterList& IOParams() const { return getParameterList()->sublist("IO"); }
  const Teuchos::ParameterList& DesignDescriptionParams() const { return getParameterList()->sublist("DESIGN DESCRIPTION"); }
  const Teuchos::ParameterList& PatSpecParams() const { return getParameterList()->sublist("PATIENT SPECIFIC"); }
  const Teuchos::ParameterList& StructuralDynamicParams() const { return getParameterList()->sublist("STRUCTURAL DYNAMIC"); }
  const Teuchos::ParameterList& Cardiovascular0DStructuralParams() const { return getParameterList()->sublist("CARDIOVASCULAR 0D-STRUCTURE COUPLING"); }
  const Teuchos::ParameterList& InverseAnalysisParams() const { return getParameterList()->sublist("INVERSE ANALYSIS"); }
  const Teuchos::ParameterList& StatInverseAnalysisParams() const { return getParameterList()->sublist("STAT INVERSE ANALYSIS"); }
  const Teuchos::ParameterList& MultiLevelMonteCarloParams() const { return getParameterList()->sublist("MULTI LEVEL MONTE CARLO"); }
  const Teuchos::ParameterList& MortarCouplingParams() const { return getParameterList()->sublist("MORTAR COUPLING"); }
  const Teuchos::ParameterList& ContactDynamicParams() const { return getParameterList()->sublist("CONTACT DYNAMIC"); }
  const Teuchos::ParameterList& FluctuatingHydrodynamicsParams() const { return getParameterList()->sublist("FLUCTUATING HYDRODYNAMICS"); }
  const Teuchos::ParameterList& BeamInteractionParams() const { return getParameterList()->sublist("BEAM INTERACTION"); }
  const Teuchos::ParameterList& BrownianDynamicsParams() const { return getParameterList()->sublist("BROWNIAN DYNAMICS"); }
  const Teuchos::ParameterList& CrosslinkingParams() const { return getParameterList()->sublist("CROSSLINKING"); }
  const Teuchos::ParameterList& ThermalDynamicParams() const { return getParameterList()->sublist("THERMAL DYNAMIC"); }
  const Teuchos::ParameterList& TSIDynamicParams() const { return getParameterList()->sublist("TSI DYNAMIC"); }
  const Teuchos::ParameterList& FluidDynamicParams() const { return getParameterList()->sublist("FLUID DYNAMIC"); }
  const Teuchos::ParameterList& TwoPhaseFlowParams() const { return getParameterList()->sublist("TWO PHASE FLOW"); }
  const Teuchos::ParameterList& LubricationDynamicParams() const { return getParameterList()->sublist("LUBRICATION DYNAMIC"); }
  const Teuchos::ParameterList& ScalarTransportDynamicParams() const { return getParameterList()->sublist("SCALAR TRANSPORT DYNAMIC"); }
  const Teuchos::ParameterList& STIDynamicParams() const { return getParameterList()->sublist("STI DYNAMIC"); }
  const Teuchos::ParameterList& FS3IDynamicParams() const { return getParameterList()->sublist("FS3I DYNAMIC"); }
  const Teuchos::ParameterList& CombustionDynamicParams() const { return getParameterList()->sublist("COMBUSTION CONTROL"); }
  const Teuchos::ParameterList& AleDynamicParams() const { return getParameterList()->sublist("ALE DYNAMIC"); }
  const Teuchos::ParameterList& FSIDynamicParams() const { return getParameterList()->sublist("FSI DYNAMIC"); }
  const Teuchos::ParameterList& FPSIDynamicParams() const { return getParameterList()->sublist("FPSI DYNAMIC"); }
  const Teuchos::ParameterList& ImmersedMethodParams() const { return getParameterList()->sublist("IMMERSED METHOD"); }
  const Teuchos::ParameterList& CellMigrationParams() const { return getParameterList()->sublist("CELL DYNAMIC"); }
  const Teuchos::ParameterList& XFEMGeneralParams() const { return getParameterList()->sublist("XFEM GENERAL"); }
  const Teuchos::ParameterList& XFluidDynamicParams() const { return getParameterList()->sublist("XFLUID DYNAMIC"); }
  const Teuchos::ParameterList& LOMAControlParams() const { return getParameterList()->sublist("LOMA CONTROL"); }
  const Teuchos::ParameterList& BIOFILMControlParams() const { return getParameterList()->sublist("BIOFILM CONTROL"); }
  const Teuchos::ParameterList& ELCHControlParams() const { return getParameterList()->sublist("ELCH CONTROL"); }
  const Teuchos::ParameterList& EPControlParams() const { return getParameterList()->sublist("CARDIAC MONODOMAIN CONTROL"); }
  const Teuchos::ParameterList& ArterialDynamicParams() const { return getParameterList()->sublist("ARTERIAL DYNAMIC"); }
  const Teuchos::ParameterList& ReducedDAirwayDynamicParams() const { return getParameterList()->sublist("REDUCED DIMENSIONAL AIRWAYS DYNAMIC"); }
  const Teuchos::ParameterList& RedAirwayTissueDynamicParams() const { return getParameterList()->sublist("COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC"); }
  const Teuchos::ParameterList& PoroelastDynamicParams() const { return getParameterList()->sublist("POROELASTICITY DYNAMIC"); }
  const Teuchos::ParameterList& PoroFluidMultiPhaseDynamicParams() const { return getParameterList()->sublist("POROFLUIDMULTIPHASE DYNAMIC"); }
  const Teuchos::ParameterList& PoroMultiPhaseScatraDynamicParams() const { return getParameterList()->sublist("POROMULTIPHASESCATRA DYNAMIC"); }
  const Teuchos::ParameterList& PoroMultiPhaseDynamicParams() const { return getParameterList()->sublist("POROMULTIPHASE DYNAMIC"); }
  const Teuchos::ParameterList& PoroScatraControlParams() const { return getParameterList()->sublist("POROSCATRA CONTROL"); }
  const Teuchos::ParameterList& ElastoHydroDynamicParams() const { return getParameterList()->sublist("ELASTO HYDRO DYNAMIC"); }
  const Teuchos::ParameterList& SSIControlParams() const { return getParameterList()->sublist("SSI CONTROL"); }
  const Teuchos::ParameterList& OptimizationControlParams() const { return getParameterList()->sublist("TOPOLOGY OPTIMIZATION CONTROL"); }
  const Teuchos::ParameterList& SearchtreeParams() const { return getParameterList()->sublist("SEARCH TREE"); }
  const Teuchos::ParameterList& StructuralNoxParams() const { return getParameterList()->sublist("STRUCT NOX"); }
  const Teuchos::ParameterList& LocaParams() const { return getParameterList()->sublist("LOCA"); }
  const Teuchos::ParameterList& CavitationParams() const { return getParameterList()->sublist("CAVITATION DYNAMIC"); }
  const Teuchos::ParameterList& ParticleParams() const { return getParameterList()->sublist("PARTICLE DYNAMIC"); }
  const Teuchos::ParameterList& LevelSetControl() const { return getParameterList()->sublist("LEVEL-SET CONTROL"); }
  const Teuchos::ParameterList& WearParams() const { return getParameterList()->sublist("WEAR"); }
  const Teuchos::ParameterList& TSIContactParams() const { return getParameterList()->sublist("TSI CONTACT"); }
  const Teuchos::ParameterList& BeamContactParams() const { return getParameterList()->sublist("BEAM CONTACT"); }
  const Teuchos::ParameterList& BeamPotentialParams() const { return getParameterList()->sublist("BEAM POTENTIAL"); }
  const Teuchos::ParameterList& SemiSmoothPlastParams() const { return getParameterList()->sublist("SEMI-SMOOTH PLASTICITY"); }
  const Teuchos::ParameterList& AcousticParams() const { return getParameterList()->sublist("ACOUSTIC DYNAMIC"); }
  const Teuchos::ParameterList& VolmortarParams() const { return getParameterList()->sublist("VOLMORTAR COUPLING"); }
  const Teuchos::ParameterList& NonlinearSolverParams() const { return getParameterList()()->sublist("NONLINEAR SOLVER"); }
  const Teuchos::ParameterList& TutorialParams()const { return getParameterList()()->sublist("TUTORIAL DYNAMIC"); };

  const Teuchos::ParameterList& SolverParams(int solverNr) const;
  const Teuchos::ParameterList& RandomFieldParams(int randomfieldNr) const;
  const Teuchos::ParameterList& RandomVariableParams(int randomvariableNr) const;

  const Teuchos::ParameterList& UMFPACKSolverParams();

  //@}

  /// @name Discretizations

  /// get access to a particular discretization
  Teuchos::RCP<DRT::Discretization> GetDis(const std::string name) const;

  /// tell number of known fields
  unsigned NumFields() const { return discretizationmap_.size(); }

  /// tell names of known fields
  std::vector<std::string> GetDisNames() const;

  /// check whether a certain discretization exists or not
  bool DoesExistDis(const std::string name) const;

  /// add a discretization to the global problem
  void AddDis(const std::string name, Teuchos::RCP<Discretization> dis);


  //@}

  /// @name Materials

  /// return pointer to materials bundled to the problem
  Teuchos::RCP<MAT::PAR::Bundle> Materials() { return materials_; }

  //@}

  std::map<std::pair<std::string,std::string>,std::map<int,int> > CloningMaterialMap()
      { return clonefieldmatmap_; }

  /// @name Spatial Functions

  /// return a specified function
  DRT::UTILS::Function& Funct(int num);

  //@}

  /// @name Spatial Functions

  /// return a specified time curve
  DRT::UTILS::TimeCurve& Curve(int num);

  //@}

  /// @name Result Tests

  /// Do the testing
  void TestAll(const Epetra_Comm& comm) { resulttest_.TestAll(comm); }

  /// add field specific result test object
  void AddFieldTest(Teuchos::RCP<ResultTest> test) { resulttest_.AddFieldTest(test); }

  //@}

  /// Return the class that handles random numbers globally
  DRT::UTILS::Random* Random() { return &random_; }

  /// Return the class that handles restart initiating -> to be extended
  DRT::UTILS::RestartManager* RestartManager() { return &restartmanager_; }

private:

  /// private default constructor to disallow creation of instances
  Problem();

  /// no copying
  Problem(const Problem&);

  /// no copying
  Problem& operator=(const Problem&);

  /// @name Private access to parameters from certain input file sections

  /// for internal use only! Use Restart(), SpatialApproximation() etc. instead !!!
  const Teuchos::ParameterList& ProblemTypeParams() const { return getParameterList()->sublist("PROBLEM TYP"); }

  /// for internal use only! Use NDim() to access the number of space dimensions
  const Teuchos::ParameterList& ProblemSizeParams() const { return getParameterList()->sublist("PROBLEM SIZE"); }
  //@}

  /// the single instance
  static std::vector<Problem* > instances_;

  /// the problem type
  PROBLEM_TYP probtype_;

  /// the restart step (given by command line or input file)
  int restartstep_;

  /// the restart time (given by command line or input file)
  double restarttime_;

  /// discretizations of this problem
  std::map<std::string,Teuchos::RCP<Discretization> > discretizationmap_;

  /// keep a list of all active element implementations
  std::vector<DRT::SingletonDestruction *> sds_;

  /// material bundle
  Teuchos::RCP<MAT::PAR::Bundle> materials_;

  /// basket of spatial function
  DRT::UTILS::FunctionManager functionmanager_;

  /// basket of time curves
  DRT::UTILS::TimeCurveManager timecurvemanager_;

  /// all test values we might have
  DRT::ResultTestManager resulttest_;

  /// map of coupled fields and corresponding material IDs (needed for cloning
  /// of discretizations)
  std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap_;

  /// nested parallelism group
  Teuchos::RCP<COMM_UTILS::NestedParGroup> npgroup_;

  /// @name File IO

  Teuchos::RCP<IO::InputControl>  inputcontrol_;
  Teuchos::RCP<IO::OutputControl> outputcontrol_;
  Teuchos::RCP<IO::ErrorFileControl> errorfilecontrol_;  ///< error file handles

  //@}

  /// handles all sorts of random numbers
  DRT::UTILS::Random random_;

  /// handles restart
  DRT::UTILS::RestartManager restartmanager_;

};

}

#endif
