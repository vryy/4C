/*----------------------------------------------------------------------*/
/*! \file

\brief internal representation of an input file

\level 1


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_GLOBAL_DATA_HPP
#define FOUR_C_GLOBAL_DATA_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_shape_function_type.hpp"
#include "4C_global_data_enums.hpp"
#include "4C_lib_utils_restart_manager.hpp"
#include "4C_utils_function_manager.hpp"
#include "4C_utils_random.hpp"
#include "4C_utils_result_test.hpp"
#include "4C_utils_std_cxx20_ranges.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_RCP.hpp>

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace IO
{
  class OutputControl;
  class InputControl;
}  // namespace IO

namespace CORE::COMM
{
  class Communicators;
}

namespace MAT
{
  namespace PAR
  {
    class Bundle;
  }
}  // namespace MAT

namespace PARTICLEENGINE
{
  class ParticleObject;
}

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class Bundle;
  }
}  // namespace CONTACT

namespace INPUT
{
  class DatFileReader;
}
namespace DRT
{
  class Discretization;
}

namespace GLOBAL
{
  /*!
   * The Problem class gathers various input parameters and provides access
   * from anywhere via the singleton Instance() function.
   *
   * This class is an old attempt to deal with parameters. The fundamental problem lies in the
   * global nature of the singleton instance. This behavior makes it very difficult to follow the
   * flow of data in the code and introduces hidden dependencies between different modules.
   * Nowadays, we know better and do not like to write new code that uses this class. Instead, try
   * to pass whatever data is needed direclty to a class or function. We work on removing
   * functionality from this class.
   *
   * Nevertheless, here is the old documentation, for as long as we'll be using this class:
   *
   * Global problem instance that keeps the discretizations
   * The global problem represents the input file passed to 4C. This class organizes the reading of
   * a dat file (utilizing the DatFileReader of course). That is way, in all but the most eccentric
   * cases there will be exactly one object of this class during a 4C run. This object contains all
   * parameters read from the input file as well as any material definitions and even all the
   * discretizations.
   *
   * <h3>Input parameters</h3>
   *
   * All input parameters are known by the global problem object. These parameters are guaranteed
   * to be valid (because they passed the validation) and are guaranteed to be there (because
   * default values have been set for all parameters missing from the input file.) This is
   * Teuchos::ParameterList magic, that just requires the list of valid parameters in the file
   * validparameters.cpp to be complete.
   *
   * The algorithms are meant to ask the global problem object for their parameters and extract
   * them from the respective parameter list directly.
   *
   * <h3>Discretizations</h3>
   *
   * The global problem object knows the discretizations defined by the input file. In order to
   * access a particular discretization you get the global problem object and ask.
   *
   * <h3>Materials</h3>
   *
   * The global problem object knows the material descriptions from the input file. These are not to
   * be confused with the material classes the elements know and work with. The global problem
   * object does not keep track of gauss point material values, all that is known here are the
   * definitions from the input file.
   *
   * <h3>Singleton behavior</h3>
   *
   * \warning This is a guru only section!
   *
   * The global problem behaves like a singleton, so there is always one instance available. But you
   * can have more than one instance of Problem. In normal situations this will not be needed. So
   * don't bother. Just call the static Instance() function to get the global instance and access
   * your discretizations.
   *
   * In the special case that you want to read more that one input file, however, you will need to
   * handle the fields from each file separately.
   *
   * One artefact that comes from using global variables together with multiple Problem objects is
   * the notion of activating of problem object. This translates to setting the global variable
   * pointers to this Problem object's internal variables. Normally, if there is just one Problem
   * object, this is done by default. If you need more that one, however, you will have to activate
   * the global problems yourself.
   */
  class Problem
  {
   public:
    /// @name Instances

    /// Disallow copying this class.
    Problem(const Problem&) = delete;

    /// Disallow copying this class.
    Problem& operator=(const Problem&) = delete;

    /// Disallow moving this class.
    Problem(Problem&&) = delete;

    /// Disallow moving this class.
    Problem& operator=(Problem&&) = delete;

    /// return an instance of this class
    static Problem* Instance(int num = 0);

    /// return number of problem instances
    static unsigned NumInstances() { return instances_.size(); }

    /// calculation done, clean up
    /*!
      There can be a variety of objects to a problem. Some of them might
      require proper cleanup. Make sure we always do it.
     */
    static void Done();

    //@}

    /// @name Input

    /// set restart step which was read from the command line
    void SetRestartStep(int r);

    void SetInputControlFile(Teuchos::RCP<IO::InputControl>& input) { inputcontrol_ = input; }

    /// manipulate problem type
    void SetProblemType(GLOBAL::ProblemType targettype);

    void SetSpatialApproximationType(CORE::FE::ShapeFunctionType shape_function_type);

    /// @name General query methods
    /// Once and for all definitions

    /// give enum of my problem type
    GLOBAL::ProblemType GetProblemType() const;

    /// give string name of my problem type
    std::string ProblemName() const;

    /// return restart step
    [[nodiscard]] int Restart() const;

    /// number of space dimensions (as specified in the input file)
    int NDim() const;

    /// current wall time in seconds with micro second precision
    static double Walltime();


    //! Return type of the basis function encoded as enum
    CORE::FE::ShapeFunctionType SpatialApproximationType() const { return shapefuntype_; }

    //! @}

    /// @name Control file

    /*!
    \brief Create control file for output and read restart data if required

    In addition, issue a warning to the screen, if no binary output will be written.

    @param[in] comm Communicator
    @param[in] inputfile File name of input file
    @param[in] prefix
    @param[in] restartkenner
    */
    void OpenControlFile(const Epetra_Comm& comm, const std::string& inputfile, std::string prefix,
        const std::string& restartkenner);

    /// control file for restart read
    Teuchos::RCP<IO::InputControl> InputControlFile() { return inputcontrol_; }

    /// control file for normal output
    Teuchos::RCP<IO::OutputControl> OutputControlFile() { return outputcontrol_; }

    /// write parameters read from input file for documentation
    void WriteInputParameters();

    //@}

    /// @name Parameters read from file

    /// Set parameters from a parameter list and return with default values.
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& parameter_list);

    /// Return a const parameter list of all of the valid parameters that
    /// this->setParameterList(...) will accept.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

    Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

    /// @name Communicators and their parallel groups

    /// set communicators
    void SetCommunicators(Teuchos::RCP<CORE::COMM::Communicators> communicators);

    /// return communicators
    Teuchos::RCP<CORE::COMM::Communicators> GetCommunicators() const;

    //@}

    /// @name Input parameter sections
    /// direct access to parameters from input file sections

    const Teuchos::ParameterList& BinningStrategyParams() const
    {
      return parameters_->sublist("BINNING STRATEGY");
    }
    const Teuchos::ParameterList& GeometricSearchParams() const
    {
      return parameters_->sublist("BOUNDINGVOLUME STRATEGY");
    }
    const Teuchos::ParameterList& IOParams() const { return parameters_->sublist("IO"); }
    const Teuchos::ParameterList& StructuralDynamicParams() const
    {
      return parameters_->sublist("STRUCTURAL DYNAMIC");
    }
    const Teuchos::ParameterList& Cardiovascular0DStructuralParams() const
    {
      return parameters_->sublist("CARDIOVASCULAR 0D-STRUCTURE COUPLING");
    }
    const Teuchos::ParameterList& MortarCouplingParams() const
    {
      return parameters_->sublist("MORTAR COUPLING");
    }
    const Teuchos::ParameterList& ContactDynamicParams() const
    {
      return parameters_->sublist("CONTACT DYNAMIC");
    }
    const Teuchos::ParameterList& BeamInteractionParams() const
    {
      return parameters_->sublist("BEAM INTERACTION");
    }
    const Teuchos::ParameterList& RveMultiPointConstraintParams() const
    {
      return getParameterList()->sublist("MULTI POINT CONSTRAINTS");
    }
    const Teuchos::ParameterList& BrownianDynamicsParams() const
    {
      return parameters_->sublist("BROWNIAN DYNAMICS");
    }
    const Teuchos::ParameterList& ThermalDynamicParams() const
    {
      return parameters_->sublist("THERMAL DYNAMIC");
    }
    const Teuchos::ParameterList& TSIDynamicParams() const
    {
      return parameters_->sublist("TSI DYNAMIC");
    }
    const Teuchos::ParameterList& FluidDynamicParams() const
    {
      return parameters_->sublist("FLUID DYNAMIC");
    }
    const Teuchos::ParameterList& LubricationDynamicParams() const
    {
      return parameters_->sublist("LUBRICATION DYNAMIC");
    }
    const Teuchos::ParameterList& ScalarTransportDynamicParams() const
    {
      return parameters_->sublist("SCALAR TRANSPORT DYNAMIC");
    }
    const Teuchos::ParameterList& STIDynamicParams() const
    {
      return parameters_->sublist("STI DYNAMIC");
    }
    const Teuchos::ParameterList& FS3IDynamicParams() const
    {
      return parameters_->sublist("FS3I DYNAMIC");
    }
    const Teuchos::ParameterList& AleDynamicParams() const
    {
      return parameters_->sublist("ALE DYNAMIC");
    }
    const Teuchos::ParameterList& FSIDynamicParams() const
    {
      return parameters_->sublist("FSI DYNAMIC");
    }
    const Teuchos::ParameterList& FPSIDynamicParams() const
    {
      return parameters_->sublist("FPSI DYNAMIC");
    }
    const Teuchos::ParameterList& ImmersedMethodParams() const
    {
      return parameters_->sublist("IMMERSED METHOD");
    }
    const Teuchos::ParameterList& CutGeneralParams() const
    {
      return parameters_->sublist("CUT GENERAL");
    }
    const Teuchos::ParameterList& XFEMGeneralParams() const
    {
      return parameters_->sublist("XFEM GENERAL");
    }
    const Teuchos::ParameterList& XFluidDynamicParams() const
    {
      return parameters_->sublist("XFLUID DYNAMIC");
    }
    const Teuchos::ParameterList& FBIParams() const
    {
      return parameters_->sublist("FLUID BEAM INTERACTION");
    }
    const Teuchos::ParameterList& LOMAControlParams() const
    {
      return parameters_->sublist("LOMA CONTROL");
    }
    const Teuchos::ParameterList& BIOFILMControlParams() const
    {
      return parameters_->sublist("BIOFILM CONTROL");
    }
    const Teuchos::ParameterList& ELCHControlParams() const
    {
      return parameters_->sublist("ELCH CONTROL");
    }
    const Teuchos::ParameterList& EPControlParams() const
    {
      return parameters_->sublist("CARDIAC MONODOMAIN CONTROL");
    }
    const Teuchos::ParameterList& ArterialDynamicParams() const
    {
      return parameters_->sublist("ARTERIAL DYNAMIC");
    }
    const Teuchos::ParameterList& ReducedDAirwayDynamicParams() const
    {
      return parameters_->sublist("REDUCED DIMENSIONAL AIRWAYS DYNAMIC");
    }
    const Teuchos::ParameterList& RedAirwayTissueDynamicParams() const
    {
      return parameters_->sublist("COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC");
    }
    const Teuchos::ParameterList& PoroelastDynamicParams() const
    {
      return parameters_->sublist("POROELASTICITY DYNAMIC");
    }
    const Teuchos::ParameterList& PoroFluidMultiPhaseDynamicParams() const
    {
      return parameters_->sublist("POROFLUIDMULTIPHASE DYNAMIC");
    }
    const Teuchos::ParameterList& PoroMultiPhaseScatraDynamicParams() const
    {
      return parameters_->sublist("POROMULTIPHASESCATRA DYNAMIC");
    }
    const Teuchos::ParameterList& PoroMultiPhaseDynamicParams() const
    {
      return parameters_->sublist("POROMULTIPHASE DYNAMIC");
    }
    const Teuchos::ParameterList& PoroScatraControlParams() const
    {
      return parameters_->sublist("POROSCATRA CONTROL");
    }
    const Teuchos::ParameterList& ElastoHydroDynamicParams() const
    {
      return parameters_->sublist("ELASTO HYDRO DYNAMIC");
    }
    const Teuchos::ParameterList& SSIControlParams() const
    {
      return parameters_->sublist("SSI CONTROL");
    }
    const Teuchos::ParameterList& SSTIControlParams() const
    {
      return parameters_->sublist("SSTI CONTROL");
    }
    const Teuchos::ParameterList& SearchtreeParams() const
    {
      return parameters_->sublist("SEARCH TREE");
    }
    const Teuchos::ParameterList& StructuralNoxParams() const
    {
      return parameters_->sublist("STRUCT NOX");
    }
    const Teuchos::ParameterList& LocaParams() const { return parameters_->sublist("LOCA"); }
    const Teuchos::ParameterList& ParticleParams() const
    {
      return parameters_->sublist("PARTICLE DYNAMIC");
    }
    const Teuchos::ParameterList& PASIDynamicParams() const
    {
      return parameters_->sublist("PASI DYNAMIC");
    }
    const Teuchos::ParameterList& LevelSetControl() const
    {
      return parameters_->sublist("LEVEL-SET CONTROL");
    }
    const Teuchos::ParameterList& WearParams() const { return parameters_->sublist("WEAR"); }
    const Teuchos::ParameterList& TSIContactParams() const
    {
      return parameters_->sublist("TSI CONTACT");
    }
    const Teuchos::ParameterList& BeamContactParams() const
    {
      return parameters_->sublist("BEAM CONTACT");
    }
    const Teuchos::ParameterList& BeamPotentialParams() const
    {
      return parameters_->sublist("BEAM POTENTIAL");
    }
    const Teuchos::ParameterList& SemiSmoothPlastParams() const
    {
      return parameters_->sublist("SEMI-SMOOTH PLASTICITY");
    }
    const Teuchos::ParameterList& ElectromagneticParams() const
    {
      return parameters_->sublist("ELECTROMAGNETIC DYNAMIC");
    }
    const Teuchos::ParameterList& VolmortarParams() const
    {
      return parameters_->sublist("VOLMORTAR COUPLING");
    }
    const Teuchos::ParameterList& MORParams() const { return parameters_->sublist("MOR"); };
    const Teuchos::ParameterList& MeshPartitioningParams() const
    {
      return parameters_->sublist("MESH PARTITIONING");
    }

    const Teuchos::ParameterList& ProblemTypeParams() const
    {
      return parameters_->sublist("PROBLEM TYP");
    }

    const Teuchos::ParameterList& ProblemSizeParams() const
    {
      return parameters_->sublist("PROBLEM SIZE");
    }

    const Teuchos::ParameterList& SolverParams(int solverNr) const;

    const Teuchos::ParameterList& UMFPACKSolverParams();

    //@}

    /// @name Discretizations

    /// get access to a particular discretization
    Teuchos::RCP<DRT::Discretization> GetDis(const std::string& name) const;

    auto DiscretizationRange() { return std_20::ranges::views::all(discretizationmap_); }

    auto DiscretizationRange() const { return std_20::ranges::views::all(discretizationmap_); }

    /// tell number of known fields
    unsigned NumFields() const { return discretizationmap_.size(); }

    /// tell names of known fields
    std::vector<std::string> GetDisNames() const;

    /// check whether a certain discretization exists or not
    bool DoesExistDis(const std::string& name) const;

    /// add a discretization to the global problem
    void AddDis(const std::string& name, Teuchos::RCP<DRT::Discretization> dis);


    //@}

    /// @name Materials

    /// return pointer to materials bundled to the problem
    Teuchos::RCP<MAT::PAR::Bundle> Materials() { return materials_; }

    // return pointer to contact constitutive law bundled to the problem
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Bundle> ContactConstitutiveLaws()
    {
      return contactconstitutivelaws_;
    }

    //@}

    /// @name Particles

    /// return reference to read in particles
    std::vector<std::shared_ptr<PARTICLEENGINE::ParticleObject>>& Particles() { return particles_; }

    //@}

    std::map<std::pair<std::string, std::string>, std::map<int, int>>& CloningMaterialMap()
    {
      return clonefieldmatmap_;
    }

    /// @name Spatial Functions

    /**
     * Get a function read from the input file by its ID @p num.
     *
     * @tparam T The type of function interface.
     */
    template <typename T>
    const T& FunctionById(int num)
    {
      return functionmanager_.template FunctionById<T>(num);
    }

    //@}

    /// @name Result Tests

    /// Do the testing
    void TestAll(const Epetra_Comm& comm) { resulttest_.TestAll(comm); }

    /// add field specific result test object
    void AddFieldTest(Teuchos::RCP<CORE::UTILS::ResultTest> test)
    {
      resulttest_.AddFieldTest(test);
    }

    CORE::UTILS::ResultTestManager& GetResultTestManager() { return resulttest_; }

    //@}

    /// Return the class that handles random numbers globally
    CORE::UTILS::Random* Random() { return &random_; }

    /// Return the class that handles restart initiating -> to be extended
    DRT::UTILS::RestartManager* RestartManager() { return &restartmanager_; }

    /**
     * Set the @p function_manager which contains all parsed functions.
     *
     * @note The parsing of functions must take place before. This calls wants a filled
     * FunctionManager.
     */
    void SetFunctionManager(CORE::UTILS::FunctionManager&& function_manager);

   private:
    /// private default constructor to disallow creation of instances
    Problem();

    /// the single instance
    static std::vector<Problem*> instances_;

    /// the problem type
    GLOBAL::ProblemType probtype_;

    /// Spatial approximation type
    CORE::FE::ShapeFunctionType shapefuntype_;

    /// the restart step (given by command line or input file)
    int restartstep_;

    /// discretizations of this problem
    std::map<std::string, Teuchos::RCP<DRT::Discretization>> discretizationmap_;

    /// material bundle
    Teuchos::RCP<MAT::PAR::Bundle> materials_;

    /// bundle containing all read-in contact constitutive laws
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Bundle> contactconstitutivelaws_;

    /// all particles that are read in
    std::vector<std::shared_ptr<PARTICLEENGINE::ParticleObject>> particles_;

    /// basket of spatial function
    CORE::UTILS::FunctionManager functionmanager_;

    /// all test values we might have
    CORE::UTILS::ResultTestManager resulttest_;

    /// map of coupled fields and corresponding material IDs (needed for cloning
    /// of discretizations)
    std::map<std::pair<std::string, std::string>, std::map<int, int>> clonefieldmatmap_;

    /// communicators
    Teuchos::RCP<CORE::COMM::Communicators> communicators_;

    /// @name File IO

    Teuchos::RCP<IO::InputControl> inputcontrol_;
    Teuchos::RCP<IO::OutputControl> outputcontrol_;

    //@}

    /// handles all sorts of random numbers
    CORE::UTILS::Random random_;

    /// handles restart
    DRT::UTILS::RestartManager restartmanager_;

    //! The central list of all paramters read from input.
    Teuchos::RCP<Teuchos::ParameterList> parameters_;
  };

}  // namespace GLOBAL

FOUR_C_NAMESPACE_CLOSE

#endif
