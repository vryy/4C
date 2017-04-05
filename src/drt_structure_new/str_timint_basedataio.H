/*-----------------------------------------------------------*/
/*!
\file str_timint_basedataio.H

\brief Input/output data container for the structural (time)
       integration

\maintainer Michael Hiermeier

\date Jan 11, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_H_

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_structure.H"

namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace IO {
  class DiscretizationWriter;
} // namespace IO
namespace STR {
namespace TIMINT {
  class ParamsRuntimeVtkOutput;

  /** \brief Input/output data container for the structural (time) integration
   *
   * This data container holds everything, which refers directly to the
   * input/output writer and the screen output.
   *
   * \author Michael Hiermeier */
  class BaseDataIO
  {
  public:
    /// constructor
    BaseDataIO();

    /// destructor
    virtual ~BaseDataIO() {};

    /// initialize the class variables
    void Init(const Teuchos::ParameterList& IOParams,
        const Teuchos::ParameterList& sDynParams,
        const Teuchos::ParameterList& xParams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// setup new class variables
    void Setup();

  protected:
    /// get the init indicator status
    virtual const bool& IsInit() const { return isinit_; };

    /// get the setup indicator status
    virtual const bool& IsSetup() const { return issetup_; };

    /// Check if Init() and Setup() have been called, yet.
    virtual void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

  public:
    /// get the binary output writer
    Teuchos::RCP<IO::DiscretizationWriter> GetMutableOutputPtr()
    {
      CheckInitSetup();
      return output_;
    };

    /// get the data container for parameters regarding VTK output at runtime
    Teuchos::RCP<const ParamsRuntimeVtkOutput> GetRuntimeVtkOutputParams() const
    {
      CheckInitSetup();
      return params_runtime_vtk_output_;
    };

    /// @name Get printing and output parameters/files (read only access)
    ///@{
    /// get the output file for energy (read only access)
    Teuchos::RCP<const std::ostream> GetEnergyFilePtr() const
    {
      CheckInitSetup();
      return energyfile_;
    };

    /// get error file pointer
    FILE* ErrorFilePtr()
    {
      CheckInitSetup();
      return errfile_;
    };

    /// Is GMSH output of displacements required?
    const bool& IsGmsh() const
    {
      CheckInitSetup();
      return gmsh_out_;
    };

    /// Shall we print the logo?
    const bool& IsLogo() const
    {
      CheckInitSetup();
      return printlogo_;
    };

    /// Shall we print the infos to errorfile #errfile_?
    const bool& IsErrorFile() const
    {
      CheckInitSetup();
      return printerrfile_;
    };

    /// Shall we print intermediate iterations during solution?
    const bool& IsPrintIntermediateIterations() const
    {
      CheckInitSetup();
      return printiter_;
    };

    /// Shall we write output every iteration?
    const bool& IsOutputEveryIter() const
    {
      CheckInitSetup();
      return outputeveryiter_;
    };

    /// Shall we write surfactant output?
    const bool& IsSurfactantOutput() const
    {
      CheckInitSetup();
      return writesurfactant_;
    };

    /// Shall we write the current state?
    const bool& IsWriteState() const
    {
      CheckInitSetup();
      return writestate_;
    };

    /// Shall we write the velocities and accelerations?
    const bool& IsWriteVelAcc() const
    {
      CheckInitSetup();
      return writevelacc_;
    };

    /// Is this the first output of the current run?
    const bool& IsFirstOutputOfRun() const
    {
      CheckInitSetup();
      return firstoutputofrun_;
    };

    /// Print infos to standard out every n step
    const int& GetPrint2ScreenEveryNStep() const
    {
      CheckInitSetup();
      return printscreen_;
    };

    /// Get the user defined filename counter (necessary for OutputEveryIter)
    const int& GetOEI_FileCounter() const
    {
      CheckInitSetup();
      return oei_filecounter_;
    };

    /// Get the output counter for OutputEveryIter
    const int& GetOEI_OutputCounter() const
    {
      CheckInitSetup();
      return outputcounter_;
    };

    /// write restart every given step. if 0, restart is not written
    const int& GetWriteRestartEveryNStep() const
    {
      CheckInitSetup();
      return writerestartevery_;
    };

    /// Write reduced Restart Data i.e no mesh for mlmc purposes
    const int& GetWriteReducedRestartEveryNStep() const
    {
      CheckInitSetup();
      return writereducedrestart_;
    };

    /// write state/stress/strain every given step
    const int& GetWriteResultsEveryNStep() const
    {
      CheckInitSetup();
      return writeresultsevery_;
    };

    /// write system energy every given step
    const int& GetWriteEnergyEveryNStep() const
    {
      CheckInitSetup();
      return writeenergyevery_;
    }

    /// get kinetic energy of system
    const double& GetKineticEnergy() const
    {
      CheckInitSetup();
      return kinergy_;
    };

    /// get internal energy of system
    const double& GetInternalEnergy() const
    {
      CheckInitSetup();
      return intergy_;
    }

    /// get external energy of system
    const double& GetExternalEnergy() const
    {
      CheckInitSetup();
      return extergy_;
    }

    /// get stress output type
    const INPAR::STR::StressType& GetStressOutputType() const
    {
      CheckInitSetup();
      return writestress_;
    }

    /// get output type of coupling stress
    const INPAR::STR::StressType& GetCouplingStressOutputType() const
    {
      CheckInitSetup();
      return writecouplstress_;
    }

    /// get strain output type
    const INPAR::STR::StrainType& GetStrainOutputType() const
    {
      CheckInitSetup();
      return writestrain_;
    }

    /// get plastic strain output type
    const INPAR::STR::StrainType& GetPlasticStrainOutputType() const
    {
      CheckInitSetup();
      return writeplstrain_;
    };
    ///@}

    /// set the flag indicator firstoutputofrun_
    void SetFirstOutputOfRun(const bool& firstoutputofrun)
    {
      CheckInitSetup();
      firstoutputofrun_ = firstoutputofrun;
    }

  protected:
    /// @name variables for internal use only
    /// @{
    ///
    bool isinit_;

    bool issetup_;
    /// @}

  private:
    /// @name Printing and output
    /// @{

    /// binary output
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    /// data container for input parameters related to VTK output at runtime
    Teuchos::RCP<ParamsRuntimeVtkOutput> params_runtime_vtk_output_;

    /// outputfile for energy
    Teuchos::RCP<std::ofstream> energyfile_;

    /// error file handle
    FILE* errfile_;

    /// Is GMSH output of displacements required?
    bool gmsh_out_;

    /// print the logo (or not)?
    bool printlogo_;

    /// print infos to error file #errfile_
    bool printerrfile_;

    /// print intermediate iterations during solution
    bool printiter_;

    /// write output every iteration (Newton, line search, load step, etc.)
    bool outputeveryiter_;

    /// write surfactant output
    bool writesurfactant_;

    /// write state on/off
    bool writestate_;

    /// write velocity and acceleration on/off
    bool writevelacc_;

    /// flag whether this output step is the first one (restarted or not)
    bool firstoutputofrun_;

    /// print infos to standard out every n steps
    int printscreen_;

    /// filename counter
    int oei_filecounter_;

    /// output counter for OutputEveryIter
    int outputcounter_;

    /// write restart every given step. if 0, restart is not written
    int writerestartevery_;

    /// Write reduced Restart Data i.e no mesh for mlmc purposes
    int writereducedrestart_;

    /// write state/stress/strain every given step
    int writeresultsevery_;

    /// write system energy every given step
    int writeenergyevery_;

    /// kinetic energy of system
    double kinergy_;

    /// internal energy of system
    double intergy_;

    /// external energy of system
    double extergy_;

    /// stress output type
    INPAR::STR::StressType writestress_;

    /// output type of coupling stress
    INPAR::STR::StressType writecouplstress_;

    /// strain output type
    INPAR::STR::StrainType writestrain_;

    /// plastic strain output type
    INPAR::STR::StrainType writeplstrain_;

    //@}
  };  // class BaseDataIO
} // namespace TIMINT
} // namespace STR



#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_H_ */
