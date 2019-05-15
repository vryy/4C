/*-----------------------------------------------------------*/
/*!

\brief Input/output data container for the structural (time) integration

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_H_

#include "../drt_inpar/inpar_structure.H"

#include "../solver_nonlin_nox/nox_nln_abstract_prepostoperator.H"

// forward declarations
namespace IO
{
  class DiscretizationWriter;
  class EveryIterationWriterInterface;
  class EveryIterationWriter;
}  // namespace IO

namespace NOX
{
  namespace LineSearch
  {
    class Generic;
  }  // namespace LineSearch
}  // namespace NOX

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace STR
{
  namespace TIMINT
  {
    class ParamsRuntimeVtkOutput;
    class ParamsRuntimeVtpOutput;

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
      virtual ~BaseDataIO(){};

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IOParams, const Teuchos::ParameterList& sDynParams,
          const Teuchos::ParameterList& xParams, Teuchos::RCP<IO::DiscretizationWriter> output);

      /// setup new class variables
      void Setup();

     protected:
      /// get the init indicator status
      virtual const bool& IsInit() const { return isinit_; };

      /// get the setup indicator status
      virtual const bool& IsSetup() const { return issetup_; };

      /// Check if Init() and Setup() have been called, yet.
      virtual void CheckInitSetup() const;

     public:
      /// get the binary output writer
      Teuchos::RCP<IO::DiscretizationWriter> GetMutableOutputPtr()
      {
        CheckInitSetup();
        return output_;
      };

      /// get the binary output writer
      Teuchos::RCP<const IO::DiscretizationWriter> GetOutputPtr() const
      {
        CheckInitSetup();
        return output_;
      }

      /// get the data container for parameters regarding VTK output at runtime
      Teuchos::RCP<const ParamsRuntimeVtkOutput> GetRuntimeVtkOutputParams() const
      {
        CheckInitSetup();
        return params_runtime_vtk_output_;
      };

      /// get the data container for parameters regarding VTK output at runtime
      Teuchos::RCP<const ParamsRuntimeVtpOutput> GetRuntimeVtpOutputParams() const
      {
        CheckInitSetup();
        return params_runtime_vtp_output_;
      };

      /// \brief return TRUE if the results shall be written for this load/time \c step
      bool WriteResultsForThisStep(const int step) const;

      /// @name Get printing and output parameters/files
      ///@{
      /// get the output file for energy
      std::ostream& GetEnergyOutputStream()
      {
        CheckInitSetup();

        if (energyfile_.is_null()) dserror("energy file stream uninitialized");

        return *energyfile_;
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

      /// Shall we write the current element volume?
      bool IsWriteCurrentEleVolume() const
      {
        CheckInitSetup();
        return writecurrentelevolume_;
      }

      /// Shall we write the jacobian to MATLAB?
      bool IsWriteJacobianToMatlab() const
      {
        CheckInitSetup();
        return writejac2matlab_;
      }

      /// Shall we compute and write the condition number?
      INPAR::STR::ConditionNumber ConditionNumberType() const
      {
        CheckInitSetup();
        return conditionnumbertype_;
      }

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

      /// get optional quantity output type
      const INPAR::STR::OptQuantityType& GetOptQuantityOutputType() const
      {
        CheckInitSetup();
        return writeoptquantity_;
      }
      ///@}

      /// set the flag indicator firstoutputofrun_
      void SetFirstOutputOfRun(const bool& firstoutputofrun)
      {
        CheckInitSetup();
        firstoutputofrun_ = firstoutputofrun;
      }

      /// Initialize and setup the every iteration output writer
      void InitSetupEveryIterationWriter(
          IO::EveryIterationWriterInterface* interface, Teuchos::ParameterList& p_nox);

      /// initialize the output of system energy
      void SetupEnergyOutputFile();

     protected:
      /// @name variables for internal use only
      ///@{
      ///
      bool isinit_;

      bool issetup_;
      ///@}

     private:
      /// @name Printing and output
      ///@{

      /// binary output
      Teuchos::RCP<IO::DiscretizationWriter> output_;

      /// additional output writer for the Newton steps
      Teuchos::RCP<IO::EveryIterationWriter> writer_every_iter_;

      /// data container for input parameters related to VTK output at runtime
      Teuchos::RCP<ParamsRuntimeVtkOutput> params_runtime_vtk_output_;

      /// data container for input parameters related to VTP output at runtime
      Teuchos::RCP<ParamsRuntimeVtpOutput> params_runtime_vtp_output_;

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

      /// write jacobian to MATLAB
      bool writejac2matlab_;

      /// flag whether this output step is the first one (restarted or not)
      bool firstoutputofrun_;

      /// flag element volume on/off
      bool writecurrentelevolume_ = false;

      /// print infos to standard out every n steps
      int printscreen_;

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

      /// stress output type
      INPAR::STR::StressType writestress_;

      /// output type of coupling stress
      INPAR::STR::StressType writecouplstress_;

      /// strain output type
      INPAR::STR::StrainType writestrain_;

      /// plastic strain output type
      INPAR::STR::StrainType writeplstrain_;

      /// optional quantity type
      INPAR::STR::OptQuantityType writeoptquantity_;

      INPAR::STR::ConditionNumber conditionnumbertype_;

      Teuchos::RCP<Teuchos::ParameterList> p_io_every_iteration_;

      ///@}
    };  // class BaseDataIO
  }     // namespace TIMINT
}  // namespace STR

namespace NOX
{
  namespace NLN
  {
    namespace Solver
    {
      namespace PrePostOp
      {
        namespace TIMINT
        {
          /*! \brief Helper class to write the output each Newton step
           *
           *  This class is an implementation of the NOX::Abstract::PrePostOperator
           *  and is used to modify the step() routine of the given NOX::Solver::Generic
           *  class.
           *  It's called by the wrapper classes NOX::Solver::PrePostOperator and
           *  NOX::PrePostOperatorVector.
           *
           *  \author Michael Hiermeier \date 03/17 */
          class WriteOutputEveryIteration : public NOX::NLN::Abstract::PrePostOperator
          {
           public:
            /// constructor
            WriteOutputEveryIteration(IO::EveryIterationWriter& every_iter_writer);

            /// destructor
            virtual ~WriteOutputEveryIteration(){};

            /// called at the very beginning of a Newton loop
            virtual void runPreSolve(const NOX::Solver::Generic& solver);

            /// called in the end of each Newton step
            virtual void runPostIterate(const NOX::Solver::Generic& solver);

            /// called before the step is reduced in a line search routine
            virtual void runPreModifyStepLength(
                const NOX::Solver::Generic& solver, const NOX::LineSearch::Generic& linesearch);

           private:
            IO::EveryIterationWriter& every_iter_writer_;
          };
        }  // namespace TIMINT
      }    // namespace PrePostOp
    }      // namespace Solver
  }        // namespace NLN
}  // namespace NOX



#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_H_ */
