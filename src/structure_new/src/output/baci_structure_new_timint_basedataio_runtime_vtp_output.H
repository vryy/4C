/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTP output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_BASEDATAIO_RUNTIME_VTP_OUTPUT_H
#define BACI_STRUCTURE_NEW_TIMINT_BASEDATAIO_RUNTIME_VTP_OUTPUT_H


#include "baci_config.H"

#include "baci_inpar_IO_runtime_vtp_output_structure.H"
#include "baci_io_visualization_parameters.H"

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class BeamRuntimeVtpOutputParams;
  }
}  // namespace DRT


namespace STR
{
  namespace TIMINT
  {
    /** \brief Input data container for VTP output at runtime for the structural (time) integration
     *
     * \author Jonas Eichinger */
    class ParamsRuntimeVtpOutput
    {
     public:
      /// destructor
      virtual ~ParamsRuntimeVtpOutput() = default;

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtp_structure_paramslist);

      /// setup new class variables
      void Setup();

      /**
       * \brief Return the container holding the general output parameters
       */
      const IO::VisualizationParameters& GetVisualizationParameters() const
      {
        return visualization_parameters_;
      }

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int OutputIntervalInSteps() const
      {
        CheckInitSetup();
        return output_interval_steps_;
      };

      [[nodiscard]] int OutputStepOffset() const
      {
        CheckInitSetup();
        return output_step_offset_;
      }

      /// whether to write output in every iteration of the nonlinear solver
      bool OutputEveryIteration() const
      {
        CheckInitSetup();
        return output_every_iteration_;
      };

      /// whether to write owner at visualization point
      bool OutputOwner() const
      {
        CheckInitSetup();
        return output_owner_;
      };

      /// whether to write orientation at visualization point
      bool OutputOrientationAndLength() const
      {
        CheckInitSetup();
        return output_orientationandlength_;
      };

      /// whether to write number of bonds at visualization point
      bool OutputNumberOfBonds() const
      {
        CheckInitSetup();
        return output_numberofbonds_;
      };

      /// whether to write number of bonds at visualization point
      bool OutputLinkingForce() const
      {
        CheckInitSetup();
        return output_linkingforce_;
      };
      /*    /// whether to write displacements
          bool OutputDisplacementState() const
          {
            CheckInitSetup();
            return output_displacement_state_;
          };*/

      /*    /// get the data container for parameters regarding beams
          Teuchos::RCP<const DRT::ELEMENTS::BeamRuntimeVtuOutputParams> GetBeamParams() const
          {
            CheckInitSetup();
            return params_runtime_vtu_output_beams_;
          };*/


     private:
      /// get the init indicator status
      const bool& IsInit() const { return isinit_; };

      /// get the setup indicator status
      const bool& IsSetup() const { return issetup_; };

      /// Check if Init() and Setup() have been called, yet.
      void CheckInitSetup() const;


     private:
      /// @name variables for internal use only
      /// @{
      ///
      bool isinit_ = false;

      bool issetup_ = false;
      /// @}

      /// @name variables controlling output
      /// @{

      //! General visualization parameters
      IO::VisualizationParameters visualization_parameters_;

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int output_interval_steps_ = -1;

      /// An offset added to the current step to shift the steps to be written
      int output_step_offset_ = 0;

      /// whether to write output in every iteration of the nonlinear solver
      bool output_every_iteration_ = false;

      /// whether to write owner at visualization point
      bool output_owner_ = false;

      /// whether to write orientation at visualization point
      bool output_orientationandlength_ = false;

      /// whether to write number of bonds at visualization point
      bool output_numberofbonds_ = false;

      /// whether to write force acting on linker
      bool output_linkingforce_ = false;

      /*    /// whether to write displacement output
          bool output_displacement_state_;*/

      /*    /// data container for input parameters related to VTU output of beams at runtime
          Teuchos::RCP<DRT::ELEMENTS::BeamRuntimeVtuOutputParams>
         params_runtime_vtu_output_beams_;*/


      //@}
    };

  }  // namespace TIMINT
}  // namespace STR

BACI_NAMESPACE_CLOSE

#endif
