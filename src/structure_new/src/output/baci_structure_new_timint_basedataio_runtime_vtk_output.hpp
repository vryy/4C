/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTK output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATAIO_RUNTIME_VTK_OUTPUT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATAIO_RUNTIME_VTK_OUTPUT_HPP

#include "baci_config.hpp"

#include "baci_io_visualization_parameters.hpp"

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class StructureRuntimeOutputParams;
    class BeamRuntimeOutputParams;
  }  // namespace ELEMENTS
}  // namespace DRT


namespace STR
{
  namespace TIMINT
  {
    /** \brief Input data container for output at runtime for the structural (time) integration
     *
     * \author Maximilian Grill */
    class ParamsRuntimeOutput
    {
     public:
      /// destructor
      virtual ~ParamsRuntimeOutput() = default;

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtk_structure_paramslist);

      /// setup new class variables
      void Setup();

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

      /// whether to write special output for structure elements
      bool OutputStructure() const
      {
        CheckInitSetup();
        return output_structure_;
      };

      /// whether to write special output for structure elements
      bool OutputBeams() const
      {
        CheckInitSetup();
        return output_beams_;
      };

      /// get the data container for parameters regarding beams
      Teuchos::RCP<const DRT::ELEMENTS::StructureRuntimeOutputParams> GetStructureParams() const
      {
        CheckInitSetup();
        return params_runtime_output_structure_;
      };

      /// get the data container for parameters regarding beams
      Teuchos::RCP<const DRT::ELEMENTS::BeamRuntimeOutputParams> GetBeamParams() const
      {
        CheckInitSetup();
        return params_runtime_output_beams_;
      };


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

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int output_interval_steps_ = -1;

      /// An offset added to the current step to shift the steps to be written
      int output_step_offset_ = 0;

      /// whether to write output in every iteration of the nonlinear solver
      bool output_every_iteration_ = false;

      /// whether to write output for structural elements
      bool output_structure_ = false;

      /// whether to write special output for beam elements
      bool output_beams_ = false;

      /// data container for input parameters related to output of structure at runtime
      Teuchos::RCP<DRT::ELEMENTS::StructureRuntimeOutputParams> params_runtime_output_structure_ =
          Teuchos::null;

      /// data container for input parameters related to output of beams at runtime
      Teuchos::RCP<DRT::ELEMENTS::BeamRuntimeOutputParams> params_runtime_output_beams_ =
          Teuchos::null;


      //@}
    };

  }  // namespace TIMINT
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif
