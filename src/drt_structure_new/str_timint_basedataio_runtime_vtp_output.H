/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTP output at runtime for the structural (time) integration

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef STR_TIMINT_BASEDATAIO_RUNTIME_VTP_OUTPUT_H_
#define STR_TIMINT_BASEDATAIO_RUNTIME_VTP_OUTPUT_H_


#include "../drt_inpar/inpar_IO_runtime_vtp_output_structure.H"

namespace Teuchos
{
  class ParameterList;
}

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
      /// constructor
      ParamsRuntimeVtpOutput();

      /// destructor
      virtual ~ParamsRuntimeVtpOutput(){};

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtp_structure_paramslist);

      /// setup new class variables
      void Setup();


      /// whether to write numeric output data binary-encoded or in ascii format
      bool WriteBinaryOutput() const
      {
        CheckInitSetup();
        return (output_data_format_ == INPAR::IO_RUNTIME_VTP_STRUCTURE::binary);
      }

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int OutputIntervalInSteps() const
      {
        CheckInitSetup();
        return output_interval_steps_;
      };

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
      bool isinit_;

      bool issetup_;
      /// @}

      /// @name variables controlling output
      /// @{

      /// data format for numeric output data
      INPAR::IO_RUNTIME_VTP_STRUCTURE::OutputDataFormat output_data_format_;

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int output_interval_steps_;

      /// whether to write output in every iteration of the nonlinear solver
      bool output_every_iteration_;

      /// whether to write owner at visualization point
      bool output_owner_;

      /// whether to write orientation at visualization point
      bool output_orientationandlength_;

      /// whether to write number of bonds at visualization point
      bool output_numberofbonds_;

      /// whether to write force acting on linker
      bool output_linkingforce_;

      /*    /// whether to write displacement output
          bool output_displacement_state_;*/

      /*    /// data container for input parameters related to VTU output of beams at runtime
          Teuchos::RCP<DRT::ELEMENTS::BeamRuntimeVtuOutputParams>
         params_runtime_vtu_output_beams_;*/


      //@}
    };

  }  // namespace TIMINT
}  // namespace STR

#endif
