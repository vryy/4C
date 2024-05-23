/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTP output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATAIO_RUNTIME_VTP_OUTPUT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATAIO_RUNTIME_VTP_OUTPUT_HPP


#include "4C_config.hpp"

#include "4C_inpar_IO_runtime_vtp_output_structure.hpp"
#include "4C_io_visualization_parameters.hpp"

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN


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

      /// whether to write owner at visualization point
      bool OutputOwner() const
      {
        CheckInitSetup();
        return output_owner_;
      };

      /// whether to write orientation at visualization point
      bool output_orientation_and_length() const
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
          bool output_displacement_state() const
          {
            CheckInitSetup();
            return output_displacement_state_;
          };*/

      /*    /// get the data container for parameters regarding beams
          Teuchos::RCP<const DRT::ELEMENTS::BeamRuntimeOutputParams> GetBeamParams() const
          {
            CheckInitSetup();
            return params_runtime_output_beams_;
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

      /// output interval regarding steps: write output every INTERVAL_STEPS steps
      int output_interval_steps_ = -1;

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

      /*    /// data container for input parameters related to output of beams at runtime
          Teuchos::RCP<DRT::ELEMENTS::BeamRuntimeOutputParams>
         params_runtime_output_beams_;*/


      //@}
    };

  }  // namespace TIMINT
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
