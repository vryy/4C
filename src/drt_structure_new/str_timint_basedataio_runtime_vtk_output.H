/*-----------------------------------------------------------------------------------------------*/
/*!
\file str_timint_basedataio_runtime_vtk_output.H

\brief input parameters related to VTK output at runtime for the structural (time) integration

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_RUNTIME_VTK_OUTPUT_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAIO_RUNTIME_VTK_OUTPUT_H_


#include "../drt_inpar/inpar_IO_runtime_vtk_output.H"

namespace Teuchos {
  class ParameterList;
}

namespace DRT {
namespace ELEMENTS {
  class StructureRuntimeVtuOutputParams;
  class BeamRuntimeVtuOutputParams;
}
}


namespace STR {
namespace TIMINT {

  /** \brief Input data container for VTK output at runtime for the structural (time) integration
   *
   * \author Maximilian Grill */
  class ParamsRuntimeVtkOutput
  {
  public:
    /// constructor
    ParamsRuntimeVtkOutput();

    /// destructor
    virtual ~ParamsRuntimeVtkOutput() {};

    /// initialize the class variables
    void Init( const Teuchos::ParameterList& IO_vtk_structure_paramslist );

    /// setup new class variables
    void Setup();


    /// whether to write numeric output data binary-encoded or in ascii format
    bool WriteBinaryOutput() const
    {
      CheckInitSetup();
      return ( output_data_format_ == INPAR::IO_RUNTIME_VTK::binary );
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
    Teuchos::RCP<const DRT::ELEMENTS::StructureRuntimeVtuOutputParams> GetStructureParams() const
    {
      CheckInitSetup();
      return params_runtime_vtu_output_structure_;
    };

    /// get the data container for parameters regarding beams
    Teuchos::RCP<const DRT::ELEMENTS::BeamRuntimeVtuOutputParams> GetBeamParams() const
    {
      CheckInitSetup();
      return params_runtime_vtu_output_beams_;
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
    bool isinit_;

    bool issetup_;
    /// @}

    /// @name variables controlling output
    /// @{

    /// data format for numeric output data
    INPAR::IO_RUNTIME_VTK::OutputDataFormat output_data_format_;

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int output_interval_steps_;

    /// whether to write output in every iteration of the nonlinear solver
    bool output_every_iteration_;

    /// whether to write output for structural elements
    bool output_structure_;

    /// whether to write special output for beam elements
    bool output_beams_;

    /// data container for input parameters related to VTU output of structure at runtime
    Teuchos::RCP<DRT::ELEMENTS::StructureRuntimeVtuOutputParams> params_runtime_vtu_output_structure_;

    /// data container for input parameters related to VTU output of beams at runtime
    Teuchos::RCP<DRT::ELEMENTS::BeamRuntimeVtuOutputParams> params_runtime_vtu_output_beams_;


    //@}
  };

}
}


#endif
