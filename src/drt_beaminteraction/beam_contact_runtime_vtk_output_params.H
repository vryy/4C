/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_contact_runtime_vtk_output_params.H

\brief data container holding all input parameters for vtk-based visualization of beam contact

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMCONTACTVTKPARAMS_H
#define BEAMCONTACTVTKPARAMS_H

#include "../drt_inpar/inpar_beamcontact.H"


namespace BEAMINTERACTION
{
  class BeamContactRuntimeVtkOutputParams
  {
   public:
    //! constructor
    BeamContactRuntimeVtkOutputParams();

    //! destructor
    virtual ~BeamContactRuntimeVtkOutputParams(){};

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();


    /// whether to write numeric output data binary-encoded or in ascii format
    bool WriteBinaryOutput() const
    {
      ThrowErrorIfNotInitAndSetup();
      return (output_data_format_ == INPAR::BEAMCONTACT::binary);
    }

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int OutputIntervalInSteps() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_interval_steps_;
    };

    /// whether to write output in every iteration of the nonlinear solver
    bool OutputEveryIteration() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_every_iteration_;
    };

    /// whether to write output for contact forces
    bool IsWriteContactForces() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_forces_;
    };

    /// whether to write output for gaps
    bool IsWriteGaps() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_gaps_;
    };


   private:
    //! returns the isinit_ flag
    inline const bool& IsInit() const { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const { return issetup_; };

    //! asserts the init and setup status
    void ThrowErrorIfNotInitAndSetup() const;

    //! asserts the init status
    void ThrowErrorIfNotInit() const;


   private:
    bool isinit_;

    bool issetup_;

    /// data format for numeric output data
    INPAR::BEAMCONTACT::OutputDataFormat output_data_format_;

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int output_interval_steps_;

    /// whether to write output in every iteration of the nonlinear solver
    bool output_every_iteration_;

    /// whether to write forces
    bool output_forces_;

    /// whether to write gaps
    bool output_gaps_;
  };

}  // namespace BEAMINTERACTION

#endif
