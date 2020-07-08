/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container for input parameters for vtk-based visualization of potential-based beam
       interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMPOTENTIALVTKPARAMS_H
#define BEAMPOTENTIALVTKPARAMS_H

#include "../drt_inpar/inpar_beampotential.H"


namespace BEAMINTERACTION
{
  /*!
   *  */
  class BeamToBeamPotentialRuntimeVtkParams
  {
   public:
    //! constructor
    BeamToBeamPotentialRuntimeVtkParams();

    //! destructor
    virtual ~BeamToBeamPotentialRuntimeVtkParams(){};

    //! initialize with the stuff coming from input file
    void Init(const Teuchos::ParameterList& beam_contact_vtk_paramslist);

    //! setup member variables
    void Setup();


    /// whether to write numeric output data binary-encoded or in ascii format
    bool WriteBinaryOutput() const
    {
      ThrowErrorIfNotInitAndSetup();
      return (output_data_format_ == INPAR::BEAMPOTENTIAL::binary);
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

    /// whether to write output for forces
    bool IsWriteForces() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_forces_;
    };

    /// whether to write output for moments
    bool IsWriteMoments() const
    {
      ThrowErrorIfNotInitAndSetup();
      return output_moments_;
    };

    /// whether to write forces/moments separately for each element pair
    bool IsWriteForcesMomentsPerElementPair() const
    {
      ThrowErrorIfNotInitAndSetup();
      return write_force_moment_per_elepair_;
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
    INPAR::BEAMPOTENTIAL::OutputDataFormat output_data_format_;

    /// output interval regarding steps: write output every INTERVAL_STEPS steps
    int output_interval_steps_;

    /// whether to write output in every iteration of the nonlinear solver
    bool output_every_iteration_;

    /// whether to write forces
    bool output_forces_;

    /// whether to write moments
    bool output_moments_;

    /// whether to write forces/moments separately for each element pair
    bool write_force_moment_per_elepair_;
  };

}  // namespace BEAMINTERACTION

#endif
