/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to solid surface output (visualization) parameters.

\level 3

\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_VTK_OUTPUT_PARAMS_H
#define BEAM_TO_SOLID_SURFACE_VTK_OUTPUT_PARAMS_H


#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_IO_runtime_vtk_output.H"


namespace BEAMINTERACTION
{
  class BeamToSolidSurfaceVtkOutputParams
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToSolidSurfaceVtkOutputParams();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceVtkOutputParams(){};

    /**
     * \brief Initialize with the stuff coming from input file.
     */
    void Init();

    /**
     * \brief Setup member variables.
     */
    void Setup();

    /**
     * \brief Return the output every iteration flag.
     */
    bool GetOutputEveryIteration() const { return output_every_iteration_; }

    /**
     * \brief If the output should be binary or ascii
     */
    bool WriteBinaryOutput() const
    {
      CheckInitSetup();
      return (output_data_format_ == INPAR::IO_RUNTIME_VTK::binary);
    }

    /**
     * \brief Output interval regarding steps: write output every INTERVAL_STEPS steps.
     */
    int OutputIntervalInSteps() const
    {
      CheckInitSetup();
      return output_interval_steps_;
    };

    /**
     * \brief Return the output flag.
     */
    bool GetOutputFlag() const { return output_flag_; }

    /**
     * \brief Return the nodal forces flag.
     */
    bool GetNodalForceOutputFlag() const { return nodal_forces_; }

    /**
     * \brief Return the averaged nodal normals flag.
     */
    bool GetAveragedNormalsOutputFlag() const { return nodal_forces_; }

    /**
     * \brief Return the mortar lambda discrete flag.
     */
    bool GetMortarLambdaDiscretOutputFlag() const { return mortar_lambda_discret_; }

    /**
     * \brief Return the mortar lambda continuous flag.
     */
    bool GetMortarLambdaContinuousOutputFlag() const { return mortar_lambda_continuous_; }

    /**
     * \brief Return the number of segments for continuous mortar output.
     */
    unsigned int GetMortarLambdaContinuousSegments() const
    {
      return mortar_lambda_continuous_segments_;
    }

    /**
     * \brief Return the nodal forces flag.
     */
    bool GetSegmentationOutputFlag() const { return segmentation_; }

    /**
     * \brief Return the integration output flag.
     */
    bool GetIntegrationPointsOutputFlag() const { return integration_points_; }

   protected:
    /**
     * \brief Checks the init and setup status.
     */
    inline void CheckInitSetup() const
    {
      if (!isinit_ or !issetup_) dserror("Call Init() and Setup() first!");
    }

    /**
     * \brief Checks the init status.
     */
    inline void CheckInit() const
    {
      if (!isinit_) dserror("Init() has not been called, yet!");
    }

    //! Flag if object is initialized.
    bool isinit_;

    //! Flag if object is set up.
    bool issetup_;

    //! Data format that will be used for the vtk files.
    INPAR::IO_RUNTIME_VTK::OutputDataFormat output_data_format_;

    //! Output interval regarding steps: write output every INTERVAL_STEPS steps.
    int output_interval_steps_;

    //! Whether to write output in every iteration of the nonlinear solver.
    bool output_every_iteration_;

    //! Flag whether or not to write output.
    bool output_flag_;

    //! Flag whether or not to output resulting nodal forces.
    bool nodal_forces_;

    //! Flag whether or not to output averaged nodal normals.
    bool averaged_normals_;

    //! Flag whether or not to output resulting discrete Lagrange multiplier values for mortar
    //! pairs.
    bool mortar_lambda_discret_;

    //! Flag whether or not to output resulting continuous Lagrange multiplier function for mortar
    //! pairs.
    bool mortar_lambda_continuous_;

    //! Number of segments to use for the continuous mortar output.
    unsigned int mortar_lambda_continuous_segments_;

    //! Flag whether or not to write segmentation data.
    bool segmentation_;

    //! Flag whether or not to write data at integration points.
    bool integration_points_;
  };

}  // namespace BEAMINTERACTION

#endif
