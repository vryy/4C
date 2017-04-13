/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_discretization_runtime_vtu_output_params.H

\brief input parameters related to VTU output at runtime for beams

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAM_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_
#define BEAM_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_


#include "../drt_inpar/inpar_IO_runtime_vtk_output_structure_beams.H"

namespace Teuchos {
  class ParameterList;
}


namespace DRT {
namespace ELEMENTS {

  /** \brief Input data container for VTU output at runtime for beams
   *
   * \author Maximilian Grill */
  class BeamRuntimeVtuOutputParams
  {
  public:
    /// constructor
    BeamRuntimeVtuOutputParams();

    /// destructor
    virtual ~BeamRuntimeVtuOutputParams() {};

    /// initialize the class variables
    void Init( const Teuchos::ParameterList& IO_vtk_structure_beams_paramslist );

    /// setup new class variables
    void Setup();


    /// whether to write triads at the Gauss points
    bool IsWriteTriadVisualizationPoints() const
    {
      CheckInitSetup();
      return write_triads_visualizationpoints_;
    };

    /// whether to use abolute positions or initial positions for the vtu geometry definition
    /// (i.e. for the visualization point coordinates)
    bool UseAbsolutePositions() const
    {
      CheckInitSetup();
      return use_absolute_positions_visualizationpoint_coordinates_;
    };

    /// whether to write material cross-section strains at the Gauss points
    bool IsWriteMaterialStrainsGaussPoints() const
    {
      CheckInitSetup();
      return write_material_crosssection_strains_gausspoints_;
    };

    /// whether to write material cross-section stresses at the Gauss points
    bool IsWriteMaterialStressesGaussPoints() const
    {
      CheckInitSetup();
      return write_material_crosssection_stresses_gausspoints_;
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

    /// whether to use abolute positions or initial positions for the vtu geometry definition
    /// (i.e. for the visualization point coordinates)
    bool use_absolute_positions_visualizationpoint_coordinates_;

    /// whether to write triads at the visualization points
    bool write_triads_visualizationpoints_;

    /// whether to write material cross-section strains at the Gauss points
    bool write_material_crosssection_strains_gausspoints_;

    /// whether to write material cross-section stresses at the Gauss points
    bool write_material_crosssection_stresses_gausspoints_;

    //@}
  };

}
}

#endif
