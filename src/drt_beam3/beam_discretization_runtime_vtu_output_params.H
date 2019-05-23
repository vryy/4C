/*-----------------------------------------------------------------------------------------------*/
/*!

\brief input parameters related to VTU output at runtime for beams

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAM_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_
#define BEAM_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_


#include "../drt_inpar/inpar_IO_runtime_vtk_output_structure_beams.H"

namespace Teuchos
{
  class ParameterList;
}


namespace DRT
{
  namespace ELEMENTS
  {
    /** \brief Input data container for VTU output at runtime for beams
     *
     * \author Maximilian Grill */
    class BeamRuntimeVtuOutputParams
    {
     public:
      /// constructor
      BeamRuntimeVtuOutputParams();

      /// destructor
      virtual ~BeamRuntimeVtuOutputParams(){};

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtk_structure_beams_paramslist);

      /// setup new class variables
      void Setup();


      /// whether to write displacements
      bool OutputDisplacementState() const
      {
        CheckInitSetup();
        return output_displacement_state_;
      };

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
      bool IsWriteInternalEnergyElement() const
      {
        CheckInitSetup();
        return write_internal_energy_element_;
      };

      /// whether to write material cross-section strains at the Gauss points
      bool IsWriteKineticEnergyElement() const
      {
        CheckInitSetup();
        return write_kinetic_energy_element_;
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

      /// whether to write material cross-section stresses at the Gauss points
      bool IsWriteSpatialStressesGaussPoints() const
      {
        CheckInitSetup();
        return write_spatial_crosssection_stresses_gausspoints_;
      };

      /// whether to write material cross-section stresses at the Gauss points
      bool IsWriteElementFilamentCondition() const
      {
        CheckInitSetup();
        return write_filament_condition_;
      };

      /// whether to write element and network orientation parameter
      bool IsWriteOrientationParamter() const
      {
        CheckInitSetup();
        return write_orientation_parameter_;
      };

      /// whether to write crosssection forces of periodic rve in x, y, and z direction
      bool IsWriteRVECrosssectionForces() const
      {
        CheckInitSetup();
        return write_rve_crosssection_forces_;
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

      /// whether to write displacement output
      bool output_displacement_state_;

      /// whether to use abolute positions or initial positions for the vtu geometry definition
      /// (i.e. for the visualization point coordinates)
      bool use_absolute_positions_visualizationpoint_coordinates_;

      /// whether to write internal (elastic) energy for each element
      bool write_internal_energy_element_;

      /// whether to write kinetic energy for each element
      bool write_kinetic_energy_element_;

      /// whether to write triads at the visualization points
      bool write_triads_visualizationpoints_;

      /// whether to write material cross-section strains at the Gauss points
      bool write_material_crosssection_strains_gausspoints_;

      /// whether to write material cross-section stresses at the Gauss points
      bool write_material_crosssection_stresses_gausspoints_;

      /// whether to write spatial cross-section stresses at the Gauss points
      bool write_spatial_crosssection_stresses_gausspoints_;

      /// whether to write beam filament condition (id, type)
      bool write_filament_condition_;

      /// whether to write element and network orientation parameter
      bool write_orientation_parameter_;

      /// whether to write crosssection forces of periodic rve in x, y, and z direction
      bool write_rve_crosssection_forces_;

      //@}
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
