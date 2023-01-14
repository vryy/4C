/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTU output at runtime for fluid

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FLD_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_
#define FLD_DISCRETIZATION_RUNTIME_VTU_OUTPUT_PARAMS_H_

namespace Teuchos
{
  class ParameterList;
}



namespace DRT
{
  namespace ELEMENTS
  {
    /** \brief Input data container for VTU output at runtime for Fluid*/
    class FluidRuntimeVtuOutputParams
    {
     public:
      /// constructor
      FluidRuntimeVtuOutputParams() = default;

      /// destructor
      virtual ~FluidRuntimeVtuOutputParams(){};

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtk_fluid_paramslist);

      /// setup new class variables
      void Setup();

      /// whether to write velocity
      bool OutputVelocityState() const
      {
        CheckInitSetup();
        return output_velocity_state_;
      };

      /// whether to write pressure
      bool OutputPressureState() const
      {
        CheckInitSetup();
        return output_pressure_state_;
      };

      /// whether to write velocity
      bool OutputAccelerationState() const
      {
        CheckInitSetup();
        return output_acceleration_state_;
      };

      /// whether to write displacements
      bool OutputDisplacementState() const
      {
        CheckInitSetup();
        return output_displacement_state_;
      };

      /// whether to write velocity
      bool OutputGridVelocityState() const
      {
        CheckInitSetup();
        return output_gridvelocity_state_;
      };

      /// whether to write the owner of elements
      bool OutputElementOwner() const
      {
        CheckInitSetup();
        return output_element_owner_;
      };

      /// whether to write the GIDs of elements
      bool OutputElementGID() const
      {
        CheckInitSetup();
        return output_element_gid_;
      };

      /// whether to write the GIDs of the nodes
      bool OutputNodeGID() const
      {
        CheckInitSetup();
        return output_node_gid_;
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

      /// whether to write velocity output
      bool output_velocity_state_ = false;

      /// whether to write pressure output
      bool output_pressure_state_ = false;

      /// whether to write acceleration output
      bool output_acceleration_state_ = false;

      /// whether to write displacement output
      bool output_displacement_state_ = false;

      /// whether to write grid velocity output
      bool output_gridvelocity_state_ = false;

      /// whether to write the owner of elements
      bool output_element_owner_ = false;

      /// whether to write the element GIDs
      bool output_element_gid_ = false;

      /// whether to write the node GIDs
      bool output_node_gid_ = false;
      //@}
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
