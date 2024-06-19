/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to output at runtime for fluid

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_FLUID_DISCRETIZATION_RUNTIME_OUTPUT_PARAMS_HPP
#define FOUR_C_FLUID_DISCRETIZATION_RUNTIME_OUTPUT_PARAMS_HPP

#include "4C_config.hpp"

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /** \brief Input data container for output at runtime for Fluid*/
    class FluidRuntimeOutputParams
    {
     public:
      /// constructor
      FluidRuntimeOutputParams() = default;

      /// destructor
      virtual ~FluidRuntimeOutputParams() = default;

      /// initialize the class variables
      void init(const Teuchos::ParameterList& IO_fluid_paramslist);

      /// setup new class variables
      void setup();

      /// whether to write velocity
      bool OutputVelocityState() const
      {
        check_init_setup();
        return output_velocity_state_;
      };

      /// whether to write pressure
      bool OutputPressureState() const
      {
        check_init_setup();
        return output_pressure_state_;
      };

      /// whether to write velocity
      bool output_acceleration_state() const
      {
        check_init_setup();
        return output_acceleration_state_;
      };

      /// whether to write displacements
      bool output_displacement_state() const
      {
        check_init_setup();
        return output_displacement_state_;
      };

      /// whether to write velocity
      bool output_grid_velocity_state() const
      {
        check_init_setup();
        return output_gridvelocity_state_;
      };

      /// whether to write the owner of elements
      bool OutputElementOwner() const
      {
        check_init_setup();
        return output_element_owner_;
      };

      /// whether to write the GIDs of elements
      bool OutputElementGID() const
      {
        check_init_setup();
        return output_element_gid_;
      };

      /// whether to write the GIDs of the nodes
      bool OutputNodeGID() const
      {
        check_init_setup();
        return output_node_gid_;
      };

     private:
      /// get the init indicator status
      const bool& is_init() const { return isinit_; };

      /// get the setup indicator status
      const bool& is_setup() const { return issetup_; };

      /// Check if init() and setup() have been called, yet.
      void check_init_setup() const;


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
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
