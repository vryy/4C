/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to output at runtime for beams

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_DISCRETIZATION_RUNTIME_OUTPUT_PARAMS_HPP
#define FOUR_C_STRUCTURE_NEW_DISCRETIZATION_RUNTIME_OUTPUT_PARAMS_HPP


#include "4C_config.hpp"

#include "4C_inpar_IO_runtime_output_structure_beams.hpp"
#include "4C_inpar_structure.hpp"

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /** \brief Input data container for output at runtime for structure*/
    class StructureRuntimeOutputParams
    {
     public:
      /// constructor
      StructureRuntimeOutputParams();

      /// destructor
      virtual ~StructureRuntimeOutputParams() = default;

      /// initialize the class variables
      void Init(const Teuchos::ParameterList& IO_vtk_structure_structure_paramslist);

      /// setup new class variables
      void Setup();

      /// whether to write displacements
      bool output_displacement_state() const
      {
        check_init_setup();
        return output_displacement_state_;
      };

      /// whether to write velocity
      bool output_velocity_state() const
      {
        check_init_setup();
        return output_velocity_state_;
      };

      /// whether to write the owner of elements
      bool output_element_owner() const
      {
        check_init_setup();
        return output_element_owner_;
      };

      /// whether to write the GIDs of elements
      bool output_element_gid() const
      {
        check_init_setup();
        return output_element_gid_;
      };

      /// whether to write the ghosting information of elements
      bool output_element_ghosting() const
      {
        check_init_setup();
        return output_element_ghosting_;
      };

      /// whether to write the GIDs of the nodes
      bool output_node_gid() const
      {
        check_init_setup();
        return output_node_gid_;
      };

      /// whether to write stress and / or strain data
      bool output_stress_strain() const
      {
        check_init_setup();
        return output_stress_strain_;
      };

      /// Return output type of Gauss point data
      Inpar::STR::GaussPointDataOutputType gauss_point_data_output() const
      {
        check_init_setup();
        return gauss_point_data_output_type_;
      }

     private:
      /// get the init indicator status
      const bool& is_init() const { return isinit_; };

      /// get the setup indicator status
      const bool& is_setup() const { return issetup_; };

      /// Check if Init() and Setup() have been called, yet.
      void check_init_setup() const;


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

      /// whether to write velocity output
      bool output_velocity_state_;

      /// whether to write the owner of elements
      bool output_element_owner_;

      /// whether to write the element GIDs
      bool output_element_gid_;

      /// whether to write the element ghosting information
      bool output_element_ghosting_;

      /// whether to write the node GIDs
      bool output_node_gid_;

      /// whether to write stress and / or strain data
      bool output_stress_strain_;

      /// Output type of Gauss point data
      Inpar::STR::GaussPointDataOutputType gauss_point_data_output_type_;
      //@}
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
