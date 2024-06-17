/*----------------------------------------------------------------------*/
/*! \file
\brief Base class to hold material parameters
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATERIAL_PARAMETER_BASE_HPP
#define FOUR_C_MATERIAL_PARAMETER_BASE_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::Mat
{
  class Material;
}  // namespace Core::Mat

namespace Core::Mat::PAR
{
  /**
   * Base class for material parameters.
   *
   * This class is used to store the parameters of a material once. Through the virtual method
   * create_material() a material instance can be created which will use the stored parameters.
   * Thus, the parameters in this class are shared across all created material instances.
   *
   * When defining a new material, a new derived class of Parameter must be created. This derived
   * class must implement the create_material() method to create the material instance of the
   * correct type.
   */
  class Parameter
  {
   public:
    /**
     * Additional data that is required to create the enclosing Parameter object.
     */
    struct Data
    {
      /**
       * An ID that is unique for each material and supplied from the input file.
       */
      int id{-1};

      /**
       * The type of the material.
       */
      Core::Materials::MaterialType type{Core::Materials::MaterialType::m_none};

      /**
       * A dynamic container for any additional parameters that a concrete derived class may
       * require.
       */
      Core::IO::InputParameterContainer parameters;
    };

    /**
     * Base class constructor taking the @p data with basic information about the material.
     */
    Parameter(Data data);

    /**
     * Virtual destructor.
     */
    virtual ~Parameter() = default;

    /**
     * The ID of the material.
     */
    [[nodiscard]] int Id() const { return data_.id; }

    /**
     * The type of the material.
     */
    [[nodiscard]] Core::Materials::MaterialType Type() const { return data_.type; }


    /**
     * Create a new material instance from the stored parameters.
     *
     * @note This method must be implemented by derived classes.
     */
    virtual Teuchos::RCP<Core::Mat::Material> create_material() = 0;

    /**
     * Access to the raw input parameters. Derived classes usually store the parameters in a more
     * convenient way. This method grants access to the raw input parameters without specific
     * knowledge of the derived class.
     */
    [[nodiscard]] const Core::IO::InputParameterContainer& raw_parameters() const
    {
      return data_.parameters;
    }


   private:
    /**
     * Data as supplied during construction.
     */
    Data data_;
  };
}  // namespace Core::Mat::PAR

FOUR_C_NAMESPACE_CLOSE

#endif
