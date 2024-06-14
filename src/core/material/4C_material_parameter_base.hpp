/*----------------------------------------------------------------------*/
/*! \file
\brief Base object to hold 'quick' access to material parameters

\level 1

*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATERIAL_PARAMETER_BASE_HPP
#define FOUR_C_MATERIAL_PARAMETER_BASE_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace Core::Mat
{
  class Material;
  namespace PAR
  {
    class Parameter;
  }
}  // namespace Core::Mat

/*----------------------------------------------------------------------*/
/* declarations */

namespace Core::Mat::PAR
{
  /*----------------------------------------------------------------------*/
  /// Base class to hold material parameters
  ///
  /// Core::Mat::PAR::Parameters is derived for the various implemented
  /// materials. These provide the 'quick' access to the read-in
  /// material parameters.
  ///
  /// For every read-in material will exist a single instance (of
  /// a derived class) of this object.
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

    /// destructor
    virtual ~Parameter() = default;

    /// (unique) material ID
    [[nodiscard]] int Id() const { return data_.id; }

    /// material type
    [[nodiscard]] Core::Materials::MaterialType Type() const { return data_.type; }

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<Core::Mat::Material> create_material() = 0;

    //! \brief return element specific or global material parameter using enum parametername which
    //! is defined in respective Mat::PAR classes
    double GetParameter(int parametername, const int EleId);

    /**
     * Access to the raw input data.
     */
    const Core::IO::InputParameterContainer& raw_parameters() const { return data_.parameters; }

   protected:
    /*! \brief
     * data structure to store all material parameters in.
     * By default all elements with the same mat share the same material properties, hence the
     * Epetra_Vector has length 1 However for elementwise material properties the Epetra_Vector
     * has EleColMap layout.
     */
    std::vector<Teuchos::RCP<Epetra_Vector>> matparams_;

   private:
    /**
     * Data as supplied during construction.
     */
    Data data_;
  };  // class Parameter

}  // namespace Core::Mat::PAR


FOUR_C_NAMESPACE_CLOSE

#endif
