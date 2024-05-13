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

#include "4C_legacy_enum_definitions_materials.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace CORE::MAT
{
  class Material;
  namespace PAR
  {
    class Material;
  }
}  // namespace CORE::MAT

/*----------------------------------------------------------------------*/
/* declarations */

namespace CORE::MAT::PAR
{
  /*----------------------------------------------------------------------*/
  /// Base object to hold 'quick' access material parameters
  ///
  /// CORE::MAT::PAR::Parameters is derived for the various implemented
  /// materials. These provide the 'quick' access to the read-in
  /// material parameters.
  ///
  /// For every read-in material will exist a single instance (of
  /// a derived class) of this object.
  class Parameter
  {
   public:
    /// construct the material object given material parameters
    Parameter(Teuchos::RCP<const CORE::MAT::PAR::Material>
            matdata  ///< read and validated material data (of 'slow' access)
    );

    /// destructor
    virtual ~Parameter() = default;

    /// (unique) material ID
    [[nodiscard]] int Id() const { return id_; }

    /// material type
    [[nodiscard]] CORE::Materials::MaterialType Type() const { return type_; }

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<CORE::MAT::Material> CreateMaterial() = 0;

    //! \brief set element specific or global material parameter using enum parametername which is
    //! defined in respective MAT::PAR classes
    void SetParameter(int parametername, Teuchos::RCP<Epetra_Vector> myparameter);

    //! \brief set element specific or global material parameter using enum parametername which is
    //! defined in respective MAT::PAR classes
    void SetParameter(int parametername, const double val, const int eleGID);

    //! \brief return element specific or global material parameter using enum parametername which
    //! is defined in respective MAT::PAR classes
    double GetParameter(int parametername, const int EleId);

    /// return matparams_
    virtual std::vector<Teuchos::RCP<Epetra_Vector>>& ReturnMatparams() { return matparams_; }

   protected:
    /*! \brief
     * data structure to store all material parameters in.
     * By default all elements with the same mat share the same material properties, hence the
     * Epetra_Vector has length 1 However for elementwise material properties the Epetra_Vector
     * has EleColMap layout.
     */
    std::vector<Teuchos::RCP<Epetra_Vector>> matparams_;

   private:
    /// material ID, as defined in input file
    int id_;

    /// material type
    CORE::Materials::MaterialType type_;

  };  // class Parameter

}  // namespace CORE::MAT::PAR


FOUR_C_NAMESPACE_CLOSE

#endif
