/*----------------------------------------------------------------------*/
/*! \file
\brief Base object to hold 'quick' access to material parameters

\level 1

*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef BACI_MAT_PAR_PARAMETER_HPP
#define BACI_MAT_PAR_PARAMETER_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include "baci_inpar_material.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace MAT
{
  class Material;
  namespace PAR
  {
    class Material;
  }

  namespace ELASTIC
  {
    class StructuralTensorStrategyBase;
  }
}  // namespace MAT

/*----------------------------------------------------------------------*/
/* declarations */

namespace MAT::PAR
{
  /*----------------------------------------------------------------------*/
  /// Base object to hold 'quick' access material parameters
  ///
  /// MAT::PAR::Parameters is derived for the various implemented
  /// materials. These provide the 'quick' access to the read-in
  /// material parameters.
  ///
  /// For every read-in material will exist a single instance (of
  /// a derived class) of this object.
  class Parameter
  {
   public:
    /// construct the material object given material parameters
    Parameter(Teuchos::RCP<const MAT::PAR::Material>
            matdata  ///< read and validated material data (of 'slow' access)
    );

    /// destructor
    virtual ~Parameter() = default;

    /// (unique) material ID
    [[nodiscard]] int Id() const { return id_; }

    /// material type
    [[nodiscard]] INPAR::MAT::MaterialType Type() const { return type_; }

    /// material name
    [[nodiscard]] std::string Name() const { return name_; }

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial() = 0;

    //! \brief set element specific or global material parameter using enum parametername which is
    //! defined in respective MAT::PAR classes
    void SetParameter(int parametername, Teuchos::RCP<Epetra_Vector> myparameter);

    //! \brief set element specific or global material parameter using enum parametername which is
    //! defined in respective MAT::PAR classes
    void SetParameter(int parametername, const double val, const int eleGID);

    //! brief extend all RCP<Epetra_Vectors> which have length one to element colmap layout
    void ExpandParametersToEleColLayout();

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
    INPAR::MAT::MaterialType type_;

    /// material name
    std::string name_;

  };  // class Parameter


  /// Extension to hold 'quick' access material parameters for anisotropy
  class ParameterAniso : public Parameter
  {
   public:
    /// construct the material object given material parameters
    ParameterAniso(Teuchos::RCP<const MAT::PAR::Material> matdata);

    /// return pointer to strategy
    const Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase>& StructuralTensorStrategy()
    {
      return structural_tensor_strategy_;
    };

   private:
    /// structural tensor strategy
    Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase> structural_tensor_strategy_;

  };  // class ParameterAniso

}  // namespace MAT::PAR


BACI_NAMESPACE_CLOSE

#endif  // MAT_PAR_PARAMETER_H
