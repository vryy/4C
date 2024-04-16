/*----------------------------------------------------------------------*/
/*! \file
\brief A condition of any kind

\level 1

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_PAR_MATERIAL_HPP
#define FOUR_C_MAT_PAR_MATERIAL_HPP

#include "baci_config.hpp"

#include "baci_inpar_container.hpp"
#include "baci_inpar_material.hpp"
#include "baci_mat_par_parameter.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


namespace MAT::PAR
{
  // forward declaration
  class Parameter;

  /// Container for read-in materials
  ///
  /// This object stores the validated material parameters as
  /// INPAR::InputParameterContainer.
  class Material : public INPAR::InputParameterContainer
  {
   public:
    /// @name life span
    //@{

    /// Standard constructor
    Material(const int id,                    ///< unique material ID
        const INPAR::MAT::MaterialType type,  ///< type of material
        const std::string name                ///< name of material
    );

    /// default constructor
    Material() = default;

    /// Copy constructor
    ///
    /// Makes a deep copy of the material parameters
    Material(const MAT::PAR::Material& old);

    /// Set pointer to readily allocated 'quick access' material parameters
    ///
    /// This function is called by the material factory MAT::Material::Factory.
    /// To circumvent more than this single major switch of material type to
    /// object, #params_ are allocated externally.
    inline void SetParameter(MAT::PAR::Parameter* matparam) { params_ = Teuchos::rcp(matparam); }

    //@}

    /// @name Query methods
    //@{

    /// Return material id
    [[nodiscard]] inline virtual int Id() const { return id_; }

    /// Return material name
    [[nodiscard]] inline virtual std::string Name() const { return name_; }

    /// Print this Condition (std::ostream << is also implemented for DRT::Condition)
    void Print(std::ostream& os) const override;

    /// Return type of condition
    [[nodiscard]] inline virtual INPAR::MAT::MaterialType Type() const { return type_; }

    /// Return quick accessible material parameter data
    ///
    /// These quick access parameters are stored in separate member #params_;
    /// whereas the originally read ones are stored in INPAR::InputParameterContainer base
    [[nodiscard]] inline MAT::PAR::Parameter* Parameter() const { return params_.get(); }

    //@}

   protected:
    /// don't want = operator
    Material operator=(const Material& old);

    /// Unique ID of this material, no second material of same ID may exist
    int id_{};

    /// Type of this material
    INPAR::MAT::MaterialType type_{};

    /// Name
    std::string name_{};

    /// Unwrapped material data for 'quick' access
    Teuchos::RCP<MAT::PAR::Parameter> params_{};
  };
}  // namespace MAT::PAR



/// out stream operator
std::ostream& operator<<(std::ostream& os, const MAT::PAR::Material& cond);


BACI_NAMESPACE_CLOSE

#endif
