/*----------------------------------------------------------------------*/
/*! \file
\brief A condition of any kind

\level 1

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_MATERIAL_INPUT_BASE_HPP
#define FOUR_C_MATERIAL_INPUT_BASE_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_material_parameter_base.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace CORE::MAT::PAR
{
  // forward declaration
  class Parameter;

  /// Container for read-in materials
  ///
  /// This object stores the validated material parameters as
  /// IO::InputParameterContainer.
  class Material : public IO::InputParameterContainer
  {
   public:
    /// @name life span
    //@{

    /// Standard constructor
    Material(const int id,                         ///< unique material ID
        const CORE::Materials::MaterialType type,  ///< type of material
        const std::string name                     ///< name of material
    );

    /// default constructor
    Material() = default;

    /// Copy constructor
    ///
    /// Makes a deep copy of the material parameters
    Material(const CORE::MAT::PAR::Material& old);

    /// Set pointer to readily allocated 'quick access' material parameters
    ///
    /// This function is called by the material factory MAT::Factory.
    /// To circumvent more than this single major switch of material type to
    /// object, #params_ are allocated externally.
    inline void SetParameter(CORE::MAT::PAR::Parameter* matparam)
    {
      params_ = Teuchos::rcp(matparam);
    }

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
    [[nodiscard]] inline virtual CORE::Materials::MaterialType Type() const { return type_; }

    /// Return quick accessible material parameter data
    ///
    /// These quick access parameters are stored in separate member #params_;
    /// whereas the originally read ones are stored in IO::InputParameterContainer base
    [[nodiscard]] inline CORE::MAT::PAR::Parameter* Parameter() const { return params_.get(); }

    //@}

   protected:
    /// don't want = operator
    Material operator=(const Material& old);

    /// Unique ID of this material, no second material of same ID may exist
    int id_{};

    /// Type of this material
    CORE::Materials::MaterialType type_{};

    /// Name
    std::string name_{};

    /// Unwrapped material data for 'quick' access
    Teuchos::RCP<CORE::MAT::PAR::Parameter> params_{};
  };
}  // namespace CORE::MAT::PAR



/// out stream operator
std::ostream& operator<<(std::ostream& os, const CORE::MAT::PAR::Material& cond);


FOUR_C_NAMESPACE_CLOSE

#endif
