/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a general solid material constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_SOLIDMATERIAL_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_SOLIDMATERIAL_HPP

#include "4C_config.hpp"

#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"

FOUR_C_NAMESPACE_OPEN


namespace MIXTURE
{
  class MixtureConstituentSolidMaterial;

  namespace PAR
  {
    class MixtureConstituentSolidMaterial : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituentSolidMaterial(const Core::Mat::PAR::Parameter::Data& matdata);
      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> create_constituent(int id) override;

      /// @name material parameters
      /// @{
      /// Id of the solid material
      const int matid_;
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Constituent for any solid material
   *
   * This constituent represents any solid material from the material toolbox. It has to
   * be paired with the Mat::Mixture and a MIXTURE::MixtureRule.
   */
  class MixtureConstituentSolidMaterial : public MIXTURE::MixtureConstituent
  {
   public:
    /// Constructor for the material given the material parameters
    explicit MixtureConstituentSolidMaterial(
        MIXTURE::PAR::MixtureConstituentSolidMaterial* params, int id);

    void pack_constituent(Core::Communication::PackBuffer& data) const override;

    void unpack_constituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    Core::Materials::MaterialType material_type() const override;

    void read_element(int numgp, Input::LineDefinition* linedef) override;

    void update(Core::LinAlg::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params,
        const int gp, const int eleGID) override;

    void evaluate(const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
        Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool evaluate_output_data(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

   private:
    /// my material parameters
    MIXTURE::PAR::MixtureConstituentSolidMaterial* params_;

    // reference to the so3 material
    Teuchos::RCP<Mat::So3Material> material_;
  };

}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif