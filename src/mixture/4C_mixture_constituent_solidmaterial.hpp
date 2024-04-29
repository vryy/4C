/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a general solid material constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_SOLIDMATERIAL_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_SOLIDMATERIAL_HPP

#include "4C_config.hpp"

#include "4C_mat_par_parameter.hpp"
#include "4C_mat_so3_material.hpp"
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
      explicit MixtureConstituentSolidMaterial(const Teuchos::RCP<MAT::PAR::Material>& matdata);
      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> CreateConstituent(int id) override;

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
   * be paired with the MAT::Mixture and a MIXTURE::MixtureRule.
   */
  class MixtureConstituentSolidMaterial : public MIXTURE::MixtureConstituent
  {
   public:
    /// Constructor for the material given the material parameters
    explicit MixtureConstituentSolidMaterial(
        MIXTURE::PAR::MixtureConstituentSolidMaterial* params, int id);

    void PackConstituent(CORE::COMM::PackBuffer& data) const override;

    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    CORE::Materials::MaterialType MaterialType() const override;

    void ReadElement(int numgp, INPUT::LineDefinition* linedef) override;

    void Update(CORE::LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params,
        const int gp, const int eleGID) override;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    void RegisterOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateOutputData(
        const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const override;

   private:
    /// my material parameters
    MIXTURE::PAR::MixtureConstituentSolidMaterial* params_;

    // reference to the so3 material
    Teuchos::RCP<MAT::So3Material> material_;
  };

}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif