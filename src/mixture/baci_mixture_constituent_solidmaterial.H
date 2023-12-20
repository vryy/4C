/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a general solid material constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_SOLIDMATERIAL_H
#define BACI_MIXTURE_CONSTITUENT_SOLIDMATERIAL_H

#include "baci_config.H"

#include "baci_mat_par_parameter.H"
#include "baci_mat_so3_material.H"
#include "baci_mixture_constituent.H"

BACI_NAMESPACE_OPEN


namespace MIXTURE
{
  class MixtureConstituent_SolidMaterial;

  namespace PAR
  {
    class MixtureConstituent_SolidMaterial : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituent_SolidMaterial(const Teuchos::RCP<MAT::PAR::Material>& matdata);
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
  class MixtureConstituent_SolidMaterial : public MIXTURE::MixtureConstituent
  {
   public:
    /// Constructor for the material given the material parameters
    explicit MixtureConstituent_SolidMaterial(
        MIXTURE::PAR::MixtureConstituent_SolidMaterial* params, int id);

    void PackConstituent(CORE::COMM::PackBuffer& data) const override;

    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    INPAR::MAT::MaterialType MaterialType() const override;

    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

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
    MIXTURE::PAR::MixtureConstituent_SolidMaterial* params_;

    // reference to the so3 material
    Teuchos::RCP<MAT::So3Material> material_;
  };

}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif  // MIXTURE_CONSTITUENT_SOLIDMATERIAL_H