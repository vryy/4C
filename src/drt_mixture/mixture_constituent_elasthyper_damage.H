/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a hyperelastic constituent with a damage process

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_ELASTHYPER_DAMAGE_H
#define BACI_MIXTURE_CONSTITUENT_ELASTHYPER_DAMAGE_H

#include "mixture_constituent_elasthyperbase.H"
#include "../drt_mat/anisotropy_extension.H"
#include "elastin_membrane_prestress_strategy.H"

namespace MAT
{
  class Anisotropy;
  namespace ELASTIC
  {
    class StructuralTensorStrategyBase;
    class IsoNeoHooke;
  }  // namespace ELASTIC
}  // namespace MAT

namespace MIXTURE
{
  class MixtureConstituent_ElastHyperDamage;

  namespace PAR
  {
    class MixtureConstituent_ElastHyperDamage
        : public MIXTURE::PAR::MixtureConstituent_ElastHyperBase
    {
     public:
      /*!
       * \brief Construct a new elastin material with a membrane
       *
       * \param matdata Material parameters
       * \param ref_mass_fraction reference mass fraction
       */
      explicit MixtureConstituent_ElastHyperDamage(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      /// default destructor
      ~MixtureConstituent_ElastHyperDamage() override = default;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::MixtureConstituent> CreateConstituent() override;

      /// @name material parameters
      /// @{
      const int damage_function_id_;
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Constituent for any hyperelastic material
   *
   * This constituent represents any hyperelastic material from the elasthyper toolbox. It has to
   * be paired with the MAT::Mixture_ElastHyper material and a MIXTURE::MixtureRule.
   */
  class MixtureConstituent_ElastHyperDamage : public MIXTURE::MixtureConstituent_ElastHyperBase
  {
   public:
    /*!
     * \brief Constructor for the material given the material parameters
     *
     * \param params Material parameters
     */
    explicit MixtureConstituent_ElastHyperDamage(
        MIXTURE::PAR::MixtureConstituent_ElastHyperDamage* params);

    /// Returns the material type enum
    INPAR::MAT::MaterialType MaterialType() const override;

    /*!
     * \brief Pack data into a char vector from this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param data (in/put) : vector storing all data to be packed into this instance.
     */
    void PackConstituent(DRT::PackBuffer& data) const override;

    /*!
     * \brief Unpack data from a char vector into this class to be called from a derived class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param position (in/out) : current position to unpack data
     * @param data (in) : vector storing all data to be unpacked into this instance.
     */
    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    /*!
     * Initialize the constituent with the parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) Parameter list for exchange of parameters
     */
    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;


    /*!
     * \brief Updates the material and all its summands
     *
     * This method is called once between each timestep after convergence.
     *
     * @param defgrd Deformation gradient
     * @param params Container for additional information
     * @param gp Gauß point
     * @param eleGID Global element identifier
     */
    void Update(LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    /*!
     * \brief Returns the current reference mass density of the constituent
     *
     * @param gp (in) : Gauss point
     * @return
     */
    double CurrentRefDensity(int gp) const override;

    /*!
     * \brief Standard evaluation of the material. This material does only support evaluation with
     * an elastic part.
     *
     * \param F Total deformation gradient
     * \param E_strain Green-Lagrange strain tensor
     * \param params Container for additional information
     * \param S_stress 2. Piola-Kirchhoff stress tensor in stress-like Voigt notation
     * \param cmat Constitutive tensor
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * \brief Evaluation of the constituent with an inelastic, external part.
     *
     * \param F Total deformation gradient
     * \param iF_in inverse inelastic (external) stretch tensor
     * \param params Container for additional information
     * \param S_stress 2. Piola Kirchhoff stress tensor in stress-like Voigt notation
     * \param cmat Constitutive tensor
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void EvaluateElasticPart(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iF_in,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

   private:
    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_ElastHyperDamage* params_;

    /// Current reference mass density
    std::vector<double> current_reference_density_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_ELASTHYPER_DAMAGE_H
