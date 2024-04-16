/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a hyperelastic constituent with a damage process and a 2D membrane material

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_ELASTHYPER_ELASTIN_MEMBRANE_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_ELASTHYPER_ELASTIN_MEMBRANE_HPP

#include "baci_config.hpp"

#include "baci_mat_anisotropy_extension.hpp"
#include "baci_mixture_constituent_elasthyperbase.hpp"
#include "baci_mixture_elastin_membrane_prestress_strategy.hpp"

BACI_NAMESPACE_OPEN

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
  class MixtureConstituent_ElastHyperElastinMembrane;

  /*!
   * \brief Anisotropy extension for elastin material.
   *
   * The anisotropy extension provides the structural tensor of the plane of the membrane, which is
   * orthogonal to the radial direction
   */
  class ElastinMembraneAnisotropyExtension : public MAT::FiberAnisotropyExtension<1>
  {
   public:
    /*!
     * \brief Constructor of the new elastin anisotropy extension
     *
     * \param structuralTensorStrategy Structural tensor strategy to compute the structural tensors
     */
    explicit ElastinMembraneAnisotropyExtension(
        const Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy);

    /*!
     * \brief This method will be called when all global data is initialized. Here we need to create
     * the membrane plane structural tensor
     */
    void OnGlobalDataInitialized() override;

    /*!
     * \brief Returns the structural tensor of the membrane plane at the Gauss point
     *
     * \param gp (in) : Gauss point
     * \return const CORE::LINALG::Matrix<3, 3>& Reference to the structural tensor of the membrane
     * plane
     */
    const CORE::LINALG::Matrix<3, 3>& GetOrthogonalStructuralTensor(int gp);

   private:
    /// Holder of the internal structural tensors
    std::vector<CORE::LINALG::Matrix<3, 3>> orthogonalStructuralTensor_;
  };

  namespace PAR
  {
    class MixtureConstituent_ElastHyperElastinMembrane
        : public MIXTURE::PAR::MixtureConstituent_ElastHyperBase
    {
     public:
      /*!
       * \brief Construct a new elastin material with a membrane
       *
       * \param matdata Material parameters
       * \param ref_mass_fraction reference mass fraction
       */
      explicit MixtureConstituent_ElastHyperElastinMembrane(
          const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> CreateConstituent(int id) override;

      /// @name material parameters
      /// @{
      const int damage_function_id_;

      /// number of summands
      const int nummat_membrane_;

      /// List of material ids of the summands
      const std::vector<int>* matids_membrane_;
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Constituent for any hyperelastic material
   *
   * This constituent represents any hyperelastic material from the elasthyper toolbox. It has to
   * be paired with the MAT::Mixture material and a MIXTURE::MixtureRule.
   */
  class MixtureConstituent_ElastHyperElastinMembrane
      : public MIXTURE::MixtureConstituent_ElastHyperBase,
        public ElastinMembraneEvaluation
  {
   public:
    /*!
     * \brief Constructor for the material given the material parameters
     *
     * \param params Material parameters
     */
    explicit MixtureConstituent_ElastHyperElastinMembrane(
        MIXTURE::PAR::MixtureConstituent_ElastHyperElastinMembrane* params, int id);

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
    void PackConstituent(CORE::COMM::PackBuffer& data) const override;

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
     * \brief Register anisotropy extensions to the global anisotropy manager
     *
     * \param anisotropy Reference to the global anisotropy manager
     */
    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

    /*!
     * Initialize the constituent with the parameters of the input line
     *
     * @param numgp (in) Number of Gauss-points
     * @param params (in/out) Parameter list for exchange of parameters
     */
    void ReadElement(int numgp, INPUT::LineDefinition* linedef) override;


    /*!
     * \brief Updates the material and all its summands
     *
     * This method is called once between each timestep after convergence.
     *
     * @param defgrd Deformation gradient
     * @param params Container for additional information
     * @param gp Gauss point
     * @param eleGID Global element identifier
     */
    void Update(CORE::LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    /*!
     * \brief Method that will be called once before the evaluation process. The elastin material
     * evaluates the prestress here.
     *
     * \param mixtureRule Mixture rule
     * \param params Container for additional information
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void PreEvaluate(
        MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID) override;


    double GetGrowthScalar(int gp) const override;

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
    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

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
    void EvaluateElasticPart(const CORE::LINALG::Matrix<3, 3>& F,
        const CORE::LINALG::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp,
        int eleGID) override;

    /*!
     * \brief Evaluation of the membrane stress only
     *
     * \param S 2. Piola Kirchhoff stress tensor of the membrane in stress-like Voigt notation
     * \param params Container for additional information
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void EvaluateMembraneStress(
        CORE::LINALG::Matrix<6, 1>& S, Teuchos::ParameterList& params, int gp, int eleGID) override;

   protected:
    /*!
     * \brief Evaluation of the stress and constitutive tensor of the membrane material
     *
     * \param F Total deformation gradient
     * \param iFin Inelastic deformation gradient
     * \param params Container for additional information
     * \param S_stress 2. Piola-Kirchhoff stress tensor in stress-like Voigt notation
     * \param cmat Constitutive tensor
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void EvaluateStressCMatMembrane(const CORE::LINALG::Matrix<3, 3>& F,
        const CORE::LINALG::Matrix<3, 3>& iFin, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp,
        int eleGID) const;

    /*!
     * \brief Evaluates the structural tensors of the radial and membrane plane direction in the
     * grown configuration
     *
     * \param Aradgr Structural tensor in radial direction in the grown configuration
     * \param Aorthgr Structural tensor in the membrane direction in the grown configuration
     * \param iFin Inverse of the inelastic deformation gradient
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void EvaluateStructuralTensorsInGrownConfiguration(CORE::LINALG::Matrix<3, 3>& Aradgr,
        CORE::LINALG::Matrix<3, 3>& Aorthgr, const CORE::LINALG::Matrix<3, 3>& iFin, int gp,
        int eleGID) const;

    /*!
     * \brief Evaluate the matrix Product \[
     *      A_{orth,gr} C_e A_{orth,gr} + A_{rad,gr}
     * \]
     *
     * \param AorthgrCeAorthgrArad (out) : desired product
     * \param Aradgr (in) : Structural tensor of the radial direction in the grown configuration
     * \param Aorthgr (in) : Structural tensor of the membrane plane direction in the grown
     * configuration
     *
     * \param Ce (in) : Elastic Cauchy-Green deformation tensor
     */
    static void EvaluateAorthgrCeAorthgrArad(CORE::LINALG::Matrix<3, 3>& AorthgrCeAorthgrArad,
        const CORE::LINALG::Matrix<3, 3>& Aradgr, const CORE::LINALG::Matrix<3, 3>& Aorthgr,
        const CORE::LINALG::Matrix<3, 3>& Ce);

    /*!
     * \brief Evaluate the matrix product \[
     *  F_{in}^{-1} A_{orth,gr} F_{in}^{-T}
     * \]
     *
     * \param iFinAorthgriFinT (out) : Desired product
     * \param iFin (in) : Inelastic inverse deformation gradient
     * \param Aorthgr (in) : Structural tensor of the membrane plane direction in grown
     * configuration
     */
    static void EvaluateiFinAorthgriFinT(CORE::LINALG::Matrix<3, 3>& iFinAorthgriFinT,
        const CORE::LINALG::Matrix<3, 3>& iFin, const CORE::LINALG::Matrix<3, 3>& Aorthgr);

    /*!
     * \brief Evaluates the prioduct \[
     *  F_{in}^{-T} A_{orth,gr}^T (A_{orth,gr} C_e A_{orth,gr} + A_{rad,gr})^{-1} A_{orth,gr}
     * F_{in}^{-1}
     * \]
     *
     * \param iFinTAorthgrTiXTAorthgriFin Desired product
     * \param AorthgrCeAorthgrArad Product \[A_{orth,gr} C_e A_{orth,gr} + A_{rad,gr}/]
     * \param iFin inverse of the inelastic part of the deformation
     * \param Aorthgr Structural tensor of the membrane plane direction in grown configuration
     */
    static void EvaluateiFinTAorthgrTiXTAorthgriFin(
        CORE::LINALG::Matrix<3, 3>& iFinTAorthgrTiXTAorthgriFin,
        const CORE::LINALG::Matrix<3, 3>& AorthgrCeAorthgrArad,
        const CORE::LINALG::Matrix<3, 3>& iFin, const CORE::LINALG::Matrix<3, 3>& Aorthgr);

   private:
    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_ElastHyperElastinMembrane* params_;

    /// Current growth factor with respect to the reference configuration
    std::vector<double> current_reference_growth_;

    /// fraction of membrane elastin material to ensure equilibrium
    std::vector<double> mue_frac_;

    /// map to membrane materials/potential summands (only IsoNeoHooke is possible)
    std::vector<Teuchos::RCP<MAT::ELASTIC::IsoNeoHooke>> potsum_membrane_;

    /// Anisotropy extension holding the structural tensor of the anisotropy
    ElastinMembraneAnisotropyExtension anisotropyExtension_;
  };

}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif
