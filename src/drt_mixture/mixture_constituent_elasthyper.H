/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a hyperelastic constituent

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_ELASTHYPER_H
#define BACI_MIXTURE_CONSTITUENT_ELASTHYPER_H

#include "../drt_mat/matpar_parameter.H"
#include "mixture_constituent_elasthyperbase.H"
#include "../drt_mat/elasthyper_service.H"
#include "mixture_prestress_strategy.H"

namespace MIXTURE
{
  class MixtureConstituent_ElastHyper;

  namespace PAR
  {
    class MixtureConstituent_ElastHyper : public MIXTURE::PAR::MixtureConstituent_ElastHyperBase
    {
     public:
      explicit MixtureConstituent_ElastHyper(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      ~MixtureConstituent_ElastHyper() override = default;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::MixtureConstituent> CreateConstituent() override;

      /// @name material parameters
      /// @{
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Constituent for any hyperelastic material
   *
   * This constituent represents any hyperelastic material from the elasthyper toolbox. It has to
   * be paired with the MAT::Mixture_ElastHyper material and a MIXTURE::MixtureRule.
   */
  class MixtureConstituent_ElastHyper : public MIXTURE::MixtureConstituent_ElastHyperBase
  {
   public:
    /// Constructor for the materiak given the material parameters
    explicit MixtureConstituent_ElastHyper(MIXTURE::PAR::MixtureConstituent_ElastHyper* params);

    /// Returns the material type enum
    INPAR::MAT::MaterialType MaterialType() const override;

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauß point
     *
     * @param F Deformation gradient
     * @param E_strain Green-Lagrange strain in strain-like Voigt notation
     * @param params Container for additional information
     * @param S_stress 2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat Constitutive tensor in Voigt notation
     */
    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * @brief Evaluates the stress and material linearization of the constituents with an inelastic
     * part of the deformation
     *
     * The total deformation is #F, which is split into two parts:
     *
     * $\boldsymbol{F} = \boldsymbol{F}_e \cdot \boldsymbol{F}_in$
     *
     * Only elastic part $\boldsymbol{F}_e$ causes stresses. The inelastic part is only needed
     * for the linearization.
     *
     * @param F Total deformation gradient
     * @param iF_in Inverse of inelastic part of the deformation
     * @param params Container for additional information
     * @param S_stress 2nd Piola-Kirchhoff stress in stress-like Voigt notation
     * @param cmat Linearization of the material tensor in Voigt notation
     * @param gp Gauss-point
     * @param eleGID Global element id
     */
    void EvaluateElasticPart(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iF_in,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * \brief Returns the current reference mass density of the constituent
     *
     * @param gp (in) : Gauss point
     * @return
     */
    double CurrentRefDensity(int gp) const override;

   private:
    MIXTURE::PAR::MixtureConstituent_ElastHyper* params_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_ELASTHYPER_H
