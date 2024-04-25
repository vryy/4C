/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a hyperelastic constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_ELASTHYPER_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_ELASTHYPER_HPP

#include "4C_config.hpp"

#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_par_parameter.hpp"
#include "4C_mixture_constituent_elasthyperbase.hpp"
#include "4C_mixture_prestress_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  class MixtureConstituentElastHyper;

  namespace PAR
  {
    class MixtureConstituentElastHyper : public MIXTURE::PAR::MixtureConstituentElastHyperBase
    {
     public:
      explicit MixtureConstituentElastHyper(const Teuchos::RCP<MAT::PAR::Material>& matdata);
      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> CreateConstituent(int id) override;

      /// @name material parameters
      /// @{
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Constituent for any hyperelastic material
   *
   * This constituent represents any hyperelastic material from the elasthyper toolbox. It has to
   * be paired with the MAT::Mixture material and a MIXTURE::MixtureRule.
   */
  class MixtureConstituentElastHyper : public MIXTURE::MixtureConstituentElastHyperBase
  {
   public:
    /// Constructor for the materiak given the material parameters
    explicit MixtureConstituentElastHyper(
        MIXTURE::PAR::MixtureConstituentElastHyper* params, int id);

    /// Returns the material type enum
    INPAR::MAT::MaterialType MaterialType() const override;

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauss point
     *
     * @param F Deformation gradient
     * @param E_strain Green-Lagrange strain in strain-like Voigt notation
     * @param params Container for additional information
     * @param S_stress 2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat Constitutive tensor in Voigt notation
     */
    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

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
    void EvaluateElasticPart(const CORE::LINALG::Matrix<3, 3>& F,
        const CORE::LINALG::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp,
        int eleGID) override;
  };

}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
