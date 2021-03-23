/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the Weickenmeier active skeletal muscle material model constituent

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_MUSCLE_WEICKENMEIER_H
#define BACI_MIXTURE_CONSTITUENT_MUSCLE_WEICKENMEIER_H

#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/anisotropy.H"
#include "../drt_mat/anisotropy_extension_default.H"
#include "../drt_mat/anisotropy_extension_provider.H"
#include "mixture_constituent.H"
#include <Teuchos_RCP.hpp>
namespace MIXTURE
{
  namespace PAR
  {
    class MixtureConstituent_Muscle_Weickenmeier : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      /// constructor
      explicit MixtureConstituent_Muscle_Weickenmeier(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::MixtureConstituent> CreateConstituent(int id) override;

      /// @name material parameters
      //@{
      //! @name passive material parameters
      const double alpha_;   ///< material parameter, >0
      const double beta_;    ///< material parameter, >0
      const double gamma_;   ///< material parameter, >0
      const double kappa_;   ///< material parameter for coupled volumetric contribution
      const double omega0_;  ///< weighting factor for isotropic tissue constituents, governs ratio
                             ///< between muscle matrix material (omega0) and muscle fibers (omegap)
                             ///< with omega0 + omegap = 1
                             //! @}

      //! @name active microstructural parameters
      //! @name stimulation frequency dependent activation contribution
      const double
          Na_;  ///< number of active motor units (MU) per undeformed muscle cross-sectional area
      const int muTypesNum_;  ///< number of motor unit (MU) types
                              ///< vectors indicating corresponding parameters of motor unit types
                              ///< e.g. slow type I fibres (index 0), fast resistant  type IIa
                              ///< fibres (index 1), fast fatigue type IIb fibres (index 2)
      const std::vector<double> I_;    ///< interstimulus interval
      const std::vector<double> rho_;  ///< fraction of motor unit types
      const std::vector<double> F_;    ///< twitch force
      const std::vector<double> T_;    ///< twitch contraction time
      //! @}

      //! @name stretch dependent activation contribution
      const double lambdaMin_;  ///< minimal active fiber stretch
      const double
          lambdaOpt_;  ///< optimal active fiber stretch related active nominal stress maximimum
      //! @}

      //! @name velocity dependent activation contribution
      //! @{
      const double dotLambdaMMin_;  ///< minimal stretch rate
      const double ke_;  ///< dimensionless constant controlling the curvature of the velocity
                         ///< dependent activation function in the eccentric case
      const double kc_;  ///< dimensionless constant controlling the curvature of the velocity
                         ///< dependent activation function in the concentric case
      //! @}
      //! @}

      //! @name prescribed activation in corresponding time intervals
      //! @{
      const int actTimesNum_;                ///< number of time boundaries
      const std::vector<double> actTimes_;   ///< time boundaries between intervals
      const int actIntervalsNum_;            ///< number of time intervals
      const std::vector<double> actValues_;  ///< scaling factor in intervals
                                             ///< (1=full activation, 0=no activation)
      //! @}
      //@}

    };  // end class MixtureConstituent_Muscle_Weickenmeier
  }     // end namespace PAR

  /*!
   * \brief Constituent for the Weickenmeier muscle material
   *
   * This constituent represents the Weickenmeier muscle material. It has to
   * be paired with the MAT::Mixture_ElastHyper material and a MIXTURE::MixtureRule.
   *
   * Reference for the material model: J. Weickenmeier, M. Itskov, E Mazza and M. Jabareen, ‘A
   * physically motivated constitutive model for 3D numerical simulation of skeletal muscles’,
   * International journal for numerical methods in biomedical engineering, vol. 30, no. 5, pp.
   * 545–562, 2014, doi: 10.1002/cnm.2618.
   */
  class MixtureConstituent_Muscle_Weickenmeier : public MIXTURE::MixtureConstituent
  {
   public:
    /// Constructor for the material given the material parameters
    explicit MixtureConstituent_Muscle_Weickenmeier(
        MIXTURE::PAR::MixtureConstituent_Muscle_Weickenmeier* params, int id);

    // Pack this class so it can be communicated
    void PackConstituent(DRT::PackBuffer& data) const override;

    // Unpack data from a char vector into this class
    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    // Material type
    INPAR::MAT::MaterialType MaterialType() const override;

    // Returns the reference mass density
    double CurrentRefDensity(int gp) const override;

    // Initialize the constituent with the parameters of the input line
    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    // Method that is executed before the first evaluate call, once for each Gauss point
    void PreEvaluate(
        MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID) override;

    // Evaluate stress and material tangent of the constituent and add to mixture stress and cmat
    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

   private:
    /*!
     *  \brief Evaluate active nominal stress Pa and its derivative w.r.t. the fiber stretch
     *
     *  \param[in]     params Container for additional information
     *  \param[in]     lambdaM Fiber stretch
     *  \param[in,out] Pa Active nominal stress
     *  \param[in,out] derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     */
    /// Evaluate active nominal stress Pa and its derivative w.r.t. the fiber stretch
    void EvaluateActiveNominalStress(
        Teuchos::ParameterList& params, const double lambdaM, double& Pa, double& derivPa);

    /*!
     *  \brief Evaluate activation level omegaa and its derivative w.r.t. the right Cauchy Green
     * tensor
     *
     *  \param[in]     params Container for additional information
     *  \param[in]     lambdaM Fiber stretch
     *  \param[in]     M Fiber structural tensor (dyadic product of fiber directions)
     *  \param[in]     Pa Active nominal stress
     *  \param[in]     derivPa Derivative of active nominal stress w.r.t. the fiber stretch
     *  \param[in,out] omegaa Activation level
     *  \param[in,out] domegaadC Derivative of the activation level w.r.t. the right Cauchy
     *                           Green tensor
     */
    void EvaluateActivationLevel(Teuchos::ParameterList& params, const double lambdaM,
        LINALG::Matrix<3, 3>& M, double Pa, double derivPa, double& omegaa,
        LINALG::Matrix<3, 3>& domegaadC);

    /*!
     *  \brief Evaluate Lambert W function with Halley's method
     *
     *  Solution of Lambert W function is functional inverse of xi = W_0*exp(W_0)
     *  Computation here restricted to principal branch W_0
     *  Use of Halley's method according to:
     *  https://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-funsolutiction/
     *
     *  \param[in]     xi Argument of Lambert W function W(xi)
     *  \param[in,out] WO Solution of principal branch of Lambert W function
     *  \param[in]     tol Tolerance for Halley's approximation
     *  \param[in]     maxiter Maximal number of iterations
     */
    void EvaluateLambert(double xi, double& W0, double tol, double maxiter);

    /*!
     *  \brief Add scaled derivative of invC*L*invC to Cmat
     *
     *  Add the following contribution to the constitutive tensor cmat(6,6) based on
     *  - the inverse of the right Cauchy-Green vector invC
     *  - the term invC*L*invC (\f$ \mathbf{C}^{-1} \mathbf{L} \mathbf{C}^{-1} \f$) with the
     *  structural tensor L (\f$ \mathbf{L} \f$)
     *
     *  \f$\text{scalar} \cdot \frac{\partial \mathbf{C}^{-1} \mathbf{L} \mathbf{C}^{-1}}{\partial
     *  \mathbf{C}} \f$
     *
     *  wherein the derivative \f$ \frac{\partial \mathbf{C}^-1 \mathbf{L} \mathbf{C}^-1}{\partial
     *  \mathbf{C}}  \f$ is computed to:  \f$ - 1/2 \cdot ( \mathbf{C}^{-1}_{ik}
     *  \mathbf{C}^{-1}_{jm} L_{mn} \mathbf{C}^{-1}_{nl} + \mathbf{C}^{-1}_{il} \mathbf{C}^{-1}_{jm}
     *  L_{mn} \mathbf{C}^{-1}_{nk} + \mathbf{C}^{-1}_{jk} \mathbf{C}^{-1}_{im} L_{mn}
     *  \mathbf{C}^{-1}_{nl} + \mathbf{C}^{-1}_{jl} \mathbf{C}^{-1}_{im} L_{mn}
     * \mathbf{C}^{-1}_{nk}) \f$
     *
     *  \param[in,out] cmat Constitutive matrix
     *  \param[in]     invC Inverse of right Cauchy Green tensor in matrix notation
     *  \param[in]     invCLinvC Product invC*L*invC with modified structual tensor and inverse of
     *                           right Cauchy Green tensor
     *  \param[in]     scalar scaling factor
     */
    void AddtoCmatDerivInvCLInvCProduct(LINALG::Matrix<6, 6>& cmat,
        const LINALG::Matrix<3, 3>& invC, const LINALG::Matrix<3, 3>& invCLinvC, double scalar);

    /// Weickenmeier material parameters
    MIXTURE::PAR::MixtureConstituent_Muscle_Weickenmeier* params_{};

    /// Total simulation time
    double t_tot_;

    /// Holder for anisotropic behavior
    MAT::Anisotropy anisotropy_;

    /// Anisotropy extension holder
    MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;
  };  // end class MixtureConstituent_Muscle_Weickenmeier

}  // end namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_MUSCLE_WEICKENMEIER_H
