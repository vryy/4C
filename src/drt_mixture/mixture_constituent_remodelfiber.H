/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a general remodel fiber constituent (basis for an explicit or implicit update
rule)

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_REMODELFIBER_H
#define BACI_MIXTURE_CONSTITUENT_REMODELFIBER_H

#include "mixture_constituent.H"
#include "growth_evolution.H"
#include "../drt_mat/anisotropy_extension.H"

// Forward declarationa
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
    class ActiveSummand;
  }  // namespace ELASTIC
}  // namespace MAT

namespace MIXTURE
{
  namespace PAR
  {
    class MixtureConstituent_RemodelFiber : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituent_RemodelFiber(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      ~MixtureConstituent_RemodelFiber() override = default;

      /// @name material parameters
      // @{

      /// Material ids of the strain energy function
      const int matid_;

      /// Decay time of Poisson (degradation) process (time constant)
      const double poisson_decay_time_;

      /// Time constant for growth of the fiber
      const double growth_constant_;

      /// deposition stretch of new fibers
      const double deposition_stretch_;

      // @}
    };
  }  // namespace PAR

  /*!
   * \brief Anisotropy extension that provides the structural tensor of the direction orthogonal to
   * the fiber direction: M = I - f*f^^T
   */
  class OrthogonalAnisotropyExtension : public MAT::FiberAnisotropyExtension<1>
  {
   public:
    /*!
     * \brief Constructor
     */
    explicit OrthogonalAnisotropyExtension();

    /*!
     * \brief Link to the anisotropy extension of the fiber direction
     *
     * \param other
     */
    void SetAnisotropyExtension(Teuchos::RCP<MAT::FiberAnisotropyExtension<1>>& other);

    ///@name Packing and Unpacking
    /// @{
    void PackAnisotropy(DRT::PackBuffer& data) const override;

    void UnpackAnisotropy(
        const std::vector<char>& data, std::vector<char>::size_type& position) override;
    /// @}

    /*!
     * \brief Returns a constant reference to the orthogonal structural tensor
     *
     * \param gp Gauss point
     * \return const LINALG::Matrix<3, 3>& Reference to the orthogonal structural tensor
     */
    const LINALG::Matrix<3, 3>& OrthogonalStructuralTensor(int gp) const;

    void OnGlobalDataInitialized() override;

   private:
    /// Anisitropy extension of the fiber summand
    Teuchos::RCP<MAT::FiberAnisotropyExtension<1>> other_;

    /// Storage for the structural tensors in orthogonal direction to the fiber direction
    std::vector<LINALG::Matrix<3, 3>> A_orth_;
  };

  /*!
   * \brief Base constituent for a remodeling fiber simulated with homogenized constrained mixtures
   */
  class MixtureConstituent_RemodelFiber : public MIXTURE::MixtureConstituent
  {
   public:
    /*!
     * \brief Construct remodel fiber by the input parameters
     *
     * \param params Input parameters
     */
    explicit MixtureConstituent_RemodelFiber(MIXTURE::PAR::MixtureConstituent_RemodelFiber* params);

    //! @copydoc MIXTURE::MixtureConstituent::PackConstituent(DRT::PackBuffer&)
    void PackConstituent(DRT::PackBuffer& data) const override;

    //! @copydoc MIXTURE::MixtureConstituent::UnpackConstituent(std::vector<char>::size_type&,const
    //! std::vector<char>&)
    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

    /*!
     * \brief Read element method to be passed to the summands
     *
     * \param numgp Number of Gauss points
     * \param linedef Line definition of the element
     */
    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /*!
     * \brief Setup called once per element
     *
     * \param params Setup Growth and remodeling parameters
     * \param eleGID global element id
     */
    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    /*!
     * \brief Extended update method.
     *
     * \param F Deformation gradient
     * \param params Container for additional information
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void Update(
        LINALG::Matrix<3, 3> const& F, Teuchos::ParameterList& params, int gp, int eleGID) override;

    /*!
     * @brief Evaluates the elastic part of the fiber
     *
     * @param F (in) : whole deformation gradient
     * @param Fextin (in) : external inelastic deformation gradient
     * @param params (in) : Container for additional information
     * @param S_stress (out) : 2nd Piola-Kirchhoff stress in stress-like Voigt notation
     * @param cmat (out) : Linearization of the material model
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    void EvaluateElasticPart(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& Fextin,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * \brief Default evaluation routine which is not implemented
     *
     * \param F total deformation gradient
     * \param E_strain Green-Lagrange strain tensor
     * \param params Container for additional information
     * \param S_stress 2. Piola-Kirchhoff stress tensor in stress-like Voigt notation
     * \param cmat Constitutitve tensor
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * \brief Returns the current reference mass density
     *
     * Due to growth, the reference mass density changes during the simulation.
     *
     * \param gp Gauss point
     * \return double current reference mass density
     */
    double CurrentRefDensity(int gp) const override;

   protected:
    /*!
     * @brief Integrates the internal growth and remodeling equations using an explicit euler
     * scheme. This method should be called once during the update call.
     *
     * @param F Deformation gradient of the finalized timestep
     * @param Fg Growth deformation gradient of the finalized timestep
     * @param dt Time increment to the next timestep
     * @param gp Current Gauss-Point
     * @param eleGID Global element id
     * @return
     */
    void UpdateGrowthAndRemodelingExpl(const LINALG::Matrix<3, 3>& F,
        const LINALG::Matrix<3, 3>& Fg, double dt, int gp, int eleGID);

    /*!
     * @brief Updates the internal variables of the fiber data
     *
     * @param dt time increment from last time step
     * @param gp Gauss point
     */
    void UpdateNewton(double dt, int gp);

    /*!
     * @brief Adds the stresses and the linearization to the vector
     *
     * @param FM (in) : deformation gradient from reference configuration to the current
     * configuration
     * @param iFinM (in) : inverse of the deformation gradient of the inelastic deformation in the
     * reference configuration
     * @param S_stress (out) : stress vector in stress notation where the 2nd Piola-Kirchhoff
     * stresses will be added
     * @param cmat (out) : linearization matrix, where the linearization of the fiber will be added
     * @param params (in) :  Container for additional information
     * @param gp (in) : Guass point
     * @param eleGID (in) : global element id
     */
    void AddStressCmat(const LINALG::Matrix<3, 3>& FM, const LINALG::Matrix<3, 3>& iFinM,
        LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, Teuchos::ParameterList& params,
        int gp, int eleGID);

    /*!
     * @brief Evaluates the first derivative of the free-energy function with respect to the right
     * Cauchy-Green deformation tensor
     * @param Ce Elastic part of the Right Cauchy-Green deformation tensor
     * @param iFin Inverse inelastic deformation gradeint
     * @param dPsidC First derivative of the free-energy function with respect to the right
     * Cauchy-Green deformation tensor
     * @param gp Gauss-Point
     * @param eleGID Global element id
     */
    void EvaluatedPsidC(const LINALG::Matrix<3, 3>& Ce, const LINALG::Matrix<3, 3>& iFin,
        LINALG::Matrix<3, 3>& dPsidC, int gp, int eleGID) const;

    /*!
     * @brief Evaluates the second derivative of the free-energy function with respect to the right
     * Cauchy-Green deformation tensor
     *
     * @param Ce Elastic part of the Right Cauchy-Green deformation tensor
     * @param iFin Inverse inelastic deformation gradient
     * @param ddPsidCdC_stress second derivative of the free-energy function with respect to the
     * right Cauchy-Green deformation tensor in stress-like Voigt notation
     * @param gp Gauss point
     * @param eleGID Global element id
     */
    void EvaluateddPsidCdC(const LINALG::Matrix<3, 3>& Ce, const LINALG::Matrix<3, 3>& iFin,
        LINALG::Matrix<6, 6>& ddPsidCdC_stress, int gp, int eleGID) const;

    /*!
     * \brief Evaluates the growth evolution equation and integrates it one timestep
     *
     * \param sig Cauchy fiber stress
     * \param gp Gauss point
     * \param eleGID Global element id
     * \return double growth evolution equation * dt
     */
    double EvaluateGrowthEvolutionEquationdt(double sig, int gp, int eleGID) const;

    /*!
     * \brief Evaluates the remodel evolution equation, integrated one timestep
     *
     * \param C Cauchy-Green deformation tensor
     * \param Fin inelastic deformation gradient
     * \param sig Cauchy fiber stress
     * \param gp Gauss point
     * \param eleGID Global element id
     * \return double remodel evolution equation * dt
     */
    double EvaluateRemodelEvolutionEquationdt(
        LINALG::Matrix<3, 3>& C, LINALG::Matrix<3, 3>& Fin, double sig, int gp, int eleGID) const;

    /*!
     * @brief Evaluates the Cauchy stress in the fiber (1D)
     *
     * @param C Right Cauchy-Green strain tensor
     * @param iFin inverse inelastic deformation gradient
     * @param gp Gauss-point
     * @param eleGID Global element id
     *
     * @return Cauchy normal stress in fiber direction (1D)
     */
    double EvaluateFiberCauchyStress(
        const LINALG::Matrix<3, 3>& C, const LINALG::Matrix<3, 3>& Fin, int gp, int eleGID) const;

    /*!
     * @brief Evaluates the elastic right Cauchy-Green deformation tensor
     *
     * @param C [in] :  Right Cauchy-Green deformation tensor
     * @param iFin [in] :  inverse inelastic deformation gradient
     * @param Ce [out] :  elastic right Cauchy-Green deformation tensor
     */
    static void EvaluateCe(
        const LINALG::Matrix<3, 3>& C, const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<3, 3>& Ce);

    /*!
     * \brief Evaluates the derivative of the Cauchy fiber stress w.r.t. the elastic Cauchy-Green
     * deformation tensor
     *
     * \param Ce Elastic Cauchy-Green deformation tensor
     * \param dsigdCe derivative of the Cauchy fiber stress w.r.t. the elastic Cauchy-Green
     * deformation tensor
     *
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void EvaluatedsigdCe(
        const LINALG::Matrix<3, 3>& Ce, LINALG::Matrix<3, 3>& dsigdCe, int gp, int eleGID) const;

    /*!
     * \brief Computes the homeostatic stress and stores it as the internal state
     */
    void UpdateSigH(int eleGID);

   private:
    /// Strain energy summand of the fiber
    Teuchos::RCP<MAT::ELASTIC::Summand> summand_;

    /// Strain energy the active summand of the fiber (May be Teuchos::null)
    Teuchos::RCP<MAT::ELASTIC::ActiveSummand> activeSummand_;

    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_RemodelFiber* params_;

    /// Vector holding the remodel deformation gradient in matrix notation
    std::vector<LINALG::Matrix<3, 3>> iFr_;

    /// Growth evolution equation
    Teuchos::RCP<GrowthEvolution> growth_evolution_;

    /// list of current densities of the fiber for each Gauss-point
    std::vector<double> cur_rho_;

    /// Current remodel parameter of the fiber for each Gauss-point
    std::vector<double> cur_lambdar_;

    /// Prestretch of the fiber for each Gauss-point
    std::vector<double> lambda_pre_;

    /// Homeostatic stress of the fiber (also deposition stress)
    std::vector<double> sig_h_;

    /// Anisotropy extension of the fiber direction
    Teuchos::RCP<MAT::FiberAnisotropyExtension<1>> fiberAnisotropyExtension_;

    /// Extension of the external fiber anisotropy extension
    OrthogonalAnisotropyExtension orthogonalAnisotropyExtension_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_REMODELFIBER_H
