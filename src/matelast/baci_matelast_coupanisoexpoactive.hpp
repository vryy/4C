/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the coupled contribution of an anisotropic active fiber material

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUPANISOEXPOACTIVE_HPP
#define FOUR_C_MATELAST_COUPANISOEXPOACTIVE_HPP

#include "baci_config.hpp"

#include "baci_mat_anisotropy_extension_default.hpp"
#include "baci_mat_anisotropy_extension_provider.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_matelast_activesummand.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for coupled contribution of an anisotropic active fiber material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_CoupAnisoExpoActive K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0
       * ADAPT_ANGLE 0 S 54000 LAMBDAM 1.4 LAMBDA0 0.8 DENS 1050
       */
      class CoupAnisoExpoActive : public MAT::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        explicit CoupAnisoExpoActive(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double k1_;
        double k2_;
        /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
        double gamma_;
        /// fiber params for the compressible case
        double k1comp_;
        double k2comp_;
        /// fiber initalization status
        int init_;
        /// adapt angle during remodeling
        bool adapt_angle_;
        /// maximum contractile stress
        double s_;
        /// stretch at maximum active force generation
        double lambdamax_;
        /// stretch at zero active force generation
        double lambda0_;
        /// total reference mass density at the beginning of the simulation
        double dens_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          dserror(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class CoupAnisoExpoActive
    }     // namespace PAR

    /*!
     * @brief Coupled anisotropic active fiber function, implemented for one possible fiber family
     * as in [1]
     *
     * This is an active anisotropic material of the most simple kind.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = \frac {s}{\rho} \left(1.0 + \frac{1}{3} \frac{\left(\lambda_m - 1.0
     *   \right)^3}{\left(\lambda_m - \lambda_0 \right)^2} \right).
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] Wilson, J.S., S. Baek, and J.D. Humphrey, Parametric study of effects of collagen
     * turnover on the natural history of abdominal aortic aneurysms. Proc. R. Soc. A, 2013.
     * 469(2150): p. 20120556.
     * </ul>
     */
    class CoupAnisoExpoActive : public ActiveSummand, public FiberAnisotropyExtensionProvider<1>
    {
     public:
      /// constructor with given material parameters
      explicit CoupAnisoExpoActive(MAT::ELASTIC::PAR::CoupAnisoExpoActive* params);

      ///@name Packing and Unpacking
      //@{

      void PackSummand(CORE::COMM::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;

      //@}

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_coupanisoexpoactive;
      }

      //@}

      /// Setup of active summand
      void Setup(int numgp, INPUT::LineDefinition* linedef) override;

      void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

      void EvaluateFirstDerivativesAniso(CORE::LINALG::Matrix<2, 1>& dPI_aniso,
          const CORE::LINALG::Matrix<3, 3>& rcg, int gp, int eleGID) override;

      void EvaluateSecondDerivativesAniso(CORE::LINALG::Matrix<3, 1>& ddPII_aniso,
          const CORE::LINALG::Matrix<3, 3>& rcg, int gp, int eleGID) override;

      /// retrieve coefficients of first, second and third derivative
      /// of summand with respect to anisotropic invariants
      /// ATTENTION: this is only the passive contribution of the fiber!
      template <typename T>
      void GetDerivativesAniso(CORE::LINALG::Matrix<2, 1, T>&
                                   dPI_aniso,  ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<3, 1, T>&
              ddPII_aniso,  ///< second derivative with respect to invariants
          CORE::LINALG::Matrix<4, 1, T>&
              dddPIII_aniso,  ///< third derivative with respect to invariants
          CORE::LINALG::Matrix<3, 3, T> const& rcg,  ///< right Cauchy-Green tensor
          int gp,                                    ///< Gauss point
          int eleGID) const;                         ///< element GID

      /// Add anisotropic principal stresses
      /// ATTENTION: this is only the passive contribution of the fiber!
      void AddStressAnisoPrincipal(const CORE::LINALG::Matrix<6, 1>&
                                       rcg,  ///< right Cauchy Green in "strain-like" Voigt notation
          CORE::LINALG::Matrix<6, 6>& cmat,  ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,  ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// Evaluates strain energy for automatic differentiation with FAD
      template <typename T>
      void EvaluateFunc(T& psi,                      ///< strain energy functions
          CORE::LINALG::Matrix<3, 3, T> const& rcg,  ///< Right Cauchy-Green tensor
          int gp,                                    ///< Gauss point
          int eleGID) const;                         ///< element GID

      /// evaluate stress and cmat
      /// ATTENTION: this is only the active contribution of the fiber!
      void AddActiveStressCmatAniso(
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< rigtht Cauchy Green tensor
          CORE::LINALG::Matrix<6, 6>& cmat,      ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,    ///< 2nd PK-stress
          int gp,                                ///< Gauss point
          int eleGID) const override;            ///< element GID


      /// evaluate stress and cmat
      /// ATTENTION: this is only the active contribution of the fiber!
      template <typename T>
      void EvaluateActiveStressCmatAniso(
          CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< rigtht Cauchy Green tensor
          CORE::LINALG::Matrix<6, 6, T>& cmat,      ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1, T>& stress,    ///< 2nd PK-stress
          int gp,                                   ///< Gauss point
          int eleGID) const;                        ///< element GID

      // add strain energy
      void AddStrainEnergy(double& psi,  ///< strain energy functions
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<6, 1>&
              glstrain,  ///< Green-Lagrange strain in strain like Voigt notation
          int gp,        ///< Gauss point
          int eleGID     ///< element GID
          ) override;

      /// @name Access methods
      //@{
      template <typename T>
      inline void GetDerivativeAnisoActive(T& dPIact) const
      {
        dPIact = dPIact_;
      };

      double GetDerivativeAnisoActive() const override { return dPIact_; };

      //@}

      /// Set fiber directions
      void SetFiberVecs(double newgamma,             ///< new angle
          const CORE::LINALG::Matrix<3, 3>& locsys,  ///< local coordinate system
          const CORE::LINALG::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        anisoprinc = true;
      };

      /*!
       * \brief Returns the reference to the used MAT::FiberAnisotropyExtension
       *
       * \return FiberAnisotropyExtension& Reference to the used MAT::AnisotropyExtension
       */
      FiberAnisotropyExtension<1>& GetFiberAnisotropyExtension() override
      {
        return anisotropyExtension_;
      }

     private:
      /// Evaluate the first derivative of the active fiber potential w.r.t active fiber stretch
      double EvaluatedPsiActive() const;

      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoExpoActive* params_;

      /// first derivative of active fiber potential w.r.t. the active fiber stretch
      double dPIact_;

      /// active fiber stretch for a given muscle tone
      double lambdaact_;

      /// Anisotropy extension that manages fibers and structural tensors
      DefaultAnisotropyExtension<1> anisotropyExtension_;
    };  // namespace PAR

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MATELAST_COUPANISOEXPOACTIVE_H
