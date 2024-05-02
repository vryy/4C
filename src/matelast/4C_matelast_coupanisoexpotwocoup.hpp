/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the passive material behaviour of cardiac muscle
according to Holzapfel and Ogden, "Constitutive modelling of passive myocardium", 2009.

\level 2
*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MATELAST_COUPANISOEXPOTWOCOUP_HPP
#define FOUR_C_MATELAST_COUPANISOEXPOTWOCOUP_HPP

#include "4C_config.hpp"

#include "4C_mat_anisotropy.hpp"
#include "4C_mat_anisotropy_extension_default.hpp"
#include "4C_mat_par_aniso.hpp"
#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for coupled passive cardiac material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_CoupAnisoExpoTwoCoup A4 18.472 B4 16.026 A6 2.481 B6 11.120 A8 0.216 B8 11.436
       * GAMMA 0.0 [INIT 1] [FIB_COMP Yes] [ADAPT_ANGLE No]
       */
      class CoupAnisoExpoTwoCoup : public MAT::PAR::ParameterAniso
      {
       public:
        /// constructor with given material parameters
        explicit CoupAnisoExpoTwoCoup(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double A4_;
        double B4_;
        double A6_;
        double B6_;
        double A8_;
        double B8_;
        /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
        double gamma_;
        /// fiber initalization status
        int init_;
        /// fibers support compression - or not
        bool fib_comp_;
        /// adapt angle during remodeling
        bool adapt_angle_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class CoupAnisoExpoTwoCoup

    }  // namespace PAR

    /*!
     * \brief Anisotropy manager for two fibers and the structural tensor of the combination of
     * those two
     */
    class CoupAnisoExpoTwoCoupAnisoExtension : public DefaultAnisotropyExtension<2>
    {
     public:
      /*!
       * \brief Constructor
       *
       * \param params Material parameters
       */
      explicit CoupAnisoExpoTwoCoupAnisoExtension(MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params);

      /*!
       * \brief Pack all data for parallel distribution and restart
       *
       * \param data data array to pack to
       */
      void PackAnisotropy(CORE::COMM::PackBuffer& data) const override;

      /*!
       * \brief Unpack data from the pack from parallel distribution and restart
       *
       * \param data data array to unpack from
       * \param position position of the data
       */
      void UnpackAnisotropy(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;

      /*!
       * \brief Notifier method when fibers are initialized.
       */
      void OnFibersInitialized() override;

      /*!
       * \brief Returns the reference to the coupled structural tensor in stress like Voigt
       * notation
       *
       * \param gp Gauss point
       * \return const CORE::LINALG::Matrix<6, 1>& Reference to the coupled structural tensor in
       * stress like Voigt notation
       */
      const CORE::LINALG::Matrix<6, 1>& GetCoupledStructuralTensor_stress(int gp) const;

      /*!
       * \brief Returns the coupled scalar product at the Gauss point
       *
       * \param gp Gauss point
       * \return double Scalar product of the two fibers
       */
      double GetCoupledScalarProduct(int gp) const;

     private:
      /// dot product fiber direction
      std::vector<double> a1a2_;

      /// mixed structural tensor (symmetric) \f$\frac{1}{2}(a1 \otimes a2 + a2 \otimes a1)\f$ in
      /// stress like Voigt notation
      std::vector<CORE::LINALG::Matrix<6, 1>> a1_a2_;
    };

    /*!
     * @brief Anisotropic cardiac material, implemented with two possible fiber families as in [1]
     *
     * This is a hyperelastic, anisotropic material for the passive response of cardiac material
     *
     * Strain energy function is given by:
     * \f[
     *   \Psi = \frac {a_4}{2 b_4} \left( e^{b_4 (IV_{\boldsymbol C} - 1)^2} - 1 \right) +
     *   \frac {a_6}{2 b_6} \left( e^{b_6 (VI_{\boldsymbol C} - 1)^2} - 1 \right) + \frac{a_8}{2b_8}
     * \left( e^{b_8 \left( VIII_{\boldsymbol C} - a_0 \cdot b_0 \right)^2} - 1\right) \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, RW Ogden, Constitutive modelling of passive myocardium: a
     * structurally based framework for material characterization
     * <li> [2] C Sansour, On the physical assumptions underlying the volumetric-isochoric split
     * and the case of anisotropy
     * </ul>
     */
    class CoupAnisoExpoTwoCoup : public Summand
    {
     public:
      /// constructor with given material parameters
      explicit CoupAnisoExpoTwoCoup(MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params);

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_coupanisoexpotwocoup;
      }
      //@}

      /*!
       * \brief Register the local anisotropy extension to the global anisotropy manager
       *
       * \param anisotropy Reference to the global anisotropy manager
       */
      void RegisterAnisotropyExtensions(Anisotropy& anisotropy) override;

      /// @name Methods for Packing and Unpacking
      ///@{
      void PackSummand(CORE::COMM::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      ///@}

      /// Add anisotropic principal stresses
      void AddStressAnisoPrincipal(
          const CORE::LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          CORE::LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

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

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params_;

      /// Special anisotropic behavior
      CoupAnisoExpoTwoCoupAnisoExtension anisotropy_extension_;
    };  // namespace PAR

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
