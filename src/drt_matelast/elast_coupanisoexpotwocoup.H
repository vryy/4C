/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the routines required for the passive behaviour of cardiac muscle
according to Holzapfel and Ogden, "Constitutive modelling of passive myocardium", 2009.

\maintainer Amadeus Gebauer

\level 2
*/
/*---------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPOTWOCOUP_H_
#define ELAST_COUPANISOEXPOTWOCOUP_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/anisotropy.H"
#include "../drt_mat/anisotropy_extension_default.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for coupled passive cardiac material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_CoupAnisoExpoTwoCoup A4 18.472 B4 16.026 A6 2.481 B6 11.120 A8 0.216 B8 11.436
      /// GAMMA 0.0 [INIT 1] [FIB_COMP Yes] [ADAPT_ANGLE No]
      class CoupAnisoExpoTwoCoup : public MAT::PAR::ParameterAniso
      {
       public:
        /// constructor with given material parameters
        explicit CoupAnisoExpoTwoCoup(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// destructor
        ~CoupAnisoExpoTwoCoup() override = default;

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
        int fib_comp_;
        /// adapt angle during remodeling
        int adapt_angle_;

        //@}

        /// create material instance of matching type with my parameters
        ///
        /// Remark:                                                  AN 03/2014
        /// Function has to be implemented, since it is a pure virtual
        /// function in the MAT::PAR::Material class.
        /// Function stays empty, since the call within the material factory is done only
        /// respective the elasthyper class
        ///
        Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; }

        /// @name Write methods
        //@{
        void SetA4(const double A4) { A4_ = A4; }
        void SetB4(const double B4) { B4_ = B4; }
        void SetA6(const double A6) { A6_ = A6; }
        void SetB6(const double B6) { B6_ = B6; }
        void SetA8(const double A8) { A8_ = A8; }
        void SetB8(const double B8) { B8_ = B8; }
        //@}

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
      void PackAnisotropy(DRT::PackBuffer& data) const override;

      /*!
       * \brief Unpüack data from the pack from parallel distribution and restart
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
       * \brief Returns the reference to the coupled structural tensor in stress like Voigt notation
       *
       * \param gp Gauss point
       * \return const LINALG::Matrix<6, 1>& Reference to the coupled structural tensor in stress
       * like Voigt notation
       */
      const LINALG::Matrix<6, 1>& GetCoupledStructuralTensor_stress(int gp) const;

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
      std::vector<LINALG::Matrix<6, 1>> A1A2_;
    };

    /*----------------------------------------------------------------------*/
    /// Anisotropic cardiac material, implemented with two possible fiber families acc [1]
    ///
    /// This is a hyperelastic, anisotropic material for the passive response of cardiac material
    ///
    /// Strain energy function is given by:
    ///\f[
    ///   \Psi = \frac {a_4}{2 b_4} \left( e^{b_4 (IV_{\boldsymbol C} - 1)^2} - 1 \right) +
    ///   \frac {a_6}{2 b_6} \left( e^{b_6 (VI_{\boldsymbol C} - 1)^2} - 1 \right) + \frac {a_8}{2
    ///   b_8} \left( e^{b_8 \left( VIII_{\boldsymbol C} - a_0 \cdot b_0 \right)^2} - 1 \right)
    ///\f]
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] GA Holzapfel, RW Ogden, Constitutive modelling of passive myocardium: a
    /// structurally based framework for material characterization
    /// </ul>#include "../drt_mat/matpar_parameter.H"
    /// <ul>
    /// <li> [2] C Sansour, On the physical assumptions underlying the volumetric-isochoric split
    /// and the case of anisotropy
    /// </ul>
    ///
    /// \author nagler
    /// \date 11/11
    class CoupAnisoExpoTwoCoup : public Summand
    {
     public:
      /// constructor with given material parameters
      explicit CoupAnisoExpoTwoCoup(MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_coupanisoexpotwocoup;
      }
      //@}

      /*!
       * \brief Register the local anisotropy extension to the global anisotropy manager
       *
       * \param anisotropy Reference to the global anisotropy manager
       */
      void RegisterAnisotropyExtensions(Anisotropy& anisotropy) override;

      /// Add anisotropic principal stresses
      void AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          LINALG::Matrix<6, 6>& cmat,                                ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,                              ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;
      /// @name Access methods
      //@{
      double A4() const { return params_->A4_; }
      double B4() const { return params_->B4_; }
      double A6() const { return params_->A6_; }
      double B6() const { return params_->B6_; }
      double A8() const { return params_->A8_; }
      double B8() const { return params_->B8_; }
      //@}

      /// Set fiber directions
      void SetFiberVecs(double newgamma,       ///< new angle
          const LINALG::Matrix<3, 3>& locsys,  ///< local coordinate system
          const LINALG::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
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
      CoupAnisoExpoTwoCoupAnisoExtension anisotropyExtension_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // COUPANISOEXPOTWOCOUP_H_
