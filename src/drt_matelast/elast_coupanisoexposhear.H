/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of classes for the exponential shear behavior for fibers

\level 3

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPOSHEAR_H_
#define ELAST_COUPANISOEXPOSHEAR_H_

#include "../drt_matelast/elast_coupanisoexpobase.H"
#include "../drt_mat/anisotropy_extension_default.H"
#include "../drt_mat/anisotropy_extension_provider.H"


namespace MAT
{
  namespace ELASTIC
  {
    /*!
     * \brief Container class for communication with the base material
     */
    class CoupAnisoExpoShearAnisotropyExtension : public BaseAnisotropyExtension,
                                                  public CoupAnisoExpoBaseInterface
    {
     public:
      /*!
       * \brief Constructor
       *
       * \param init_mode initialization mode of the fibers
       * \param gamma Angle of the fiber if they are given in a local coordinate system
       * \param adapt_angle boolean, whether the fiber is subject to growth and remodeling
       * \param structuralTensorStrategy Strategy to compute the structural tensor
       * \param fiber_ids Ids of the fiber to be used for shear behavior
       */
      CoupAnisoExpoShearAnisotropyExtension(int init_mode, std::array<int, 2> fiber_ids);

      ///@name Packing and Unpacking
      /// @{
      void PackAnisotropy(DRT::PackBuffer& data) const override;

      void UnpackAnisotropy(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      /// @}

      double GetScalarProduct(int gp) const override;
      const LINALG::Matrix<3, 3>& GetStructuralTensor(int gp) const override;
      const LINALG::Matrix<6, 1>& GetStructuralTensor_stress(int gp) const override;

      /*!
       * /copydoc
       *
       * The coupling structural tensor and the scalar product will be computed here
       */
      void OnGlobalDataInitialized() override;

     protected:
      void OnGlobalElementDataInitialized() override;
      void OnGlobalGPDataInitialized() override;

     private:
      /**
       * Scalar products of the fibers at the Gauss points
       */
      std::vector<double> scalarProducts_;

      /**
       * Coupling structural tensor of the fibers in stress like Voigt notation at the Gauss points
       */
      std::vector<LINALG::Matrix<6, 1>> structuralTensors_stress_;

      /**
       * Coupling structural tensor of the fibers at the Gauss points
       */
      std::vector<LINALG::Matrix<3, 3>> structuralTensors_;

      /// Flag whether fibers are initialized
      bool isInitialized_;

      /// Initialization mode
      const int init_mode_;

      /// Fiber ids to be used to build the shear behavior
      const std::array<int, 2> fiber_ids_;
    };

    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for coupled contribution of a anisotropic exponential fiber
      /// material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_CoupAnisoExpoShear K1 10.0 K2 1.0 K1COMP 0.0 K2COMP 1.0 INIT 0 FIBER_IDS 1 2
      class CoupAnisoExpoShear : public MAT::PAR::Parameter,
                                 public MAT::ELASTIC::PAR::CoupAnisoExpoBase
      {
       public:
        /// standard constructor
        explicit CoupAnisoExpoShear(Teuchos::RCP<MAT::PAR::Material> matdata);

        Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

        /// @name material parameters
        //@{
        /// Ids of the fiber for the shear behavior
        std::array<int, 2> fiber_id_;
        //@}

      };  // class CoupAnisoExpoShear

    }  // namespace PAR

    /*!
     * \brief Exponential shear behavior between two fibers
     *
     * The strain energy function of this summand is
     * \[\psi = \frac{a_{fs}}{2b_{fs}} \left[ \exp( b_{fs} (I_{8fs} - \boldsymbol{f} \cdot
     * \boldsymbol{s})^2 ) - 1 \right]\]
     */
    class CoupAnisoExpoShear : public MAT::ELASTIC::CoupAnisoExpoBase
    {
     public:
      /// constructor with given material parameters
      explicit CoupAnisoExpoShear(MAT::ELASTIC::PAR::CoupAnisoExpoShear* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_coupanisoexposhear;
      }

      //@}

      /// @name Methods for Packing and Unpacking
      ///@{
      void PackSummand(DRT::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      ///@}

      /*!
       * \brief Register the anisotropy extension to the global anisotropy manager
       *
       * \param anisotropy anisotropy manager
       */
      void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

      /// Set fiber directions
      void SetFiberVecs(double newgamma,       ///< new angle
          const LINALG::Matrix<3, 3>& locsys,  ///< local coordinate system
          const LINALG::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

     protected:
      const CoupAnisoExpoBaseInterface& GetCoupAnisoExpoBaseInterface() const override
      {
        return anisotropyExtension_;
      }

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoExpoShear* params_;

      /// Internal ansotropy information
      MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension anisotropyExtension_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ELAST_COUPANISOEXPOSHEAR_H_
