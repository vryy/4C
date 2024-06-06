/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of classes for the exponential shear behavior for fibers

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_COUPANISOEXPOSHEAR_HPP
#define FOUR_C_MATELAST_COUPANISOEXPOSHEAR_HPP

#include "4C_config.hpp"

#include "4C_mat_anisotropy_extension_default.hpp"
#include "4C_mat_anisotropy_extension_provider.hpp"
#include "4C_matelast_coupanisoexpobase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace Elastic
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
      void pack_anisotropy(Core::Communication::PackBuffer& data) const override;

      void unpack_anisotropy(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      /// @}

      double GetScalarProduct(int gp) const override;
      const Core::LinAlg::Matrix<3, 3>& get_structural_tensor(int gp) const override;
      const Core::LinAlg::Matrix<6, 1>& get_structural_tensor_stress(int gp) const override;

      /*!
       * /copydoc
       *
       * The coupling structural tensor and the scalar product will be computed here
       */
      void on_global_data_initialized() override;

     protected:
      void on_global_element_data_initialized() override;
      void on_global_gp_data_initialized() override;

     private:
      /**
       * Scalar products of the fibers at the Gauss points
       */
      std::vector<double> scalar_products_;

      /**
       * Coupling structural tensor of the fibers in stress like Voigt notation at the Gauss points
       */
      std::vector<Core::LinAlg::Matrix<6, 1>> structural_tensors_stress_;

      /**
       * Coupling structural tensor of the fibers at the Gauss points
       */
      std::vector<Core::LinAlg::Matrix<3, 3>> structural_tensors_;

      /// Flag whether fibers are initialized
      bool is_initialized_{};

      /// Initialization mode
      const int init_mode_;

      /// Fiber ids to be used to build the shear behavior
      const std::array<int, 2> fiber_ids_;
    };

    namespace PAR
    {
      /*!
       * @brief material parameters for coupled contribution of a anisotropic exponential fiber
       * material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_CoupAnisoExpoShear K1 10.0 K2 1.0 K1COMP 0.0 K2COMP 1.0 INIT 0 FIBER_IDS 1 2
       */
      class CoupAnisoExpoShear : public Core::Mat::PAR::Parameter,
                                 public Mat::Elastic::PAR::CoupAnisoExpoBase
      {
       public:
        /// standard constructor
        explicit CoupAnisoExpoShear(const Teuchos::RCP<Core::Mat::PAR::Material>& matdata);

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<Core::Mat::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "Mat::Elastic::Summand::Factory.");
          return Teuchos::null;
        };

        /// @name material parameters
        //@{
        /// Ids of the fiber for the shear behavior
        std::array<int, 2> fiber_id_{};
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
    class CoupAnisoExpoShear : public Mat::Elastic::CoupAnisoExpoBase
    {
     public:
      /// constructor with given material parameters
      explicit CoupAnisoExpoShear(Mat::Elastic::PAR::CoupAnisoExpoShear* params);

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType MaterialType() const override
      {
        return Core::Materials::mes_coupanisoexposhear;
      }

      //@}

      /// @name Methods for Packing and Unpacking
      ///@{
      void PackSummand(Core::Communication::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      ///@}

      /*!
       * \brief Register the anisotropy extension to the global anisotropy manager
       *
       * \param anisotropy anisotropy manager
       */
      void register_anisotropy_extensions(Mat::Anisotropy& anisotropy) override;

      /// Set fiber directions
      void SetFiberVecs(double newgamma,             ///< new angle
          const Core::LinAlg::Matrix<3, 3>& locsys,  ///< local coordinate system
          const Core::LinAlg::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

     protected:
      const CoupAnisoExpoBaseInterface& get_coup_aniso_expo_base_interface() const override
      {
        return anisotropy_extension_;
      }

     private:
      /// my material parameters
      Mat::Elastic::PAR::CoupAnisoExpoShear* params_;

      /// Internal ansotropy information
      Mat::Elastic::CoupAnisoExpoShearAnisotropyExtension anisotropy_extension_;
    };

  }  // namespace Elastic
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
