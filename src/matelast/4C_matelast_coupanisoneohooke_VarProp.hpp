/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the coupled anisotropic neo-Hooke material with one fiber direction
with space-time varying coefficients

\level 3
*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MATELAST_COUPANISONEOHOOKE_VARPROP_HPP
#define FOUR_C_MATELAST_COUPANISONEOHOOKE_VARPROP_HPP

#include "4C_config.hpp"

#include "4C_mat_par_aniso.hpp"
#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace Elastic
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for anisotropic contribution of a neo-Hooke material with one
       * fiber direction with space-time varying coefficients
       *
       * <h3>Input line</h3>
       * MAT 1 CoupAnisoNeoHooke_VarProp C 100 GAMMA 35.0 INIT 0 ADAPT_ANGLE 0
       */
      class CoupAnisoNeoHookeVarProp : public Mat::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        CoupAnisoNeoHookeVarProp(const Core::Mat::PAR::Parameter::Data& matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double c_;
        /// Where the activation comes from: 0=scatra , >0 Id for FUNCT
        int sourceactiv_;
        /// azimute angle of spherical coordinates
        double gamma_;
        /// polar angle of spherical coordinates
        double theta_;
        /// fiber initalization status
        int init_;
        /// adapt angle during remodeling
        bool adapt_angle_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<Core::Mat::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "Mat::Elastic::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class CoupAnisoNeoHooke_VarProp

    }  // namespace PAR

    /*!
     * @brief Coupled anisotropic exponential fiber function, implemented for one possible fiber
     * family as in [1]
     *
     * This is a hyperelastic, anisotropic material of the most simple kind. The material parameter
     * can vary in space and time
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = c(X,t) \left(IV_{\boldsymbol C}-1\right)
     * \f]
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, TC Gasser, A viscoelastic model for fiber-reinforced composites at
     * finite strains: ... 2000
     * </ul>
     */
    class CoupAnisoNeoHookeVarProp : public Summand
    {
     public:
      /// empty constructor
      //    CoupAnisoNeoHooke_VarProp();

      /// constructor with given material parameters
      CoupAnisoNeoHookeVarProp(Mat::Elastic::PAR::CoupAnisoNeoHookeVarProp* params);

      ///@name Packing and Unpacking
      //@{

      void pack_summand(Core::Communication::PackBuffer& data) const override;

      void unpack_summand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;

      //@}

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType material_type() const override
      {
        return Core::Materials::mes_coupanisoneohooke_varprop;
      }

      //@}

      /// Setup of summand
      void setup(int numgp, Input::LineDefinition* linedef) override;

      /// Add anisotropic principal stresses
      void add_stress_aniso_principal(
          const Core::LinAlg::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          Core::LinAlg::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          Core::LinAlg::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// Set fiber directions
      void set_fiber_vecs(const double newgamma,     ///< new angle
          const Core::LinAlg::Matrix<3, 3>& locsys,  ///< local coordinate system
          const Core::LinAlg::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Get fiber directions
      void get_fiber_vecs(
          std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

      /// Setup of patient-specific materials
      void setup_aaa(Teuchos::ParameterList& params, const int eleGID) override { return; }

      /// Indicator for formulation
      void specify_formulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        anisoprinc = true;
        return;
      };

     private:
      /// my material parameters
      Mat::Elastic::PAR::CoupAnisoNeoHookeVarProp* params_;

      /// fiber direction
      Core::LinAlg::Matrix<3, 1> a_;
      /// structural tensors in voigt notation for anisotropy
      Core::LinAlg::Matrix<6, 1> structural_tensor_;
    };

  }  // namespace Elastic
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
