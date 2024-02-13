/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric contribution of a anisotropic exponential fiber
material

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MATELAST_ISOANISOEXPO_HPP
#define BACI_MATELAST_ISOANISOEXPO_HPP

#include "baci_config.hpp"

#include "baci_mat_par_parameter.hpp"
#include "baci_matelast_summand.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of a anisotropic exponential fiber
       * material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_IsoAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE
       * 0
       */
      class IsoAnisoExpo : public MAT::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        IsoAnisoExpo(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double k1_;
        double k2_;
        /// angle between circumferential and fiber direction
        double gamma_;
        /// fiber params for the compressible case
        double k1comp_;
        double k2comp_;
        /// fiber initalization status
        int init_;
        /// adapt angle during remodeling
        int adapt_angle_;
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
      };  // class IsoAnisoExpo

    }  // namespace PAR

    /*!
     * @brief Isochoric anisotropic exponential fiber function, implemented for one possible fiber
     * family [1]
     *
     * This is a hyperelastic, anisotropic material of the most simple kind.
     *
     * Strain energy function is given by
     * \f[
     *    \Psi = \frac {k_1}{2 k_2} \left(e^{k_2 (\overline{IV}_{\boldsymbol C}-1)^2 }-1 \right).
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] G.A. Holzapfel, T.C. Gasser, R.W. Ogden: A new constitutive framework for arterial
     * wall mechanics and a comparative study of material models, J. of Elasticity 61 (2000) 1-48.
     * </ul>
     */
    class IsoAnisoExpo : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoAnisoExpo(MAT::ELASTIC::PAR::IsoAnisoExpo* params);

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
        return INPAR::MAT::mes_isoanisoexpo;
      }

      //@}

      /// Setup of summand
      void Setup(int numgp, INPUT::LineDefinition* linedef) override;

      /// Add anisotropic modified stresses
      void AddStressAnisoModified(
          const CORE::LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          const CORE::LINALG::Matrix<6, 1>& icg,  ///< inverse of right Cauchy Green Tensor
          CORE::LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          double I3,                              ///< third principal invariant
          int gp,                                 ///< Gauss point
          int eleGID,                             ///< element GID
          Teuchos::ParameterList& params          ///< Container for additional information
          ) override;

      /// retrieve coefficients of first, second and third derivative
      /// of summand with respect to anisotropic invariants
      virtual void GetDerivativesAniso(
          CORE::LINALG::Matrix<2, 1>& dPI_aniso,  ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<3, 1>&
              ddPII_aniso,  ///< second derivative with respect to invariants
          CORE::LINALG::Matrix<4, 1>&
              dddPIII_aniso,  ///< third derivative with respect to invariants
          double I4,          ///< fourth invariant
          int gp,             ///< Gauss point
          int eleGID);        ///< element GID

      /// Set fiber directions
      void SetFiberVecs(const double newgamma,       ///< new angle
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
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        anisomod = true;
        return;
      };

     protected:
      /// my material parameters
      MAT::ELASTIC::PAR::IsoAnisoExpo* params_;

      /// fiber direction
      CORE::LINALG::Matrix<3, 1> a_;
      /// structural tensors in voigt notation for anisotropy
      CORE::LINALG::Matrix<6, 1> A_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MATELAST_ISOANISOEXPO_H
