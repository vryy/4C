/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a pow-like anisotropic material

\level 3
*/
/*---------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOPOW_H_
#define ELAST_COUPANISOPOW_H_

#include "matelast_summand.H"
#include "mat_par_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for anisochoric contribution of a pow like material with one
       * fiber direction
       *
       * <h3>Input line</h3>
       * MAT 1 CoupAnisoPow C 1.0 D 2.0 [ GAMMA 35.0 INIT 0 ADAPT_ANGLE 0]
       */
      class CoupAnisoPow : public MAT::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        CoupAnisoPow(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// stiffness factor
        double k_;
        /// exponential factor for I4
        double d1_;
        /// exponential factor for ((I4)^{d1}-1)
        double d2_;
        /// fiber number (1,2,3,...) used later as FIBER1,FIBER2,FIBER3,...
        int fibernumber_;
        /// Deformation threshold where fibers are active.
        double activethres_;
        /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
        double gamma_;
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
      };  // class CoupAnisoPow

    }  // namespace PAR

    /*!
     * @brief Coupled anisotropic pow-like fiber function, implemented for one possible fiber family
     * as in [1]
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = K \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2}.
     * \f]
     *
     * The corresponding derivatives are
     * \f[
     *   \frac{d\Psi}{d IV_{\boldsymbol C}} = K\ D2\ D1\ (IV_{\boldsymbol C})^{D1-1}
     *   \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2-1}.
     * \f]
     * and
     * \f[
     *   \frac{d^2\Psi}{d^2 IV_{\boldsymbol C}} = K\ D2\ (D2-1)\ \left( D1\ (IV_{\boldsymbol
     *   C})^{D1-1} \right)^2 \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2-2}\ +\ K\ D2\ D1
     *   (D1-1)\ (IV_{\boldsymbol C})^{D1-2}\ \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2-1}
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, Nonlinear solid mechanics 2004
     * </ul>
     */
    class CoupAnisoPow : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupAnisoPow(MAT::ELASTIC::PAR::CoupAnisoPow* params);

      ///@name Packing and Unpacking
      //@{
      void PackSummand(DRT::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;
      //@}

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_coupanisoneohooke;
      }

      /// Setup of summand
      void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

      /// Add anisotropic principal stresses
      void AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          LINALG::Matrix<6, 6>& cmat,                                ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,                              ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// Set fiber directions
      void SetFiberVecs(const double newgamma,  ///< new angle
          const LINALG::Matrix<3, 3>& locsys,   ///< local coordinate system
          const LINALG::Matrix<3, 3>& defgrd    ///< deformation gradient
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
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoPow* params_;

      /// fiber direction
      LINALG::Matrix<3, 1> a_;

      /// structural tensors in Voigt notation for anisotropy
      LINALG::Matrix<6, 1> A_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // COUPANISOPOW_H_
