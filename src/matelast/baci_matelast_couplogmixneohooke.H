/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a mixed-logarithmic neo-Hooke material according to example 1 in
"Stupkiewicz. An ALE formulation for implicit time integration of quasi-steady state wear problems.
Comput. Methods Appl. Mech. Engrg. 260, 2013"

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MATELAST_COUPLOGMIXNEOHOOKE_H
#define BACI_MATELAST_COUPLOGMIXNEOHOOKE_H

#include "baci_matelast_summand.H"
#include "baci_mat_par_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for mixed logarithmic neo-Hooke material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_CoupLogMixNeoHooke MODE YN C1 1.0 C2 0.3 (C1 = Young, c2 = nue)
       * MAT 1 ELAST_CoupLogMixNeoHooke MODE Lame C1 1.0 C2 1.0
       */
      class CoupLogMixNeoHooke : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupLogMixNeoHooke(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double mue_;
        /// Lame's constant
        double lambda_;

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
      };  // class CoupLogNeoHooke

    }  // namespace PAR

    /*!
     * @brief Logarithmic mixed neo-Hooke material:
     *
     * This is a compressible, hyperelastic, isotropic material.
     *
     * The material strain energy density function is
     * \f[
     *    \Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3) - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
     * + \frac{\lambda}{2} \big( \sqrt{I\!I\!I_{\boldsymbol{C}}} -1 \big)^2 \f]
     */
    class CoupLogMixNeoHooke : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupLogMixNeoHooke(MAT::ELASTIC::PAR::CoupLogMixNeoHooke* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_couplogmixneohooke;
      }

      /// add shear modulus equivalent
      void AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
          double& shearmod                  ///< variable to add upon
      ) const override;

      //@}

      /*!
       * @brief Main material call to determine 2nd PK stress and material constitutive tensor
       *
       * The material strain energy density function is shown above.
       *
       * Here is \f$I_{\boldsymbol{C}}\f$ the first principal invariant
       * of the right Cauchy--Green strain tensor \f$\boldsymbol{C}\f$
       * and \f$I\!I\!I_{\boldsymbol{C}}\f$ its third.
       * The isochoric part is proportional to \f$\mu\f$,
       * whereas the \f$\lambda\f$-proportional part constitutes the volumetric contribution.
       *
       * The parameters are the shear modulus. In #params_->parmode_==0 it is set directly,
       * in #params_->parmode_==1 it is computed by
       * \f[
       * \mu = \frac{E}{2(1+\nu)}
       * \f]
       * and Lame's coefficient, again for #params_->parmode_==1, otherwise set by user,
       * \f[
       * \lambda = \left\{\begin{array}{ll}
       *      \frac{E \nu}{(1+\nu) (1-2\nu)} & \text{if $\nu \neq 1/2$}
       *   \\ 0                              & \text{else}
       * \end{array}\right.
       * \f]
       */

      // add strain energy
      void AddStrainEnergy(double& psi,  ///< strain energy function
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<3, 1>&
              modinv,  ///< modified invariants of right Cauchy-Green tensor
          const CORE::LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                      ///< Gauss point
          int eleGID                                   ///< element GID
          ) override;

      void AddDerivativesPrincipal(
          CORE::LINALG::Matrix<3, 1>& dPI,    ///< first derivative with respect to invariants
          CORE::LINALG::Matrix<6, 1>& ddPII,  ///< second derivative with respect to invariants
          const CORE::LINALG::Matrix<3, 1>&
              prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,     ///< Gauss point
          int eleGID  ///< element GID
          ) override;

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4) override
      {
        dserror("not implemented");
      }

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        isoprinc = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupLogMixNeoHooke* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // MATELAST_COUPLOGMIXNEOHOOKE_H
