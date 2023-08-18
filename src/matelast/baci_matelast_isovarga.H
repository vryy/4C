/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for Varga's material

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MATELAST_ISOVARGA_H
#define BACI_MATELAST_ISOVARGA_H

#include "baci_mat_par_parameter.H"
#include "baci_matelast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters of Varga's material
       *
       * <h3>Input line</h3>
       * MAT 1 ELAST_IsoVarga MUE 1.0 BETA 1.0
       */
      class IsoVarga : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoVarga(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double mue_;
        /// 'Anti-modulus'
        double beta_;

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

      };  // class IsoVarga

    }  // namespace PAR

    /*!
     * \brief Isochoric Varga's material according to [1], [2]
     *
     * This is a compressible, hyperelastic, isotropic material
     * of the most simple kind.
     *
     * The material strain energy density function is
     * \f[
     * \Psi = \underbrace{(2\mu-\beta)}_{\displaystyle\alpha} \Big(\bar{\lambda}_1 +
     * \bar{\lambda}_2 + \bar{\lambda}_3 - 3\Big)
     *      + \beta \Big(\frac{1}{\bar{\lambda}_1} + \frac{1}{\bar{\lambda}_2} +
     *      \frac{1}{\bar{\lambda}_3} - 3\Big)
     * \f]
     * which was taken from [1] Eq (6.129) and [2] Eq (1.3).
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     * <li> [2] JM Hill and DJ Arrigo, "New families of exact solutions for
     *          finitely deformed incompressible elastic materials",
     *          IMA J Appl Math, 54:109-123, 1995.
     * </ul>
     */
    class IsoVarga : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoVarga(MAT::ELASTIC::PAR::IsoVarga* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_isovarga; }

      /// add shear modulus equivalent
      void AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
          double& shearmod                  ///< variable to add upon
      ) const override;

      //@}

      /// Answer if coefficients with respect to principal stretches are provided
      bool HaveCoefficientsStretchesModified() override { return true; }

      /// Add coefficients with respect to modified principal stretches (or zeros)
      void AddCoefficientsStretchesModified(
          CORE::LINALG::Matrix<3, 1>&
              modgamma,  ///< [\bar{\gamma}_1, \bar{\gamma}_2, \bar{\gamma}_3]
          CORE::LINALG::Matrix<6, 1>&
              moddelta,  ///< [\bar{\delta}_11, \bar{\delta}_22, \bar{\delta}_33,
                         ///< \bar{\delta}_12,\bar{\delta}_23, \bar{\delta}_31]
          const CORE::LINALG::Matrix<3, 1>&
              modstr  ///< modified principal stretches, [\bar{\lambda}_1,
                      ///< \bar{\lambda}_2, \bar{\lambda}_3]
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
        return;
      };

     private:
      /// Varga material parameters
      MAT::ELASTIC::PAR::IsoVarga* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // MATELAST_ISOVARGA_H
