/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric one-term Ogden material.

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MATELAST_ISOOGDEN_HPP
#define BACI_MATELAST_ISOOGDEN_HPP

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
      class IsoOgden : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoOgden(Teuchos::RCP<MAT::PAR::Material> matdata);

        //! @name material parameters
        //! @{
        double mue_;    ///< shear modulus
        double alpha_;  ///< nonlinearity parameter
        //! @}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          dserror(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class IsoOgden
    }     // namespace PAR

    /*!
     * \brief Isochoric part of the one-term Ogden material, see Holzapfel [1] or Ogden [2]
     *
     * This is a hyperelastic isotropic material with compression-tension asymmetry expressed in
     * terms of the modified principal stretches. Amongst other applications, it can be used to
     * model human brain tissue [3]. In contrast to the fully incompressible formulation in [1],
     * here, the compressible formulation, i.e. the isochoric contribution of a decoupled
     * strain-energy function is implemented. Further, here the number of considered terms N is set
     * to one, s.t. the formulation reduces to a one-term Ogden model.
     *
     * The strain-energy function is hence given in terms of the modified principal stretches as:
     * \f[
     *   \Psi_{iso}=\frac{2\mu}{\alpha^2}\,(\bar{\lambda}_1^\alpha+\bar{\lambda}_2^\alpha+\bar{\lambda}_3^\alpha-3)
     * \f]
     *
     * References:
     * [1] G. A. Holzapfel, 'Nonlinear solid mechanics', Wiley, pp. 235-236, 2000.
     * [2] R. W. Ogden, 'Large deformation isotropic elasticity: on the correlation of theory and
     * experiment for compressible rubberlike solids', Proc. R. Soc. Lond. A, vol. 326, pp. 565-584,
     * 1972, doi: 10.1098/rspa.1972.0096.
     * [3] S. Budday et. al., 'Mechanical characterization of human brain tissue', Acta
     * Biomaterialia, vol. 48, pp. 319-340, 2017, doi: 10.1016/j.actbio.2016.10.036.
     */
    class IsoOgden : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoOgden(MAT::ELASTIC::PAR::IsoOgden* params);

      /// Provide the material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_isoogden; }

      /// Answer if coefficients with respect to modified principal stretches are provided
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

      /// Specify the formulation as isochoric in terms of modified principal invariants
      void SpecifyFormulation(bool& isoprinc, bool& isomod, bool& anisoprinc, bool& anisomod,
          bool& viscogeneral) override
      {
        isomod = true;
      };

     private:
      /// one-term Ogden material parameters
      MAT::ELASTIC::PAR::IsoOgden* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MATELAST_ISOOGDEN_H