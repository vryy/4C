/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric contribution of an isotropic general power-type
material in terms of the second Cauchy-Green invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef ELAST_ISO2POW_H_
#define ELAST_ISO2POW_H_

#include "mat_par_parameter.H"
#include "matelast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of a general power material
       *
       *  <h3>Input line</h3>
       *  MAT 1 ELAST_Iso2Pow C 1 D 1
       */
      class Iso2Pow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Iso2Pow(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double c_;
        int d_;

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
      };  // class Iso2Pow

    }  // namespace PAR

    /*!
     * @brief Isochoric general power material
     *
     * This is the isochoric part of a variable order hyperelastic, isotropic
     * material depending on the seccond invariant of the right Cauchy-Green tensor.
     *
     * Strain energy function is given by
     * \f[
     *    \Psi = C (\overline{II}_{\boldsymbol{C}}-3)^D.
     * \f]
     */
    class Iso2Pow : public Summand
    {
     public:
      /// constructor with given material parameters
      Iso2Pow(MAT::ELASTIC::PAR::Iso2Pow* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_iso2pow; }

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,       ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                ///< Gauss point
          const int eleGID                       ///< element GID
          ) override;

      // Add derivatives with respect to modified invariants.
      void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddmodPII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          int gp,                              ///< Gauss point
          const int eleGID                     ///< element GID
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
        isomod = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::Iso2Pow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ISO2POW_H_
