/*----------------------------------------------------------------------*/
/*! \file
\brief
Definition of classes for the isochoric contribution of a Yeoh-type material

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef ELAST_ISOYEOH_H_
#define ELAST_ISOYEOH_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of a Yeoh material
       *
       *  <h3>Input line</h3>
       *  MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200
       */
      class IsoYeoh : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoYeoh(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Shear modulus
        double c1_;
        double c2_;
        double c3_;

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
      };  // class IsoYeoh
    }     // namespace PAR

    /*!
     * @brief Isochoric Yeoh material according to [1]
     *
     * This is the isochoric part of a third order hyperelastic, isotropic
     * material depending on the modified trace of the right Cauchy-Green tensor.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2 (\overline{I}_{\boldsymbol{C}}-3)^2 +
     *          C3 (\overline{I}_{\boldsymbol{C}}-3)^3.
     * \f]
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     * </ul>
     */
    class IsoYeoh : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoYeoh(MAT::ELASTIC::PAR::IsoYeoh* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_isoyeoh; }

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,       ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                ///< Gauss point
          int eleGID                             ///< element GID
          ) override;

      // Add derivatives with respect to modified invariants.
      void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddmodPII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          int gp,                              ///< Gauss point
          int eleGID                           ///< element GID
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
        isomod = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::IsoYeoh* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ISOYEOH_H_
