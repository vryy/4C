/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the isochoric contribution of a viscous rate dependent material
law, modified from Pioletti, 1997

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef VISCO_ISORATEDEP_H_
#define VISCO_ISORATEDEP_H_

#include "matpar_parameter.H"
#include "elast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for isochoric contribution of a frequency independent
       * viscoelastic material
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_IsoRateDep N 1
       */
      class IsoRateDep : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoRateDep(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double n_;

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
      };  // class IsoRateDep

    }  // namespace PAR

    /*!
     * @brief Isochoric general power material
     *
     * This is the isochoric part of a viscohyperelastic, isotropic
     * material depending on the modified invariants of the rate of the right
     * Cauchy-Green tensor.
     *
     * Strain energy function is given by
     * \f[
     *   \Psi = n \overline{J}_2 (\overline{I}_1 -3).
     * \f]
     * (n = \eta)
     */
    class IsoRateDep : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoRateDep(MAT::ELASTIC::PAR::IsoRateDep* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_isoratedep; }

      //@}

      /// Add modified coeffiencts.
      void AddCoefficientsViscoModified(
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          LINALG::Matrix<8, 1>& modmu,   ///< necassary coefficients for piola-kirchhoff-stress
          LINALG::Matrix<33, 1>& modxi,  ///< necassary coefficients for viscosity tensor
          LINALG::Matrix<7, 1>& modrateinv, Teuchos::ParameterList& params, int gp,
          int eleGID) override;

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
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      void SpecifyViscoFormulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
          ) override
      {
        isovisco = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::IsoRateDep* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ISORATEDEP_H_
