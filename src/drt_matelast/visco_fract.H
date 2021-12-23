/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a viscous material contribution, calculated according to an
FSLS-model

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef VISCO_FRACT_H_
#define VISCO_FRACT_H_

#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for viscous contribution according the FSLS-model
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_Fract TAU 0.1 ALPHA 0.5 BETA 1
       */
      class Fract : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Fract(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double tau_;
        double alpha_;
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

        /// @name Write methods
        //@{
        void SetTau(const double tau) { tau_ = tau; }
        void SetAlpha(const double alpha) { alpha_ = alpha; }
        void SetBeta(const double beta) { beta_ = beta; }
        //@}

      };  // class Fract

    }  // namespace PAR


    /*!
     * @brief Material Viscofract
     *
     * This material offers a viscous and hyperelastic part. The model consists
     * of one spring in parallel to one sequential branch of a spring and a springpot.
     *
     * A springpot is between a spring and a dashpot. The parameter alpha regulates
     * how much damping is introduced.
     * Alpha=0, means the springpot is a spring
     * Alpha=1, means the springpot is a dashpot; this is equal to the GenMax Material
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] Adolfson and Enelund (2003): Fractional Derivative Visocelasticity at
     *          Large Deformations
     * </ul>
     */
    class Fract : public Summand
    {
     public:
      /// constructor with given material parameters
      Fract(MAT::ELASTIC::PAR::Fract* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_fract; }

      //@}

      /// Read material parameters
      void ReadMaterialParametersVisco(double& tau,  ///< relaxation parameter tau
          double& beta,                              ///< emphasis of viscous to elastic part
          double& alpha,  ///< fractional order derivative (just for visoc_fract)
          std::string&
              solve  //!< variant of the solution of the evolution integral (just for genmax)
          ) override;

      /// @name Access methods
      //@{
      double Tau() const { return params_->tau_; }
      double Alpha() const { return params_->alpha_; }
      double Beta() const { return params_->beta_; }

      //@}

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
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
        viscofract = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::Fract* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // FRACT_H_
