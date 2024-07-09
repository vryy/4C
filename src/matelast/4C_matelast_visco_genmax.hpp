/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a viscous material contribution, calculated according to an SLS-model

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_VISCO_GENMAX_HPP
#define FOUR_C_MATELAST_VISCO_GENMAX_HPP

#include "4C_config.hpp"

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
       * @brief material parameters for viscous contribution according the SLS-model
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE OST
       * MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE CONVOL
       */
      class GenMax : public Core::Mat::PAR::Parameter
      {
       public:
        /// standard constructor
        GenMax(const Core::Mat::PAR::Parameter::Data& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double tau_;
        double beta_;
        std::string solve_;

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
      };  // class GenMax

    }  // namespace PAR

    /*!
     * @brief Material Viscogenmax
     *
     * This material offers a viscous and hyperelastic part. The model consists
     * of one spring in parallel to one sequential branch of a spring and a dashpot.
     *
     * The hyperelasic part is possibly any hyperelastic law of the Hyperelastic
     * toolbox. The materials of the hyperelastic toolbox are composed
     * of (Helmholtz free energy density) potentials.  Effectively, we want
     * \f[
     *  \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
     * \f]
     *
     * in which the individual \f$\Psi_i\f$ is implemented as #Mat::Elastic::Summand.
     *
     *
     * The stress tensor of the viscous part is calculated from the evolution equation
     * \f[
     *  \dot{\boldsymbol{Q}} = 1 / \tau \boldsymbol{Q} = \beta \dot{\boldsymbol{S}}
     * \f]
     * with elastic stress tensor S. S contains isochoric and volumetric contribution,
     * so Q is also isochoric and volumetric.
     * The viscous elasicity tensor is calculated from \f$\boldsymbol{Q}\f$ with
     * derivation to \f$\boldsymbol{C}\f$.
     * The viscous effect can be applied to any part of the SEF (isotropic coupled,
     * isotropic isochoric, isotropic volumetric, anisotropic).
     *
     * <h3>References</h3>
     * <ul>
     * <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     * <li> [2] Bul-Brunon et.al., "Numerical identification method for the non-linear
     *         viscoelastic compressible behavior of soft tissue using uniaxial tensile
     *         tests and image registration - Application to rat lung parenchyma, 2013
     * </ul>
     */

    class GenMax : public Summand
    {
     public:
      /// constructor with given material parameters
      GenMax(Mat::Elastic::PAR::GenMax* params);

      /// @name Access material constants
      //@{

      /// material type
      Core::Materials::MaterialType material_type() const override
      {
        return Core::Materials::mes_genmax;
      }

      //@}

      /// Read material parameters
      void read_material_parameters_visco(double& tau,  ///< relaxation parameter tau
          double& beta,                                 ///< emphasis of viscous to elastic part
          double& alpha,  ///< fractional order derivative (just for visoc_fract)
          std::string&
              solve  //!< variant of the solution of the evolution integral (just for genmax)
          ) override;

      /// Indicator for formulation
      void specify_formulation(
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
      void specify_visco_formulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
          ) override
      {
        viscogenmax = true;
        return;
      };


     private:
      /// my material parameters
      Mat::Elastic::PAR::GenMax* params_;
    };


  }  // namespace Elastic
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
