/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an isochoric coupled viscous material with pseudo-potential
representing the collagen and elastin matrix surrounding the myocardial fiber (chappelle12)

\level 2
*/
/*---------------------------------------------------------------------*/
#ifndef BACI_MATELAST_VISCO_COUPMYOCARD_H
#define BACI_MATELAST_VISCO_COUPMYOCARD_H

#include "baci_mat_par_parameter.H"
#include "baci_matelast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for viscous part of myocardial matrix
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_CoupMyocard N 1
       */
      class CoupMyocard : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupMyocard(const Teuchos::RCP<MAT::PAR::Material>& matdata);

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
      };  // class CoupMyocard
    }     // namespace PAR

    /*!
     * @brief Isochoric coupled viscous material with pseudo-potential
     *
     * Strain energy function is given by
     * \f[
     *   \Psi_v = \eta/2 tr(\dot{E}^2) = \eta/8 tr(\dot{C}^2).
     * \f]
     *
     * Viscous second Piola-Kirchhoff stress
     * \f[
     *   S_v =  2 \frac{\partial \Psi_v}{\partial \dot{C}} = \eta/2 \dot{C}.
     * \f]
     *
     * Viscous constitutive tensor
     * \f[
     *   C_v =  4 \frac{\partial^2 W_v}{\partial \dot{C} \partial \dot{C}} = \eta I^\#,
     * \f]
     *
     * with
     *
     * \f[
     *   I^\#_{ijkl} = \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
     * \f]
     */
    class CoupMyocard : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupMyocard(MAT::ELASTIC::PAR::CoupMyocard* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_coupmyocard; }

      //@}

      /// Add modified coeffiencts.
      void AddCoefficientsViscoPrincipal(
          const CORE::LINALG::Matrix<3, 1>& prinv,  ///< invariants of right Cauchy-Green tensor
          CORE::LINALG::Matrix<8, 1>& mu,   ///< necassary coefficients for piola-kirchhoff-stress
          CORE::LINALG::Matrix<33, 1>& xi,  ///< necassary coefficients for viscosity tensor
          CORE::LINALG::Matrix<7, 1>& rateinv, Teuchos::ParameterList& params, int gp,
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
        isoprinc = true;
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
      MAT::ELASTIC::PAR::CoupMyocard* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // MATELAST_VISCO_COUPMYOCARD_H
