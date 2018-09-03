/*----------------------------------------------------------------------*/
/*!
\file visco_coupmyocard.H

\brief
Isochoric coupled viscous material with pseudo-potential representing the collagen and
elastin matrix surrounding the myocardial fiber (chappelle12)
The input line should read
  MAT 1 VISCO_CoupMyocard N 1

\level 2

<pre>
\maintainer Martin Pfaller
            pfaller@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef VISCO_COUPMYOCARD_H_
#define VISCO_COUPMYOCARD_H_

#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for viscous part of myocardial matrix
      ///
      /// <h3>Input line</h3>
      /// MAT 1 VISCO_CoupMyocard N 1
      class CoupMyocard : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupMyocard(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~CoupMyocard() { ; }

        /// @name material parameters
        //@{

        /// material parameters
        double n_;

        //@}

        /// create material instance of matching type with my parameters
        ///
        /// Remark:                                                  AN 03/2014
        /// Function has to be implemented, since it is a pure virtual
        /// function in the MAT::PAR::Material class.
        /// Function stays empty, since the call within the material factory is done only
        /// respective the elasthyper class
        ///
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{
        void SetN(const double n) { n_ = n; }
        //@}

      };  // class CoupMyocard

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /*!
     *  Isochoric coupled viscous material with pseudo-potential
     *
     *  Strain energy function is given by
     * \f[
     *    \Psi_v = \eta/2 tr(\dot{E}^2) = \eta/8 tr(\dot{C}^2).
     * \f]
     *
     *  Viscous second Piola-Kirchhoff stress
     * \f[
     *    S_v =  2 \frac{\partial \Psi_v}{\partial \dot{C}} = \eta/2 \dot{C}.
     * \f]
     *
     *  Viscous constitutive tensor
     * \f[
     *    C_v =  4 \frac{\partial^2 W_v}{\partial \dot{C} \partial \dot{C}} = \eta I^\#,
     * \f]
     *
     *  with
     *
     * \f[
     *    I^\#_{ijkl} = \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
     * \f]
     *
     *  \author pfaller
     *  \date 04/15
     */
    class CoupMyocard : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupMyocard(MAT::ELASTIC::PAR::CoupMyocard* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupmyocard; }

      //@}

      /// Add modified coeffiencts.
      virtual void AddCoefficientsViscoPrincipal(
          const LINALG::Matrix<3, 1>& prinv,  ///< invariants of right Cauchy-Green tensor
          LINALG::Matrix<8, 1>& mu,           ///< necassary coefficients for piola-kirchhoff-stress
          LINALG::Matrix<33, 1>& xi,          ///< necassary coefficients for viscosity tensor
          LINALG::Matrix<7, 1>& rateinv, Teuchos::ParameterList& params, const int eleGID);

      /// @name Access methods
      //@{
      double N() const { return params_->n_; }

      //@}

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
      )
      {
        isoprinc = true;
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      virtual void SpecifyViscoFormulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
      )
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

#endif  // COUPMYOCARD_H_
