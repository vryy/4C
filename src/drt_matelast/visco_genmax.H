/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required to read the material parameters
and add the contribution of the viscous part, calculated according to a
SLS-Model.
Depending on the other stresses this viscous part is added, always using
the same tau and beta.
The input line should read
  MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE OST

\level 1

\maintainer Fabian Braeu

*----------------------------------------------------------------------*/
#ifndef VISCO_GENMAX_H_
#define VISCO_GENMAX_H_

#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for viscous contribution according the SLS-model
      ///
      /// <h3>Input line</h3>
      /// MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE OST
      /// MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE CONVOL
      class GenMax : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        GenMax(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~GenMax() { ; }

        /// @name material parameters
        //@{

        /// material parameters
        double tau_;
        double beta_;
        std::string solve_;

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
        void SetTau(const double tau) { tau_ = tau; }
        void SetBeta(const double beta) { beta_ = beta; }
        void SetSolve(const std::string solve) { solve_ = solve; }
        //@}

      };  // class GenMax

    }  // namespace PAR


    /*----------------------------------------------------------------------*/
    /*!
     *  Material Viscogenmax
     *
     *  This material offers a viscous and hyperelastic part. The model consists
     *  of one spring in parallel to one sequential branch of a spring and a dashpot.
     *
     *  The hyperelasic part is possibly any hyperelastic law of the Hyperelastic
     *  toolbox. The materials of the hyperelastic toolbox are composed
     *  of (Helmholtz free energy density) potentials.  Effectively, we want
     * \f[
     *   \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
     * \f]
     *
     *  in which the individual \f$\Psi_i\f$ is implemented as #MAT::ELASTIC::Summand.
     *
     *
     *  The stress tensor of the viscous part is calculated from the evolution equation
     * \f[
     *   \dot{\boldsymbol{Q}} = 1 / \tau \boldsymbol{Q} = \beta \dot{\boldsymbol{S}}
     * \f]
     *  with elastic stress tensor S. S contains isochoric and volumetric contribution,
     *  so Q is also isochoric and volumetric.
     *  The viscous elasicity tensor is calculated from \f$\boldsymbol{Q}\f$ with
     *  derivation to \f$\boldsymbol{C}\f$.
     *  The viscous effect can be applied to any part of the SEF (isotropic coupled,
     *  isotropic isochoric, isotropic volumetric, anisotropic).
     *
     *  <h3>References</h3>
     *  <ul>
     *  <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
     *  <li> [2] Bul-Brunon et.al., "Numerical identification method for the non-linear
     *          viscoelastic compressible behavior of soft tissue using uniaxial tensile
     *          tests and image registration - Application to rat lung parenchyma, 2013
     *  </ul>
     *
     *  \author bel, birzle
     *  \date 09/14
     */

    class GenMax : public Summand
    {
     public:
      /// constructor with given material parameters
      GenMax(MAT::ELASTIC::PAR::GenMax* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_genmax; }

      //@}

      /// Read material parameters
      virtual void ReadMaterialParametersVisco(double& tau,  ///< relaxation parameter tau
          double& beta,   ///< emphasis of viscous to elastic part
          double& alpha,  ///< fractional order derivative (just for visoc_fract)
          std::string&
              solve  //!< variant of the solution of the evolution integral (just for genmax)
      );

      /// @name Access methods
      //@{
      double Tau() const { return params_->tau_; }
      double Beta() const { return params_->beta_; }
      std::string Solve() const { return params_->solve_; }

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
        viscogenmax = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::GenMax* params_;
    };


  }  // namespace ELASTIC
}  // namespace MAT

#endif  // GENMAX_H_
