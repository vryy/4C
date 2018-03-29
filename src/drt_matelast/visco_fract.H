/*----------------------------------------------------------------------*/
/*!
\file visco_fract.H
\brief
This file contains the routines required to read the material parameters
and add the contribution of the viscous part, calculated according to a
FSLS-Model.
Depending on the other stresses this viscous part is added, always using
the same tau, alpha and beta.

\level 2

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>

*----------------------------------------------------------------------*/
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

    /*----------------------------------------------------------------------*/
    /// material parameters for viscous contribution according the FSLS-model
    ///
    /// <h3>Input line</h3>
    /// MAT 1 VISCO_Fract TAU 0.1 ALPHA 0.5 BETA 1

    class Fract : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      Fract(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Fract() { ; }

      /// @name material parameters
      //@{

      /// material parameters
      double tau_;
      double alpha_;
      double beta_;

      //@}

      /// create material instance of matching type with my parameters
      ///
      /// Remark:                                                  AN 03/2014
      /// Function has to be implemented, since it is a pure virtual
      /// function in the MAT::PAR::Material class.
      /// Function stays empty, since the call within the material factory is done only
      /// respective the elasthyper class
      ///
      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

      /// @name Write methods
      //@{
      void SetTau(const double tau) { tau_ = tau; }
      void SetAlpha(const double alpha) { alpha_ = alpha; }
      void SetBeta(const double beta) { beta_ = beta; }
      //@}

    };  // class Fract

  }  // namespace PAR


/*----------------------------------------------------------------------*/
/*!
 *  Material Viscofract
 *
 *  This material offers a viscous and hyperelastic part. The model consists
 *  of one spring in parallel to one sequential branch of a spring and a springpot.
 *
 *  A springpot is between a spring and a dashpot. The parameter alpha regulates
 *  how much damping is introduced.
 *  Alpha=0, means the springpot is a spring
 *  Alpha=1, means the springpot is a dashpot; this is equal to the GenMax Material
 *
 *  <h3>References</h3>
 *  <ul>
 *  <li> [1] Adolfson and Enelund (2003): Fractional Derivative Visocelasticity at
 *           Large Deformations
 *  </ul>
 *
 *  \author birzle
 *  \date 12/17
 */

  class Fract : public Summand
  {

  public:

    /// constructor with given material parameters
    Fract(MAT::ELASTIC::PAR::Fract* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_fract; }

    //@}

    /// Read material parameters
    virtual void ReadMaterialParametersVisco(
      double& tau, ///< relaxation parameter tau
      double& beta, ///< emphasis of viscous to elastic part
      double& alpha, ///< fractional order derivative (just for visoc_fract)
      std::string& solve //!< variant of the solution of the evolution integral (just for genmax)
      );

    /// @name Access methods
    //@{
    double Tau() const { return params_->tau_; }
    double Alpha() const { return params_->alpha_; }
    double Beta() const { return params_->beta_; }

    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< general indicator, if one viscoelastic formulation is used
      )
    {
      viscogeneral = true;
      return;
    };

    /// Indicator for the chosen viscoelastic formulations
    virtual void SpecifyViscoFormulation(
        bool& isovisco, ///< global indicator for isotropic, splitted and viscous formulation
        bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
        bool& viscogeneralizedgenmax, ///< global indicator for viscoelastic contribution according to the generalized Maxwell Model
        bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
    )
    {
      viscofract = true;
      return;
    };


  private:

    /// my material parameters
    MAT::ELASTIC::PAR::Fract* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // FRACT_H_
