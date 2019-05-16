/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite
element analysis", Cambridge, 1997.
The input line should read
  MAT 1 ELAST_VolSussmanBathe KAPPA 100

\level 1

<pre>
\maintainer Fabian Braeu
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_VOLSUSSMANBATHE_H_
#define ELAST_VOLSUSSMANBATHE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for volumetric contribution \f$\Psi=\kappa(J-1)^2\f$
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_VolSussmanBathe KAPPA 10000
      class VolSussmanBathe : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        VolSussmanBathe(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~VolSussmanBathe() { ; }

        /// @name material parameters
        //@{

        /// Dilatation modulus
        double kappa_;

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
        void SetKappa(const double kappa) { kappa_ = kappa; }
        //@}

      };  // class VolSussmanBathe

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Volumetric SussmanBathe material according to [1].
    ///
    /// Strain energy function is given by
    ///
    ///
    ///\f[
    ///   \Psi = \frac \kappa 2 (J-1)^2
    ///\f]
    ///
    ///
    /// [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
    ///     Journal of Applied Mechanics, 2000
    /// \author rausch/tk
    /// \date 04/09
    class VolSussmanBathe : public Summand
    {
     public:
      /// constructor with given material parameters
      VolSussmanBathe(MAT::ELASTIC::PAR::VolSussmanBathe* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_volsussmanbathe;
      }

      //@}

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          const int eleGID                       ///< element GID
      );


      // Add derivatives with respect to modified invariants.
      virtual void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddmodPII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const int eleGID                     ///< element GID
      );

      /// Add third derivative w.r.t. J
      virtual void Add3rdVolDeriv(const LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3);


      /// @name Access methods
      //@{
      double Kappa() const { return params_->kappa_; }

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
        isomod = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::VolSussmanBathe* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // VOLSUSSMANBATHE_H_
