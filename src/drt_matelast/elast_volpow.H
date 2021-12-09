/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines for a volumetric power law.
The resultant pressure depends on the prefactor and pressure.
Should be combined with non-stiff volumetric law (eg. Ogden).
The input line should read
  MAT 1 ELAST_VolPow A 100 EXPON 5

\level 1


*----------------------------------------------------------------------*/
#ifndef ELAST_VOLPOW_H_
#define ELAST_VOLPOW_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for volumetric power law contribution \f$\Psi = \frac{a}{expon-1}
      /// J^(-expon+1) + aJ\f$
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_VolPow A 100 EXPON 5
      class VolPow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        VolPow(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// Prefactor
        double a_;
        /// Exponent
        double expon_;

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
        void SetA(const double a) { a_ = a; }
        void SetExpon(const double expon) { expon_ = expon; }
        //@}

      };  // class VolPow

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Volumetric Power-law material.
    ///
    /// Strain energy function is given by
    ///
    ///
    ///\f[
    ///  \Psi = \frac{a}{expon-1} J^(-expon+1) + aJ
    ///\f]
    ///
    /// results in the following pressure (stress)
    ///
    ///\f[
    /// p = -a (J^(-expon) -1)
    ///\f]
    ///
    /// \author birzle
    /// \date 11/15
    class VolPow : public Summand
    {
     public:
      /// constructor with given material parameters
      VolPow(MAT::ELASTIC::PAR::VolPow* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_volpow; }

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

      /// Add third derivative w.r.t. J
      void Add3rdVolDeriv(const LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3) override;


      /// @name Access methods
      //@{
      double A() const { return params_->a_; }
      double Expon() const { return params_->expon_; }
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
        isomod = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::VolPow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // VOLPOW_H_
