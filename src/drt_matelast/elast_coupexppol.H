/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for a strain energy function,
which is exponential according to Weickenmeier_2014 and contains a linear
(in I_1 and J) polynom in the exponent.
The input line should read
  MAT 1 ELAST_CoupExpPol A 600. B 2. C 5.

\level 1


*----------------------------------------------------------------------*/
#ifndef ELAST_EXPPOL_H_
#define ELAST_EXPPOL_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for compressible soft tissue material
      ///
      /// <h3>Input line</h3>
      ////  MAT 1 ELAST_CoupExpPol A 600. B 2. C 5.

      class CoupExpPol : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        CoupExpPol(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{
        double a_;
        /// constant for linear part of I_1
        double b_;
        /// constant for linear part of J
        double c_;
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
        void SetB(const double b) { b_ = b; }
        void SetC(const double c) { c_ = c; }
        //@}

      };  // class CoupExpPol

    }  // namespace PAR

    /* ----------------------------------------------------------------*/
    /// coupled compressible, isotropic material
    /// according to [1] with linear summands in exponent
    ///
    /// This is a compressible, hyperelastic, isotropic material.
    ///
    /// The material strain energy density function is
    ///\f[
    /// \Psi = a \exp[ b(I_1 - 3) - (2b + c)ln{J} + c(J-1) ] - a
    ///\f]
    /// More details at #AddCoefficientsPrincipal()
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] Weickenmeier, Jabareen "Elastic-viscoplastic modeling of soft biological
    ///           tissues using a mixed finite element formulation based on the relative
    ///           deformation gradient", 2014
    /// </ul>
    ///
    /// \author abirzle
    /// \date 07/14
    class CoupExpPol : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupExpPol(MAT::ELASTIC::PAR::CoupExpPol* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_coupexppol; }

      // add strain energy
      void AddStrainEnergy(double& psi,       ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                ///< Gauss point
          int eleGID                             ///< element GID
          ) override;


      void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,                             ///< Gauss point
          int eleGID                          ///< element GID
          ) override;
      /// @name Access methods
      //@{
      double A() const { return params_->a_; }
      double B() const { return params_->b_; }
      double C() const { return params_->c_; }

      //@}

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        isoprinc = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupExpPol* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // EXPPOl_H_
