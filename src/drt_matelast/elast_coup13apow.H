/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required to calculate the coupled
power-type material in invariant 1 multiplicative invariant 3^(-a).
The input line should read
  MAT 1 ELAST_Coup13aPow C 1 D 1 A 0.1

\level 1

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/
#ifndef ELAST_COUP13APOW_H_
#define ELAST_COUP13APOW_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for coupled contribution of a multiplicative coupled I1 and I3 power
      /// material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_Coup13aPow C 1 D 1 A 1
      class Coup13aPow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Coup13aPow(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~Coup13aPow() { ; }

        /// @name material parameters
        //@{

        /// material parameters
        double c_;  // constant
        int d_;     // exponent of all
        double a_;  // negative exponent of I3

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
        void SetC(const double c) { c_ = c; }
        void SetD(const int d) { d_ = d; }
        void SetA(const double a) { a_ = a; }
        //@}

      };  // class Coup13aPow

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Coupled general power material
    ///
    /// This is a summand of a variable order hyperelastic, coupled
    /// material depending on the first and third invariant of the right
    /// Cauchy-Green tensor. The third invariant has a additional negative exponent.
    ///
    /// Strain energy function is given by
    ///\f[
    ///   \Psi = c (I_{\boldsymbol{C}}*(III_{\boldsymbol{C}}^(-a))-3)^d.
    ///\f]
    ///
    ///
    /// \author birzle
    /// \date 11/16
    class Coup13aPow : public Summand
    {
     public:
      /// constructor with given material parameters
      Coup13aPow(MAT::ELASTIC::PAR::Coup13aPow* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coup13apow; }

      //@}

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>&
              glstrain,  ///< Green-Lagrange strain in strain like Voigt notation
          int gp,        ///< Gauss point
          int eleGID     ///< element GID
      );

      virtual void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,                             ///< Gauss point
          int eleGID                          ///< element GID
      );

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      virtual void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
      {
        dserror("not implemented");
      }

      /// @name Access methods
      //@{
      double C() const { return params_->c_; }
      int D() const { return params_->d_; }
      double A() const { return params_->a_; }

      //@}

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
      )
      {
        isoprinc = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::Coup13aPow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // COUP13APOW_H_
