/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the routines required to calculate the volumetric contribution
of a general power-type material.

\level 2

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/
#ifndef ELAST_COUP3POW_H_
#define ELAST_COUP3POW_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for volumtric contribution of a general power material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_Coup3Pow C 1 D 1
      class Coup3Pow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Coup3Pow(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~Coup3Pow() { ; }

        /// @name material parameters
        //@{

        /// material parameters
        double c_;
        int d_;

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
        //@}

      };  // class Coup3Pow

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Volumtric general power material
    ///
    /// This is a summand of a variable order hyperelastic, volumetric
    /// material depending on the first invariant of the right Cauchy-Green tensor.
    ///
    /// Strain energy function is given by
    ///\f[
    ///   \Psi = C (I_3^(1/3)-1)^D.
    ///\f]
    ///
    ///
    /// \author birzle
    /// \date 07/17
    class Coup3Pow : public Summand
    {
     public:
      /// constructor with given material parameters
      Coup3Pow(MAT::ELASTIC::PAR::Coup3Pow* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coup3pow; }

      //@}

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1> glstrain,  ///< Green-Lagrange strain
          const int eleGID                      ///< element GID
      );

      virtual void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const int eleGID                    ///< element GID
      );

      /// @name Access methods
      //@{
      double C() const { return params_->c_; }
      int D() const { return params_->d_; }

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
      MAT::ELASTIC::PAR::Coup3Pow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // COUP3POW_H_
