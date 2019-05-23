/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required to calculate the isochoric contribution
of a Yeoh-type material.
The input line should read
  MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200

\level 1

\maintainer Fabian Braeu

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOYEOH_H_
#define ELAST_ISOYEOH_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for isochoric contribution of a Yeoh material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200
      class IsoYeoh : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoYeoh(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~IsoYeoh() { ; }

        /// @name material parameters
        //@{

        /// Shear modulus
        double c1_;
        double c2_;
        double c3_;

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
        void SetC1(const double c1) { c1_ = c1; }
        void SetC2(const double c2) { c2_ = c2; }
        void SetC3(const double c3) { c3_ = c3; }
        //@}

      };  // class IsoYeoh

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Isochoric Yeoh material according to [1]
    ///
    /// This is the isochoric part of a third order hyperelastic, isotropic
    /// material depending on the modified trace of the right Cauchy-Green tensor.
    ///
    /// Strain energy function is given by
    ///\f[
    ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2 (\overline{I}_{\boldsymbol{C}}-3)^2 +
    ///          C3 (\overline{I}_{\boldsymbol{C}}-3)^3.
    ///\f]
    ///
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
    /// </ul>
    ///
    /// \author rausch/tk
    /// \date 04/09
    class IsoYeoh : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoYeoh(MAT::ELASTIC::PAR::IsoYeoh* params);

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoyeoh; }

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


      /// @name Access methods
      //@{
      double C1() const { return params_->c1_; }
      double C2() const { return params_->c2_; }
      double C3() const { return params_->c3_; }

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
        isomod = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::IsoYeoh* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ISOYEOH_H_
