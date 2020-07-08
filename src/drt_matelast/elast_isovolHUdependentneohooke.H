/*----------------------------------------------------------------------*/
/*! \file
\brief  This file contains the routines required to calculate the isochoric and volumetric
contribution of the HU dependent neo hooke material. E.g used for inclusion of calcification in AAA
simulations

\level 2


*----------------------------------------------------------------------*/
#ifndef ELAST_ISOVOLHUDEPENDENTNEOHOOKE_H_
#define ELAST_ISOVOLHUDEPENDENTNEOHOOKE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// This file contains the routines required to calculate the
      /// isochoric and volumetric contribution of the HU dependent Neo
      /// Hooke material. The Hounsfield Units (HU) of e.g. CT data are
      /// used to determine the material stiffness for each finite
      /// element. The input line for the consideration of
      /// calcifications in AAA simulation reads:
      ///
      /// MAT 20 ELAST_IsoVolHUDependentNeoHooke ALPHA_MAX 8.929E6 CT_MIN 30.0 CT_MAX 600.0 NUE 0.49
      /// BETA -2.0
      ///
      /// An Neo Hooke material is chosen for the isochoric contribution
      /// with varying stiffness parameters according to the HU and the
      /// underlaying function given in Marc Hirschvogels Semesterarbeit
      /// (available in LNM library) The volumetric contribution is
      /// modeled by an Ogden-Simo-Miehe type SEF  \f$\Psi=\frac
      /// {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$  (see: Doll and
      /// Schweizerhof, On the Development of Volumetric Strain Energy
      /// Functions, Journal of Applied Mechanics, 2000). The Poisson
      /// ratio NUE has to be given in the material input line and the
      /// bulk modulus for each element is then calculated according to
      /// the isochoric stiffness.
      ///
      /// Strain energy function is given by
      /// \f[
      ///   \Psi = \alpha (\overline{I}_1 -3)
      ///          + \frac{\kappa}{\beta^2} \big(\beta log(J) + J^{-\beta} -1 \big)
      /// \f]
      ///
      /// with
      /// \f[
      ///   \kappa = \frac{2 \alpha}{1-2 \nu}
      /// \f]

      class IsoVolHUDependentNeoHooke : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoVolHUDependentNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~IsoVolHUDependentNeoHooke() { ; }

        /// @name material parameters
        //@{

        ///
        double alphamax_;
        ///
        double ctmin_;
        ///
        double ctmax_;
        ///
        double nue_;
        ///
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
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

      };  // class IsoVolHUDependentNeoHooke

    }  // namespace PAR


    class IsoVolHUDependentNeoHooke : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoVolHUDependentNeoHooke(MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params);

      ///@name Packing and Unpacking
      //@{

      virtual void PackSummand(DRT::PackBuffer& data) const;

      virtual void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position);

      //@}

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_isovolHUdependentneohooke;
      }

      //@}

      /// Setup of summand
      virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                ///< Gauss point
          const int eleGID                       ///< element GID
      );

      // Add derivatives with respect to modified invariants.
      virtual void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddmodPII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          int gp,                              ///< Gauss point
          const int eleGID                     ///< element GID
      );
      /// @name Access methods
      //@{
      double ALPHA_MAX() const { return params_->alphamax_; }
      double CT_MIN() const { return params_->ctmin_; }
      double CT_MAX() const { return params_->ctmax_; }
      double NUE() const { return params_->nue_; }
      double BETA() const { return params_->beta_; }
      //@}

      /// @name Write methods
      //@{
      void SetALPHA_MAX(const double alphamax) { params_->alphamax_ = alphamax; }
      void SetCT_MIN(const double ctmin) { params_->ctmin_ = ctmin; }
      void SetCT_MAX(const double ctmax) { params_->ctmax_ = ctmax; }
      void SetNUE(const double nue) { params_->nue_ = nue; }
      void SetBETA(const double beta) { params_->beta_ = beta; }
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
      MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params_;

      /// factor for material parameter
      double alpha_;
    };

  }  // namespace ELASTIC
}  // namespace MAT


#endif
