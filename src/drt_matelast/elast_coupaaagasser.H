/*----------------------------------------------------------------------*/
/*!
\file elast_coupaaagasser.H
\brief

This file contains the routines required to calculate the isochoric contribution
of the aaagasser material and the corresponding volumetric contribution.

MAT 20 ELAST_CoupAAAGasser CLUM 2.62E3 CMED 1.98E3 CABLUM 1.73E3 NUE 0.49 BETA -2.0

<pre>
maintainer: Andreas Maier
            a.maier@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPAAAGASSER_H_
#define ELAST_COUPAAAGASSER_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// This material provides the isochoric and(!) the volumetric
    /// contribution for a AAA thromubs. The input line should read
    /// like this:
    ///     MAT 21 ELAST_CoupAAAGasser CLUM 2.62E3 CMED 1.98E3
    ///     CABLUM 1.73E3 NUE 0.49 BETA -2.0
    ///
    /// An Ogden type material
    /// \f$\Psi=c\underset{i=1}{\overset{3}{\sum}}(\lambda_{i}^{4}-1)\f$ is chosen for the isochoric
    /// contribution with spatially varying stiffness parameters
    /// governed by luminal thrombus layer stiffness CLUM, the medial
    /// layer stiffness CMED and the abluminal thrombus layer
    /// stiffness CABLUM (see: Gasser
    /// et al., Failure properties of intraluminal thrombus in
    /// abdominal aortic aneurysm under static and pulsating
    /// mechanical loads, J Vasc Surg, 2008). The volumetric
    /// contribution is modeled by an Ogden-Simo_Miehe type SEF:
    /// \f$\Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$ (see  Doll and Schweizerhof, On the Development of Volumetric Strain
    /// Energy Functions, Journal of Applied Mechanics, 2000).


    class CoupAAAGasser : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      CoupAAAGasser(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupAAAGasser() { ; }

      /// @name material parameters
      //@{

      /// stiffness parameter (luminal)
      const double clum_;
      /// stiffness parameter (medial)
      const double cmed_;
      /// stiffness parameter (abluminal)
      const double cablum_;
      /// poisson ration
      const double nue_;
      /// numerical parameter (should be -2.0 (Doll and Schwizerhof)
      /// or 9.0(Holzapfel and Ogden))
      const double beta_;

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

    };  // class CoupAAAGasser

  }  // namespace PAR


  class CoupAAAGasser : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupAAAGasser(MAT::ELASTIC::PAR::CoupAAAGasser* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupaaagasser; }

    //@}

    /// Setup of patient-specific stuff
    virtual void SetupAAA(Teuchos::ParameterList& params);

    // add strain energy
    virtual void AddStrainEnergy(
        double& psi,                      ///< strain energy function
        const LINALG::Matrix<3,1>& prinv, ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& modinv ///< modified invariants of right Cauchy-Green tensor
    );

    /// Add principal coefficients.
    void AddCoefficientsPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

    /// @name Access methods
    //@{
    double CLUM() const { return params_->clum_; }
    double CMED() const { return params_->cmed_; }
    double CABLUM() const { return params_->cablum_; }
    double NUE() const { return params_->nue_; }
    double BETA() const { return params_->beta_; }
    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isomodvisco, ///< global indicator for isotropic splitted viscos formulation
        bool& viscogenmax  ///< global indicator for viscous contribution according the SLS-Model
        )
    {
      isoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupAAAGasser* params_;

    /// for Gasser AAA-ILT material
    double normdist_;

    /// factor for material parameter
    double cele_;

  };

}  // namespace ELAST
}  // namespace MAT



#endif
