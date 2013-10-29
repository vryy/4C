/*----------------------------------------------------------------------*/
/*!
\file elast_isovolHUdependentneohooke.H
\brief

This file contains the routines required to calculate the isochoric and volumetric contribution
of the HU dependent neo hooke material. E.g used for inclusion of
calcification in AAA simulations

<pre>
maintainer: Andreas Maier
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOVOLHUDEPENDENTNEOHOOKE_H_
#define ELAST_ISOVOLHUDEPENDENTNEOHOOKE_H_

// Linalg matrices used instead; included in elast_summand
//#include <Epetra_SerialDenseVector.h>
//#include <Epetra_SerialDenseMatrix.h>
// included in elast_summand
//#include "../drt_inpar/inpar_material.H"
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
    /// MAT 20 ELAST_IsoVolHUDependentNeoHooke ALPHA_MAX 8.929E6 CT_MIN 30.0 CT_MAX 600.0 NUE 0.49 BETA -2.0
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

    class IsoVolHUDependentNeoHooke : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoVolHUDependentNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class IsoVolHUDependentNeoHooke

  }  // namespace PAR


  class IsoVolHUDependentNeoHooke : public Summand
  {

  public:

    /// empty constructor
    IsoVolHUDependentNeoHooke();

    /// constructor with given material parameters
    IsoVolHUDependentNeoHooke(MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isovolHUdependentneohooke; }

    //@}

    /// Setup of summand
    virtual void Setup(DRT::INPUT::LineDefinition* linedef);

    /// Add principal calcified coefficients for IsoVolHUDependentNeoHooke material.
    void AddCoefficientsPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
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
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isomodvisco ///< global indicator for isotropic splitted viscos formulation
        )
    {
      isoprinc = true;
      return;
    };


  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params_;

    /// factor for material parameter
    double alpha_;
  };

}  // namespace ELAST
}  // namespace MAT


#endif
