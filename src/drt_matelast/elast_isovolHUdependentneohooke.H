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
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


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


    /// @name Access material constants
    //@{
    
    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isovolHUdependentneohooke; }

    //@}

    /// Add principal coefficients.
    /// Nothing to do for this (decoupled) material!
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };
    
    /// Add principal calcified coefficients for IsoVolHUDependentNeoHooke material.
    void AddCoefficientsPrincCalcified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        LINALG::Matrix<3,1>& matparam,
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

    /// Add principal coefficients for IsoVolAAAGasser material.
    void AddCoefficientsPrincGasser(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        LINALG::Matrix<3,1>& matparam,
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// Add modified coefficients.
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };


    /// Add anisotropic principal coefficients.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsPrincipalAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };

    /// Add anisotropic modified coefficients.
    /// Nothing to do for this (isotropic) material!
    void AddCoefficientsModifiedAniso(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<15,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<6,1>& inv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };


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
  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT


#endif  // CCADISCRET
#endif
