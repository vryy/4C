/*----------------------------------------------------------------------*/
/*!
\file elast_isovolaaagasser.H
\brief

This file contains the routines required to calculate the isochoric contribution
of the aaagasser material and the corresponding volumetric contribution.

MAT 20 ELAST_IsoVolAAAGasser CLUM 2.62E3 CMED 1.98E3 CABLUM 1.73E3 NUE 0.49 BETA -2.0

<pre>
maintainer: Andreas Maier
            a.maier@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOVOLAAAGASSER_H_
#define ELAST_ISOVOLAAAGASSER_H_
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
    /// This material provides the isochoric and(!) the volumetric
    /// contribution for a AAA thromubs. The input line should read
    /// like this:
    ///     MAT 21 ELAST_IsoVolAAAGasser CLUM 2.62E3 CMED 1.98E3
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


    class IsoVolAAAGasser : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoVolAAAGasser(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoVolAAAGasser() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class IsoVolAAAGasser

  }  // namespace PAR


  class IsoVolAAAGasser : public Summand
  {

  public:

    /// empty constructor
    IsoVolAAAGasser();

    /// constructor with given material parameters
    IsoVolAAAGasser(MAT::ELASTIC::PAR::IsoVolAAAGasser* params);

    /// @name Access material constants
    //@{
    
    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isovolaaagasser; }

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
    
   
    /// Add principal coefficients for IsoVolAAAGasser material.
    void AddCoefficientsPrincGasser(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        LINALG::Matrix<3,1>& matparam,
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

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
    double CLUM() const { return params_->clum_; }
    double CMED() const { return params_->cmed_; }
    double CABLUM() const { return params_->cablum_; }
    double NUE() const { return params_->nue_; }
    double BETA() const { return params_->beta_; }
    //@}

    /// @name Write methods
    //@{
    //void SetCLUM(const double clum) { params_->clum_ = clum; }
    //void SetCMED(const double cmed) { params_->cmed_ = cmed; }
    //void SetCABLUM(const double cablum) { params_->cablum_ = cablum; }
    //void SetNUE(const double nue) { params_->nue_ = nue; }
    //void SetBETA(const double beta) { params_->beta_ = beta; }

    //@}
  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoVolAAAGasser* params_;


  };

}  // namespace ELAST
}  // namespace MAT



#endif  // CCADISCRET
#endif
