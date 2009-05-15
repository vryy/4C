/*----------------------------------------------------------------------*/
/*!
\file elast_volsussmanbathe.H
\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite 
element analysis", Cambridge, 1997.

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            {rausch,kloeppel}@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_VOLSUSSMANBATHE_H_
#define ELAST_VOLSUSSMANBATHE_H_
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
namespace ELAST
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for volumetric contribution \f$\Psi=\kappa(J-1)^2\f$
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_VolSussmanBathe KAPPA 100
    class VolSussmanBathe : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      VolSussmanBathe(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~VolSussmanBathe() { ; }

      /// @name material parameters
      //@{

      /// Dilatation modulus
      const double kappa_;

      //@}

    };  // class VolSussmanBathe

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Volumetric SussmanBathe material according to [1]
  /// Strain energy function is given by
  ///
  ///
  /// \f[
  ///   \Psi = \kappa (J-1)^2
  /// \f]
  ///
  ///
  /// [1] Doll, S. and Schweizerhof, K. On the Development of Volumetric Strain Energy Functions
  ///     Journal of Applied Mechanics, 2000
  /// 
  class VolSussmanBathe : public Summand
  {

  public:

    /// empty constructor
    VolSussmanBathe();

    /// constructor with given material parameters
    VolSussmanBathe(MAT::ELAST::PAR::VolSussmanBathe* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_volsussmanbathe; }

    //@}

    /// Add principal coeffiencts.
    /// Nothing to do for this (decoupled) material!
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;

      return;
    };
    
    /// Add modified coeffiencts.
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    );

  private:

    /// my material parameters
    MAT::ELAST::PAR::VolSussmanBathe* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // VOLSUSSMANBATHE_H_
