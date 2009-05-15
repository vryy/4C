/*----------------------------------------------------------------------*/
/*!
\file elast_isomooneyrivlin.H
\brief

This file contains the routines required to calculate the isochoric contribution 
of a Mooney-Rivlin-type material

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOMOONEYRIVLIN_H_
#define ELAST_ISOMOONEYRIVLIN_H_
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
    /// material parameters for isochoric contribution of a Mooney-Rivlin type material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoMooneyRivlin C1 100 C2 50  
    class IsoMooneyRivlin : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoMooneyRivlin(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoMooneyRivlin() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      const double c1_;
      const double c2_;

      //@}

    };  // class IsoMooneyRivlin

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric Mooney-Rivline material according to [1]
  ///
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2 (\overline{II}_{\boldsymbol{C}}-3).
  ///\f]
  ///
  /// 
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch,tk
  /// \date 04/09
  class IsoMooneyRivlin : public Summand
  {

  public:

    /// empty constructor
    IsoMooneyRivlin();

    /// constructor with given material parameters
    IsoMooneyRivlin(MAT::ELAST::PAR::IsoMooneyRivlin* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isomooneyrivlin; }

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
    MAT::ELAST::PAR::IsoMooneyRivlin* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  
