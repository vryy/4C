/*----------------------------------------------------------------------*/
/*!
\file elast_couplogneohooke.H
\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite 
element analysis", Cambridge, 1997.

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_LOGNEOHOOKE_H_
#define ELAST_LOGNEOHOOKE_H_
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
    /// material parameters for logarithmic neo-Hooke material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_CoupLogNeoHooke YOUNG 1.044E7 NUE 0.3
    class CoupLogNeoHooke : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
     CoupLogNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupLogNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double nue_;

      //@}

    };  // class CoupLogNeoHooke

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Logarithmic neo-Hooke material according to [1]
  ///
  /// This is a compressible, hyperelastic, isotropic material 
  /// of the most simple kind.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] J Bonet and RD Wood, "Nonlinear continuum mechanics for finite 
  ///          element analysis", Cambridge, 1997.
  /// <li> [2] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author bborn
  /// \date 04/09
  class CoupLogNeoHooke : public Summand
  {

  public:

    /// empty constructor
    CoupLogNeoHooke();

    /// constructor with given material parameters
    CoupLogNeoHooke(MAT::ELASTIC::PAR::CoupLogNeoHooke* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_couplogneohooke; }

    //@}

    /// Main material call to determine  2nd PK stress and material constitutive tensor
    ///
    /// The material strain energy density function is
    ///\f[
    /// \Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3)
    ///       - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
    ///      + \frac{\lambda}{2} \big( \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \big)^2
    ///\f]
    /// which was taken from [1].
    ///
    /// Here is \f$I_{\boldsymbol{C}}\f$ the first principal invariant
    /// of the right Cauchy--Green strain tensor \f$\boldsymbol{C}\f$ 
    /// and \f$I\!I\!I_{\boldsymbol{C}}\f$ its third.
    /// The isochoric part is proportional to \f$\mu\f$,
    /// whereas the \f$\lambda\f$-proportional part constitutes the volumetric contribution.
    ///
    /// The parameters are the shear modulus
    ///\f[
    /// \mu = \frac{E}{2(1+\nu)}
    ///\f]
    /// and Lame's coefficient
    ///\f[
    /// \lambda = \left\{\begin{array}{ll}
    ///      \frac{E \nu}{(1+\nu) (1-2\nu)} & \text{if $\nu \neq 1/2$}
    ///   \\ 0                              & \text{else}
    /// \end{array}\right.
    ///\f]
    ///
    /// The 2nd Piola--Kirchhoff stress is
    ///\f[
    /// \boldsymbol{S} = \mu \big( \boldsymbol{1} - \boldsymbol{C}^{-1} \big)
    ///                + \lambda \, \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \,  \boldsymbol{C}^{-1}
    ///\f]
    ///
    /// The material constitutive 4-tensor
    /// \f$\boldsymbol{C}_\text{m}=C_{IJKL}\boldsymbol{E}^I\otimes\boldsymbol{E}^J\otimes\boldsymbol{E}^K\otimes\boldsymbol{E}^L\f$
    /// is determined by
    ///\f[
    /// C_{IJKL} 
    /// = \lambda (\boldsymbol{C}^{-1})_{IJ} \, (\boldsymbol{C}^{-1})_{KL}
    /// + 2\big(\mu-\lambda \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})\big)
    ///   \, \frac{1}{2} \big( \boldsymbol{C}^{-1})_{IK} \, (\boldsymbol{C}^{-1})_{JL}
    ///                        + (\boldsymbol{C}^{-1})_{IL} \, (\boldsymbol{C}^{-1})_{JK} \big)
    ///\f]
    /// or
    ///\f[
    /// \boldsymbol{C}_\text{m}
    /// = \lambda \boldsymbol{C}^{-1} \otimes  \boldsymbol{C}^{-1}
    /// + 2\big(\mu-\lambda \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})\big)
    ///   \, \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
    ///\f]
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );
    
    /// Add modified coeffiencts.
    /// Nothing to do for this (coupled) material!
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;
    };
        
    /// Add anisotropic principal coeffiencts.
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
        
    /// Add anisotropic modified coeffiencts.
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
    

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupLogNeoHooke* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // LOGNEOHOOKE_H_
