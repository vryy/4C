/*----------------------------------------------------------------------*/
/*!
\file elast_couplogneohooke.H
\brief
This file contains the routines required for Saint-Venant-Kirchhoff material

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_SVK_H_
#define ELAST_SVK_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    class CoupSVK : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
     CoupSVK(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CoupSVK() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      double mue_;
      /// Lame's constant
      double lambda_;

      //@}

      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

    };

  }  // namespace PAR

  class CoupSVK : public Summand
  {

  public:

    /// constructor with given material parameters
    CoupSVK(MAT::ELASTIC::PAR::CoupSVK* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_coupSVK; }

    //@}

    virtual void AddCoefficientsPrincipal(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    );

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isomodvisco, ///< global indicator for isotropic splitted viscos formulation
        bool& viscogenmax_  ///< global indicator for viscous contribution according the SLS-Model
      )
    {
      isoprinc = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::CoupSVK* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_SVK_H_
