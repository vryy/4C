/*----------------------------------------------------------------------*/
/*!
\file visco_isoratedep.H
\brief
This file contains the routines required to calculate the isochoric contribution
of a viscos rate dependent material law, modified from Pioletti,1997

<pre>
Maintainer: Anna Birzle
            anna.birzle@tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef VISCO_ISORATEDEP_H_
#define VISCO_ISORATEDEP_H_

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a frequency independent viscoelastic material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 VISCO_IsoRateDep N 1
    class IsoRateDep : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoRateDep(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoRateDep() { ; }

      /// @name material parameters
      //@{

      /// material parameters
      double n_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetN(const double n) { n_ = n; }
      //@}

    };  // class IsoRateDep

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric general power material
  ///
  /// This is the isochoric part of a viscohyperelastic, isotropic
  /// material depending on the modified invariants of the rate of the right
  /// Cauchy-Green tensor.
  ///
  /// Strain energy function is given by
  ///\f[
  ///   \Psi = n \overline{J}_2 (\overline{I}_1 -3).
  ///\f]
  /// (n = \eta)
  ///
  /// \author birzle
  /// \date 09/13
  class IsoRateDep : public Summand
  {

  public:

    /// empty constructor
    IsoRateDep();

    /// constructor with given material parameters
    IsoRateDep(MAT::ELASTIC::PAR::IsoRateDep* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoratedep; }

    //@}

    /// Add modified coeffiencts.
    virtual void AddCoefficientsViscoModified(
        const LINALG::Matrix<3,1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
    	LINALG::Matrix<8,1>& modmy,		///<necassary coefficients for piola-kirchhoff-stress
    	LINALG::Matrix<33,1>& modxi,	///<necassary coefficients for viscosity tensor
        LINALG::Matrix<7,1>& modrateinv,
        Teuchos::ParameterList& params
    );

    /// @name Access methods
    //@{
    double N() const { return params_->n_; }

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
      isomod = true;
      isomodvisco = true;
      return;
    };

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::IsoRateDep* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ISORATEDEP_H_
