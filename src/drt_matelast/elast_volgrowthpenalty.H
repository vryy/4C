/*----------------------------------------------------------------------*/
/*!
\file elast_volgrowthpenalty.H
\brief

<pre>
Maintainer: Fabian Bräu
            braeu@lnm.mw.tum.de
            089/289 15236
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_VOLGROWTHPENALTY_H_
#define ELAST_VOLGROWTHPENALTY_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for volumetric contribution
    /// \f$\Psi=\epsilon \left( J - \hat{v} \right)^n\f$
    /// \hat{v}: desired volumetric change
    /// \epsilon: penalty parameter
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_VolGrowthPenalty EPSILON 1. EXPONENT n
    class VolGrowthPenalty : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      VolGrowthPenalty(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~VolGrowthPenalty() { ; }

      /// @name material parameters
      //@{

      /// Penalty parameter
      double eps_;
      /// Exponent
      double n_;

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

      /// @name Write methods
      //@{
      void SetEpsilon(const double eps) { eps_ = eps; }
      void SetExponent(const double n) { n_ = n; }
      //@}

    };  // class VolGrowthPenalty

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Volumetric Growth Penalty
  ///
  /// Strain energy function is given by
  ///
  ///
  ///\f[
  ///   \Psi=\epsilon \left( J - \hat{v} \right)
  ///\f]
  ///
  ///
  /// \author braeu
  /// \date 2015
  class VolGrowthPenalty : public Summand
  {

  public:

    /// constructor with given material parameters
    VolGrowthPenalty(MAT::ELASTIC::PAR::VolGrowthPenalty* params);

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_volgrowthpenalty; }

    //@}

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    // Add derivatives with respect to modified invariants.
    virtual void AddDerivativesModified(
        LINALG::Matrix<3,1>& dPmodI,       ///< first derivative with respect to modified invariants
        LINALG::Matrix<6,1>& ddmodPII,     ///< second derivative with respect to modified invariants
        const LINALG::Matrix<3,1>& modinv, ///< modified invariants of right Cauchy-Green tensor
        const int eleGID                   ///< element GID
        );

    /// @name Access methods
    //@{
    double Epsilon() const { return params_->eps_; }
    double Exponent() const { return params_->n_; }
    //@}

    /// Indicator for formulation
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& isovisco, ///< global indicator for isotropic viscous formulation
        bool& viscogenmax  ///< global indicator for viscous contribution according the SLS-Model
        )
    {
      isomod = true;
      return;
    };


  private:

    /// my material parameters
   MAT::ELASTIC::PAR::VolGrowthPenalty* params_;

   /// penalty parameter
   double eps_;

   /// exponent in strain energy function
   double n_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // VOLGROWTHPENALTY_H_
