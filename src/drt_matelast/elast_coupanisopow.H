/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the routines required for the anisotropic material using a pow-like
material

\level 3

\maintainer Amadeus Gebauer
*/
/*---------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOPOW_H_
#define ELAST_COUPANISOPOW_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for anisochoric contribution of a
      /// pow like material with one fiber direction
      ///
      ///
      /// <h3>Input line</h3>
      /// MAT 1 CoupAnisoPow C 1.0 D 2.0 [ GAMMA 35.0 INIT 0 ADAPT_ANGLE 0]
      class CoupAnisoPow : public MAT::PAR::ParameterAniso
      {
       public:
        /// standard constructor
        CoupAnisoPow(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~CoupAnisoPow() { ; }

        /// @name material parameters
        //@{

        /// stiffness factor
        double k_;
        /// exponential factor for I4
        double d1_;
        /// exponential factor for ((I4)^{d1}-1)
        double d2_;
        /// fiber number (1,2,3,...) used later as FIBER1,FIBER2,FIBER3,...
        int fibernumber_;
        /// Deformation threshold where fibers are active.
        double activethres_;
        /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
        double gamma_;
        /// fiber initalization status
        int init_;
        /// adapt angle during remodeling
        int adapt_angle_;

        //@}

        /// create material instance of matching type with my parameters
        ///
        /// Remark:                                                  AN 03/2014
        /// Function has to be implemented, since it is a pure virtual
        /// function in the MAT::PAR::Material class.
        /// Function stays empty, since the call within the material factory is done only
        /// respective the elasthyper class
        ///
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{
        //      void SetK(const double k) { k_ = k; };
        //      void SetD1(const double d) { d_ = d; };
        //@}

      };  // class CoupAnisoPow

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Coupled anisotropic pow-like fiber function, implemented for one possible fiber family
    /// [1]
    ///
    /// Strain energy function is given by
    ///\f[
    ///   \Psi = K \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2}.
    ///\f]
    ///
    /// The corresponding derivatives are
    ///
    ///\f[
    ///   \frac{d\Psi}{d IV_{\boldsymbol C}} = K\ D2\ D1\ (IV_{\boldsymbol C})^{D1-1}
    ///   \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2-1}.
    ///\f]
    ///
    ///\f[
    ///   \frac{d^2\Psi}{d^2 IV_{\boldsymbol C}} = K\ D2\ (D2-1)\ \left( D1\ (IV_{\boldsymbol
    ///   C})^{D1-1} \right)^2 \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2-2}\ +\ K\ D2\ D1
    ///   (D1-1)\ (IV_{\boldsymbol C})^{D1-2}\ \left((IV_{\boldsymbol C})^{D1}-1\right)^{D2-1}
    ///\f]
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] GA Holzapfel, Nonlinear solid mechanics 2004
    /// </ul>
    ///
    /// \author nagler
    /// \date 07/13
    class CoupAnisoPow : public Summand
    {
     public:
      /// constructor with given material parameters
      CoupAnisoPow(MAT::ELASTIC::PAR::CoupAnisoPow* params);

      ///@name Packing and Unpacking
      //@{

      virtual void PackSummand(DRT::PackBuffer& data) const;

      virtual void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position);

      //@}


      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_coupanisoneohooke;
      }

      //@}

      /// Setup of summand
      virtual void Setup(DRT::INPUT::LineDefinition* linedef);

      /// Add anisotropic principal stresses
      virtual void AddStressAnisoPrincipal(
          const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          LINALG::Matrix<6, 6>& cmat,       ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,     ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,       ///< additional parameters for computation of material properties
          const int eleGID  ///< element GID
      );

      /// @name Access material constants
      //@{
      /// Get stiffness
      //    double C() const { return params_->c_; };
      /// Get exponent factor
      //    double D() const { return params_->d_; };

      //@}

      /// Set fiber directions
      virtual void SetFiberVecs(const double newgamma,  ///< new angle
          const LINALG::Matrix<3, 3>& locsys,           ///< local coordinate system
          const LINALG::Matrix<3, 3>& defgrd            ///< deformation gradient
      );

      /// Get fiber directions
      virtual void GetFiberVecs(
          std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
      );

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
      )
      {
        anisoprinc = true;
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoPow* params_;

      /// fiber direction
      LINALG::Matrix<3, 1> a_;

      /// structural tensors in Voigt notation for anisotropy
      LINALG::Matrix<6, 1> A_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // COUPANISOPOW_H_
