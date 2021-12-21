/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for an isotropic general power-type material in terms of the first
Cauchy-Green invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef ELAST_COUP1POW_H_
#define ELAST_COUP1POW_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for isochoric contribution of a general power material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_Coup1Pow C 1 D 1
      class Coup1Pow : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        Coup1Pow(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double c_;
        int d_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          dserror(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };

        /// @name Write methods
        //@{
        void SetC(const double c) { c_ = c; }
        void SetD(const int d) { d_ = d; }
        //@}

      };  // class Coup1Pow

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Isochoric general power material
    ///
    /// This is a summand of a variable order hyperelastic, isotropic
    /// material depending on the first invariant of the right Cauchy-Green tensor.
    ///
    /// Strain energy function is given by
    ///\f[
    ///   \Psi = C (I_{\boldsymbol{C}}-3)^D.
    ///\f]
    ///
    ///
    /// \author rausch
    /// \date 03/11
    class Coup1Pow : public Summand
    {
     public:
      /// constructor with given material parameters
      Coup1Pow(MAT::ELASTIC::PAR::Coup1Pow* params);

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::mes_coup1pow; }

      void CopyStatInvAnaMatParams(std::vector<Teuchos::RCP<Epetra_Vector>> input) override;

      /// Add parameters of elasthyper-summand for stat inverse analysis
      void SetStatInvAnaSummandMatParams() override;

      /// Add parameters of elasthyper-summand for stat inverse analysis
      void AddElastOptParams(std::map<std::string, int>* pnames) override;

      //@}

      // add strain energy
      void AddStrainEnergy(double& psi,       ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          int gp,                                ///< Gauss point
          int eleGID                             ///< element GID
          ) override;

      void AddDerivativesPrincipal(
          LINALG::Matrix<3, 1>& dPI,          ///< first derivative with respect to invariants
          LINALG::Matrix<6, 1>& ddPII,        ///< second derivative with respect to invariants
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          int gp,                             ///< Gauss point
          int eleGID                          ///< element GID
          ) override;

      /// add the derivatives of a coupled strain energy functions associated with a purely
      /// isochoric deformation
      void AddCoupDerivVol(
          const double j, double* dPj1, double* dPj2, double* dPj3, double* dPj4) override
      {
        dserror("not implemented");
      }

      /// @name Access methods
      //@{
      double C() const { return params_->c_; }
      int D() const { return params_->d_; }

      //@}

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        isoprinc = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::Coup1Pow* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // COUP1POW_H_
