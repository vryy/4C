/*----------------------------------------------------------------------*/
/*! \file
\brief the input line should read
  MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0

\level 1

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/
#ifndef ELAST_COUPANISOEXPOBASE_H_
#define ELAST_COUPANISOEXPOBASE_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/anisotropy_extension_base.H"
#include "../drt_mat/anisotropy_extension_provider.H"


namespace MAT
{
  namespace ELASTIC
  {
    // forward declarations
    class CoupAnisoExpoBase_TestSuite;

    /*!
     * \brief Pure abstract class to define the interface to the implementation of specific
     * materials
     *
     */
    class CoupAnisoExpoBaseInterface
    {
     public:
      ~CoupAnisoExpoBaseInterface() = default;

      /*!
       * \brief Returns the scalar product of the fibers
       *
       * \param gp Gauss point
       * \return double
       */
      virtual double GetScalarProduct(int gp) const = 0;

      /*!
       * \brief Returns the structural tensor that should be used
       *
       * \param gp Gauss point
       * \return const LINALG::Matrix<3, 3>&
       */
      virtual const LINALG::Matrix<3, 3>& GetStructuralTensor(int gp) const = 0;

      /*!
       * \brief Returns the structural tensor in stress like Voigt notation that should be used
       *
       * \param gp Gauss point
       * \return const LINALG::Matrix<6, 1>&
       */
      virtual const LINALG::Matrix<6, 1>& GetStructuralTensor_stress(int gp) const = 0;
    };

    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for coupled contribution of a anisotropic exponential fiber material
      ///
      /// <h3>Input line</h3>
      /// MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0
      /// ADAPT_ANGLE 0
      class CoupAnisoExpoBase
      {
        // This class is a friend of the corresponding unit test class
        friend class MAT::ELASTIC::CoupAnisoExpoBase_TestSuite;

       public:
        /// standard constructor
        explicit CoupAnisoExpoBase(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// @name material parameters
        //@{

        /// fiber params
        double k1_;
        double k2_;
        /// angle between circumferential and fiber direction (used for cir, axi, rad nomenclature)
        double gamma_;
        /// fiber params for the compressible case
        double k1comp_;
        double k2comp_;
        /// fiber initalization status
        int init_;

        //@}

        /// @name Write methods
        //@{
        void SetK1(const double k1) { k1_ = k1; }
        void SetK2(const double k2) { k2_ = k2; }
        void SetK1COMP(const double k1comp) { k1comp_ = k1comp; }
        void SetK2COMP(const double k2comp) { k2comp_ = k2comp; }
        //@}

       private:
        /// Constructor only used for unit testing
        CoupAnisoExpoBase();

      };  // class CoupAnisoExpoBase

    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /// Coupled anisotropic exponential fiber function, implemented for one possible fiber family
    /// [1]
    /// This is a hyperelastic, anisotropic material
    /// of the most simple kind.
    ///
    /// Strain energy function is given by
    ///\f[
    ///   \Psi = \frac {k_1}{2 k_2} \left(e^{k_2 (IV_{\boldsymbol C}-1)^2 }-1 \right).
    ///\f]
    ///
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] G.A. Holzapfel, T.C. Gasser, R.W. Ogden: A new constitutive framework for arterial
    /// wall mechanics and a comparative study of material models, J. of Elasticity 61 (2000) 1-48.
    /// </ul>
    ///
    /// \author st
    /// \date 03/12

    class CoupAnisoExpoBase : public Summand
    {
     public:
      /// constructor with given material parameters
      explicit CoupAnisoExpoBase(MAT::ELASTIC::PAR::CoupAnisoExpoBase* params);

      /*!
       * \brief Evaluate first derivative of the strain energy function with respect to the
       * anisotropic invariants.
       *
       * \param dPI_aniso (out) : First derivatives of the strain energy function with respect to
       * the anisotropic invariants
       *
       * \param C (in) : Cauchy Green deformation tensor
       *
       * \param gp (in) : Gauss point
       *
       * \param eleGID (in) : global element id
       */
      void EvaluateFirstDerivativesAniso(LINALG::Matrix<2, 1>& dPI_aniso,
          const LINALG::Matrix<3, 3>& C, int gp, int eleGID) override;

      /*!
       * \brief Evaluate second derivative of the strain energy function with respect to the
       * anisotropic invariants.
       *
       * \param ddPI_aniso (out) : Second derivatives of the strain energy function with respect to
       * the anisotropic invariants
       *
       * \param C (in) : Cauchy Green deformation tensor
       *
       * \param gp (in) : Gauss point
       *
       * \param eleGID (in) : global element id
       */
      void EvaluateSecondDerivativesAniso(LINALG::Matrix<3, 1>& ddPI_aniso,
          const LINALG::Matrix<3, 3>& C, int gp, int eleGID) override;

      /// retrieve coefficients of first, second and third derivative
      /// of summand with respect to anisotropic invariants
      ///
      /// The derivatives of the summand
      /// \f$\Psi(IV_{\boldsymbol{C},\boldsymbol{a}},V_{\boldsymbol{C},\boldsymbol{a}})\f$ in which
      /// the principal anisotropic invariants are the arguments are defined as following:
      ///
      /// First derivatives:
      ///
      ///\f[
      /// dPI_{0,aniso} = \frac{\partial \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      ///\f[
      /// dPI_{1,aniso} = \frac{\partial \Psi}{\partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      /// Second derivatives:
      ///\f[
      /// ddPII_{0,aniso} = \frac{\partial^2 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}^2} ;
      ///\f]
      ///\f[
      /// ddPII_{1,aniso} = \frac{\partial^2 \Psi}{\partial V_{\boldsymbol{C},\boldsymbol{a}}^2} ;
      ///\f]
      ///\f[
      /// ddPII_{2,aniso} = \frac{\partial^2 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
      /// \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      /// Third derivatives:
      ///\f[
      /// dddPIII_{0,aniso} = \frac{\partial^3 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
      /// \partial IV_{\boldsymbol{C},\boldsymbol{a}} \partial IV_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      /// dddPIII_{1,aniso} = \frac{\partial^3 \Psi}{\partial V_{\boldsymbol{C},\boldsymbol{a}}
      /// \partial V_{\boldsymbol{C},\boldsymbol{a}} \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      /// dddPIII_{2,aniso} = \frac{\partial^3 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
      /// \partial IV_{\boldsymbol{C},\boldsymbol{a}} \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      /// dddPIII_{3,aniso} = \frac{\partial^3 \Psi}{\partial IV_{\boldsymbol{C},\boldsymbol{a}}
      /// \partial V_{\boldsymbol{C},\boldsymbol{a}} \partial V_{\boldsymbol{C},\boldsymbol{a}}} ;
      ///\f]
      template <typename T>
      void GetDerivativesAniso(
          LINALG::Matrix<2, 1, T>& dPI_aniso,      ///< first derivative with respect to invariants
          LINALG::Matrix<3, 1, T>& ddPII_aniso,    ///< second derivative with respect to invariants
          LINALG::Matrix<4, 1, T>& dddPIII_aniso,  ///< third derivative with respect to invariants
          LINALG::Matrix<3, 3, T> const& rcg,      ///< right Cauchy-Green tensor
          int gp,                                  ///< Gauss point
          int eleGID) const;                       ///< element GID

      /// Add anisotropic principal stresses
      void AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,  ///< right Cauchy Green Tensor
          LINALG::Matrix<6, 6>& cmat,                                ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,                              ///< 2nd PK-stress
          Teuchos::ParameterList&
              params,  ///< additional parameters for computation of material properties
          int gp,      ///< Gauss point
          int eleGID   ///< element GID
          ) override;

      /// add strain energy
      void AddStrainEnergy(double& psi,       ///< strain energy functions
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>&
              glstrain,  ///< Green-Lagrange strain in strain like Voigt notation
          int gp,        //< Gauss point
          int eleGID     ///< element GID
          ) override;

      /// Evaluates strain energy for automatic differentiation with FAD
      template <typename T>
      void EvaluateFunc(T& psi,                ///< strain energy functions
          LINALG::Matrix<3, 3, T> const& rcg,  ///< Right Cauchy-Green tensor
          int gp,                              ///< Gauss point
          int eleGID) const;                   ///< element GID


      /// @name Access methods
      ///@{
      double K1() const { return params_->k1_; }
      double K2() const { return params_->k2_; }
      double K1COMP() const { return params_->k1comp_; }
      double K2COMP() const { return params_->k2comp_; }
      ///@}

      /// Set fiber directions
      void SetFiberVecs(double newgamma,       ///< new angle
          const LINALG::Matrix<3, 3>& locsys,  ///< local coordinate system
          const LINALG::Matrix<3, 3>& defgrd   ///< deformation gradient
          ) override;

      /// Set fiber directions
      void SetFiberVecs(const LINALG::Matrix<3, 1>& fibervec  ///< new fiber vector
          ) override;

      /// Get fiber directions
      void GetFiberVecs(
          std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
          ) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
          ) override
      {
        anisoprinc = true;
      };

     protected:
      /*!
       * \brief Get the structural tensor interface from the derived materials
       *
       * \return const CoupAnisoExpoBaseInterface& Interface that computes structural tensors and
       * scalar products
       */
      virtual const CoupAnisoExpoBaseInterface& GetCoupAnisoExpoBaseInterface() const = 0;

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::CoupAnisoExpoBase* params_;
    };

  }  // namespace ELASTIC
}  // namespace MAT

#endif  // ELAST_COUPANISOEXPOBASE_H_
