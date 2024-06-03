/*----------------------------------------------------------------------*/
/*! \file
\brief growth law for the volumetric growth materials. This file contains routines needed for the
calculation of the volumetric growth parameter theta. These so called growthlaws are to be used for
the isovolumetric split growth materials MAT::Growth. Upon request a growthlaw delivers a growth
factor \f$\vartheta\f$ and its derivative wrt. \f$\frac{\partial \vartheta}{\partial C}\f$.

\level 2

*/

/*----------------------------------------------------------------------*/


#ifndef FOUR_C_MAT_GROWTH_LAW_HPP
#define FOUR_C_MAT_GROWTH_LAW_HPP

#include "4C_config.hpp"

#include "4C_inpar_material.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/* forward declarations */
namespace MAT
{
  class Growth;
  class So3Material;

  /*----------------------------------------------------------------------*/
  /*! \class GrowthLaw
      \brief GrowthLaw base class

      It provides the interfaces Evaluate(theta,dthetadC,...) to compute
      the growth factor \f$\vartheta\f$ and its derivative
      \f$\frac{\partial \vartheta}{\partial C}\f$ and
      EvaluateNonLinMass(...,linmass_disp,...) to account for the dynamic change
      of density.

      All growth laws can be modified by a function depending on scalar concentrations
      resulting from scalar transport equations.

     \author kehl
     \date 6/2015
   */
  class GrowthLaw
  {
   public:
    //! constructors
    GrowthLaw();

    //! explicit constructor
    explicit GrowthLaw(CORE::MAT::PAR::Parameter* params);

    //! destructor
    virtual ~GrowthLaw() = default;
    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     * \param In
     * gp - Gauss point
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void Evaluate(double* theta, const double& thetaold,
        CORE::LINALG::Matrix<6, 1>* dthetadC, MAT::Growth& matgrowth,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<6, 1>* glstrain,
        const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) = 0;

    /*! \brief derivative of theta wrt parameters
     *
     * \param Out
     * theta - theta diff. wrt. parameters given in params
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matelastic - an elastic material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 07/2015
     */
    virtual void EvaluatePDeriv(double* theta, const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        const int eleGID) = 0;

    /// calculate growth part of deformation gradient
    virtual void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) = 0;

    /// density scaling of specific growth law
    virtual double DensityScale(const double theta) = 0;

    /// density derivative scaling of specific growth law
    virtual double DensityDerivScale(const double theta) = 0;

    //@}

    //! @name Query methods
    //@{

    //! material type
    virtual CORE::Materials::MaterialType MaterialType() const = 0;

    //! return material parameters
    CORE::MAT::PAR::Parameter* Parameter() { return params_; }

    //! Return whether material has a varying material density
    virtual bool VaryingDensity() const = 0;

    //@}

   private:
    //! material parameters
    CORE::MAT::PAR::Parameter* params_;
  };

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawDyn
        \brief Common parameters for dynamic growth laws

        \author kehl
        \date 6/2015
   */
    class GrowthLawDyn : public CORE::MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      GrowthLawDyn(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// tolerance for local Newton iteration
      const double abstol_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

    };  // class Growth

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawDyn
   *  \brief GrowthLawDyn base class
   *
   *  Base functionality of a dynamic growth law, i.e. solution of a generic
   *  ode \f$\dot\vartheta = f_{rhs}(\vartheta)\f$, via implicit timestep discretization
   *  and application of some Newton-like method (see Diss. Susanna Tinkl).
   *
   *  The specific rhs-functionality \f$f_{rhs}\f$ must be provided by
   *  specializiations of this class. Therefore the pure virtual interface-functions
   *  EvaluateGrowthFuntion(...), evaluate_growth_function_deriv_theta(...) and
   *  evaluate_growth_function_deriv_c(...) are provided.
   *
   *  \author kehl
   *  \date 6/2015
   */
  class GrowthLawDyn : public GrowthLaw
  {
   public:
    //! construct empty material object
    GrowthLawDyn();
    explicit GrowthLawDyn(CORE::MAT::PAR::Parameter* params);


    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(double* theta, const double& thetaold, CORE::LINALG::Matrix<6, 1>* dthetadC,
        MAT::Growth& matgrowth, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    /*! \brief derivative of theta wrt parameters
     *
     * \param Out
     * theta - theta diff. wrt. parameters given in params
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matelastic - an elastic material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 07/2015
     */
    void EvaluatePDeriv(double* theta, const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        const int eleGID) override;

    //@}

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    virtual void evaluate_growth_function(
        double& growthfunc, const double growthtrig, const double theta) = 0;

    virtual void evaluate_growth_trigger(double& growthtrig,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<6, 1>& Cdachvec, const CORE::LINALG::Matrix<3, 1>& direction,
        const double& theta, const double& consttrig) = 0;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    virtual void evaluate_growth_function_deriv_theta(double& dgrowthfunctheta, double growthtrig,
        double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const CORE::LINALG::Matrix<3, 3>& F_g, const CORE::LINALG::Matrix<3, 1>& direction) = 0;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    virtual void evaluate_growth_function_deriv_c(
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const CORE::LINALG::Matrix<3, 1>& direction) = 0;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override = 0;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    //! return the material parameters
    MAT::PAR::GrowthLawDyn* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawDyn*>(MAT::GrowthLaw::Parameter());
    }

    //! add parameter to parameter list
    void add_params_to_parameter_list(Teuchos::ParameterList& params, const double& theta)
    {
      FOUR_C_THROW(
          "You should not get here, "
          "function is only intended for GrowthLawIso so far, feel free to implement it for others "
          "if you need to");
    }

    //! Return whether material has a varying material density
    bool VaryingDensity() const override { return true; }
  };

  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawStatic
     \brief GrowthLawDyn base class

     Interface class for a static growth law, i.e. solution of a generic
     expression \f$\vartheta = f\f$.

     The specific rhs-functionality \f$f\f$ must be provided by
     specializiations of this class.

     \author kehl
     \date 6/2015
   */
  class GrowthLawStatic : public GrowthLaw
  {
   public:
    //! construct empty material object
    GrowthLawStatic();
    explicit GrowthLawStatic(CORE::MAT::PAR::Parameter* params);

    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(double* theta, const double& thetaold, CORE::LINALG::Matrix<6, 1>* dthetadC,
        MAT::Growth& matgrowth, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) override = 0;

    /*! \brief derivative of theta wrt parameters
     *
     * \param Out
     * theta - theta diff. wrt. parameters given in params
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matelastic - an elastic material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 07/2015
     */
    void EvaluatePDeriv(double* theta, const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        const int eleGID) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override = 0;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    /// material type
    CORE::Materials::MaterialType MaterialType() const override = 0;

    //! Return whether material has a varying material density
    bool VaryingDensity() const override { return true; }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawAnisoStrain
         \brief Common parameters for growth law

         \author hirschvogel
         \date 4/2017
     */
    class GrowthLawAnisoStrain : public GrowthLawDyn
    {
     public:
      /// standard constructor
      GrowthLawAnisoStrain(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// growth time scale
      const double tau_;
      /// reverse growth time scale
      const double taurev_;
      /// minimal growth stretch
      const double thetamin_;
      /// maximal growth stretch
      const double thetamax_;
      /// growth law exponent
      const double gamma_;
      /// reverse growth law exponent
      const double gammarev_;
      /// critical myofiber stretch initiating longitudinal growth
      const double lambdacrit_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawAnisoStrain
   */
  class GrowthLawAnisoStrain : public GrowthLawDyn
  {
   public:
    /// construct empty material object
    GrowthLawAnisoStrain();
    explicit GrowthLawAnisoStrain(MAT::PAR::GrowthLawAnisoStrain* params);

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    void evaluate_growth_function(
        double& growthfunc, const double growthtrig, const double theta) override;

    void evaluate_growth_trigger(double& growthtrig,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<6, 1>& Cdachvec, const CORE::LINALG::Matrix<3, 1>& direction,
        const double& theta, const double& consttrig) override;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void evaluate_growth_function_deriv_theta(double& dgrowthfunctheta, double growthtrig,
        double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const CORE::LINALG::Matrix<3, 3>& F_g,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void evaluate_growth_function_deriv_c(CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double growthtrig, double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_aniso_strain;
    }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStrain* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawAnisoStrain*>(MAT::GrowthLawDyn::Parameter());
    }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawAnisoStress
         \brief Common parameters for growth law

         \author hirschvogel
         \date 4/2017
     */
    class GrowthLawAnisoStress : public GrowthLawDyn
    {
     public:
      /// standard constructor
      GrowthLawAnisoStress(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// growth time scale
      const double tau_;
      /// reverse growth time scale
      const double taurev_;
      /// minimal growth stretch
      const double thetamin_;
      /// maximal growth stretch
      const double thetamax_;
      /// growth law exponent
      const double gamma_;
      /// reverse growth law exponent
      const double gammarev_;
      /// critical cavity pressure initiating transverse growth
      const double pcrit_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawAnisoStress
   */
  class GrowthLawAnisoStress : public GrowthLawDyn
  {
   public:
    /// construct empty material object
    GrowthLawAnisoStress();
    explicit GrowthLawAnisoStress(MAT::PAR::GrowthLawAnisoStress* params);

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    void evaluate_growth_function(
        double& growthfunc, const double growthtrig, const double theta) override;

    void evaluate_growth_trigger(double& growthtrig,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<6, 1>& Cdachvec, const CORE::LINALG::Matrix<3, 1>& direction,
        const double& theta, const double& consttrig) override;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void evaluate_growth_function_deriv_theta(double& dgrowthfunctheta, double growthtrig,
        double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const CORE::LINALG::Matrix<3, 3>& F_g,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void evaluate_growth_function_deriv_c(CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double growthtrig, double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_aniso_stress;
    }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStress* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawAnisoStress*>(MAT::GrowthLawDyn::Parameter());
    }
  };


  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawAnisoStrainConstTrig
         \brief Common parameters for growth law

         \author hirschvogel
         \date 5/2017
     */
    class GrowthLawAnisoStrainConstTrig : public GrowthLawAnisoStrain
    {
     public:
      /// standard constructor
      GrowthLawAnisoStrainConstTrig(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw() override;

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawAnisoStrainConstTrig
   */
  class GrowthLawAnisoStrainConstTrig : public GrowthLawAnisoStrain
  {
   public:
    /// construct empty material object
    GrowthLawAnisoStrainConstTrig();
    explicit GrowthLawAnisoStrainConstTrig(MAT::PAR::GrowthLawAnisoStrainConstTrig* params);

    //! @name Evaluation methods for the ode's rhs
    //@{


    void evaluate_growth_trigger(double& growthtrig,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<6, 1>& Cdachvec, const CORE::LINALG::Matrix<3, 1>& direction,
        const double& theta, const double& consttrig) override;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void evaluate_growth_function_deriv_theta(double& dgrowthfunctheta, double growthtrig,
        double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const CORE::LINALG::Matrix<3, 3>& F_g,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void evaluate_growth_function_deriv_c(CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double growthtrig, double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    //@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_aniso_strain_const_trig;
    }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStrainConstTrig* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawAnisoStrainConstTrig*>(MAT::GrowthLawDyn::Parameter());
    }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawAnisoStressConstTrig
         \brief Common parameters for growth law

         \author hirschvogel
         \date 5/2017
     */
    class GrowthLawAnisoStressConstTrig : public GrowthLawAnisoStress
    {
     public:
      /// standard constructor
      GrowthLawAnisoStressConstTrig(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw() override;

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawAnisoStressConstTrig
   */
  class GrowthLawAnisoStressConstTrig : public GrowthLawAnisoStress
  {
   public:
    /// construct empty material object
    GrowthLawAnisoStressConstTrig();
    explicit GrowthLawAnisoStressConstTrig(MAT::PAR::GrowthLawAnisoStressConstTrig* params);

    //! @name Evaluation methods for the ode's rhs
    //@{


    void evaluate_growth_trigger(double& growthtrig,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<6, 1>& Cdachvec, const CORE::LINALG::Matrix<3, 1>& direction,
        const double& theta, const double& consttrig) override;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void evaluate_growth_function_deriv_theta(double& dgrowthfunctheta, double growthtrig,
        double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const CORE::LINALG::Matrix<3, 3>& F_g,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void evaluate_growth_function_deriv_c(CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double growthtrig, double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    //@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_aniso_stress_const_trig;
    }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStressConstTrig* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawAnisoStressConstTrig*>(MAT::GrowthLawDyn::Parameter());
    }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawIsoStress
         \brief Common parameters for linear growth law

         \author kehl
         \date 6/2015
     */
    class GrowthLawIsoStress : public GrowthLawDyn
    {
     public:
      /// standard constructor
      GrowthLawIsoStress(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// maximal growth stretch
      const double thetaplus_;
      /// growth law parameter kthetaplus
      const double kthetaplus_;
      /// growth law parameter mthetaplus
      const double mthetaplus_;
      /// minimal growth stretch
      const double thetaminus_;
      /// growth law parameter kthetaminus
      const double kthetaminus_;
      /// growth law parameter mthetaminus
      const double mthetaminus_;
      /// homeostatic value for mandelstress
      const double hommandel_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawIsoStress
     \brief GrowthLaw based on a linear growth function

     The rhs-function \f$ f_{rhs}\f$ (see \sa GrowthLawDyn) is implemented from
     'Himpel,G. Kuhl, E. Menzel, A. & Steinmann, P. Computational modelling
     of isotropic multiplicative growth, Computer Modeling in Engineering
     and Sciences, 2005, 8, 119-134'. See also Diss. Susanna Tinkl.

     \author kehl
     \date 6/2015
   */
  class GrowthLawIsoStress : public GrowthLawDyn
  {
   public:
    /// construct empty material object
    GrowthLawIsoStress();
    explicit GrowthLawIsoStress(MAT::PAR::GrowthLawIsoStress* params);

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    void evaluate_growth_function(
        double& growthfunc, const double growthtrig, const double theta) override;

    void evaluate_growth_trigger(double& growthtrig,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<6, 1>& Cdachvec, const CORE::LINALG::Matrix<3, 1>& direction,
        const double& theta, const double& consttrig) override;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void evaluate_growth_function_deriv_theta(double& dgrowthfunctheta, double growthtrig,
        double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const CORE::LINALG::Matrix<3, 3>& F_g,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * growthtrig - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void evaluate_growth_function_deriv_c(CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double growthtrig, double theta, const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const CORE::LINALG::Matrix<3, 1>& direction) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_iso_stress;
    }

    //! material parameters
    MAT::PAR::GrowthLawIsoStress* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawIsoStress*>(MAT::GrowthLawDyn::Parameter());
    }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawAC : public CORE::MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      GrowthLawAC(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// maximal growth stretch
      const int Sc1_;
      /// growth law parameter kthetaplus
      const double alpha_;
      /// growth law parameter mthetaplus
      const int Sc2_;
      /// minimal growth stretch
      const double beta_;
      //@}


      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();


    };  // class Growth

  }  // namespace PAR

  class GrowthLawAC : public GrowthLawStatic
  {
   public:
    /// construct empty material object
    GrowthLawAC();
    explicit GrowthLawAC(MAT::PAR::GrowthLawAC* params);



    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(double* theta, const double& thetaold, CORE::LINALG::Matrix<6, 1>* dthetadC,
        MAT::Growth& matgrowth, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_ac;
    }

    MAT::PAR::GrowthLawAC* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawAC*>(MAT::GrowthLawStatic::Parameter());
    }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawACRadial : public GrowthLawAC
    {
     public:
      /// standard constructor
      GrowthLawACRadial(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawACRadial : public GrowthLawStatic
  {
   public:
    /// construct empty material object
    GrowthLawACRadial();
    explicit GrowthLawACRadial(MAT::PAR::GrowthLawAC* params);

    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(double* theta, const double& thetaold, CORE::LINALG::Matrix<6, 1>* dthetadC,
        MAT::Growth& matgrowth, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_ac_radial;
    }

    MAT::PAR::GrowthLawACRadial* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawACRadial*>(MAT::GrowthLawStatic::Parameter());
    }
  };


  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawACRadialRefConc : public GrowthLawAC
    {
     public:
      /// standard constructor
      GrowthLawACRadialRefConc(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawACRadialRefConc : public GrowthLawStatic
  {
   public:
    /// construct empty material object
    GrowthLawACRadialRefConc();
    explicit GrowthLawACRadialRefConc(MAT::PAR::GrowthLawAC* params);

    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(double* theta, const double& thetaold, CORE::LINALG::Matrix<6, 1>* dthetadC,
        MAT::Growth& matgrowth, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_ac_radial_refconc;
    }

    MAT::PAR::GrowthLawACRadialRefConc* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawACRadialRefConc*>(MAT::GrowthLawStatic::Parameter());
    }
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawConst : public CORE::MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      GrowthLawConst(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      //! enum for mapping between material parameter and entry in Parameter::matparams_
      enum Matparamnames
      {
        thetarate,
        first = thetarate,
        last = thetarate
      };
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  class GrowthLawConst : public GrowthLawStatic
  {
   public:
    /// construct empty material object
    GrowthLawConst();
    explicit GrowthLawConst(MAT::PAR::GrowthLawConst* params);

    //! @name Evaluation methods
    //@{

    /*! \brief evaluate the volumetric growth factor
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matgrowth - a reference to a growth material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(double* theta, const double& thetaold, CORE::LINALG::Matrix<6, 1>* dthetadC,
        MAT::Growth& matgrowth, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    /*! \brief derivative of theta wrt parameters
     *
     * \param Out
     * theta - theta diff. wrt. parameters given in params
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matelastic - an elastic material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 07/2015
     */
    void EvaluatePDeriv(double* theta, const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        const int eleGID) override;

    /// calculate growth part of deformation gradient
    void CalcFg(const double& theta, const double& thetaold, const int& gp,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 1>& refdir,
        const std::vector<CORE::LINALG::Matrix<3, 1>>& curdir,
        const std::vector<CORE::LINALG::Matrix<3, 3>>& histdefgrd,
        CORE::LINALG::Matrix<3, 3>& F_g) override;

    /// density scaling of specific growth law
    double DensityScale(const double theta) override;

    /// density derivative scaling of specific growth law
    double DensityDerivScale(const double theta) override;

    //@}

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_growth_const;
    }

    MAT::PAR::GrowthLawConst* Parameter()
    {
      return dynamic_cast<MAT::PAR::GrowthLawConst*>(MAT::GrowthLawStatic::Parameter());
    }
  };
}  // namespace MAT



FOUR_C_NAMESPACE_CLOSE

#endif
