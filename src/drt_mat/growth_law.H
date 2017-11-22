/*----------------------------------------------------------------------*/
/*!
\file growth_law.H

\brief growth law for the volumetric growth materials. This file contains routines needed for the calculation of the volumetric growth
parameter theta. These so called growthlaws are to be used for the isovolumetric
split growth materials MAT::Growth. Upon request a growthlaw delivers a growth factor \f$\vartheta\f$ and its derivative
wrt. \f$\frac{\partial \vartheta}{\partial C}\f$.

\level 2

<pre>
   \maintainer Marc Hirschvogel, Moritz Thon
               hirschvogel@mhpc.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10363
</pre>
*/

/*----------------------------------------------------------------------*/


#ifndef GROWTH_LAW_H
#define GROWTH_LAW_H

/* forward declarations */
namespace MAT
{

class Growth;

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
    explicit GrowthLaw(MAT::PAR::Parameter* params);

    //! destructor
    virtual ~GrowthLaw() {;}

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
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void Evaluate(double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        const double& consttrig,
        Teuchos::ParameterList& params,
        const int eleGID ) = 0;

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
    virtual void EvaluatePDeriv(
      double* theta,
      const double& thetaold,
      Teuchos::RCP<MAT::So3Material> matelastic,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      Teuchos::ParameterList& params,
      const int eleGID
    ) = 0;

    ///calculate growth part of deformation gradient
    virtual void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    ) = 0;

    ///density scaling of specific growth law
    virtual double DensityScale(const double theta) = 0;

    ///density derivative scaling of specific growth law
    virtual double DensityDerivScale(const double theta) = 0;

    //@}

    //! @name Query methods
    //@{

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    //! return material parameters
    MAT::PAR::Parameter* Parameter() {return params_;}

    // add parameter to parameter list
    virtual void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) = 0;

    //! Return whether material has a varying material density
    virtual bool VaryingDensity() const = 0;

    //@}

  private:

    //! material parameters
    MAT::PAR::Parameter* params_;

  };

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class GrowthLawDyn
        \brief Common parameters for dynamic growth laws

        \author kehl
        \date 6/2015
   */
    class GrowthLawDyn : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawDyn(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawDyn() { ; }

      /// @name material parameters
      //@{
      /// tolerance for local Newton iteration
      const double abstol_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

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
   *  EvaluateGrowthFuntion(...), EvaluateGrowthFunctionDerivTheta(...) and
   *  EvaluateGrowthFunctionDerivC(...) are provided.
   *
   *  \author kehl
   *  \date 6/2015
   */
  class GrowthLawDyn : public GrowthLaw
  {
  public:

    //! construct empty material object
    GrowthLawDyn();
    explicit GrowthLawDyn(MAT::PAR::Parameter* params);

    //! destructor
    virtual ~GrowthLawDyn(){;}


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
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void Evaluate
    (
        double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        const double& consttrig,
        Teuchos::ParameterList& params,
        const int eleGID
    );

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
    virtual void EvaluatePDeriv
    (
      double* theta,
      const double& thetaold,
      Teuchos::RCP<MAT::So3Material> matelastic,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    //@}

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        const double growthtrig,
        const double theta
    ) = 0;

    virtual void EvaluateGrowthTrigger
    (
        double& growthtrig,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<6,1>& Cdachvec,
        const LINALG::Matrix<3,1>& direction,
        const double& theta,
        const double& consttrig
    ) = 0;

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const LINALG::Matrix<3,3>& F_g,
        const LINALG::Matrix<3,1>& direction
    ) = 0;

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const LINALG::Matrix<3,1>& direction
    ) = 0;

    ///calculate growth part of deformation gradient
    virtual void CalcFg
    (
        const double& theta,
        const double& thetaold,
        const int& gp,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        LINALG::Matrix<3,3>& F_g
    ) = 0;

    ///density scaling of specific growth law
    virtual double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    virtual double DensityDerivScale(const double theta);

    //@}

    //! return the material parameters
    MAT::PAR::GrowthLawDyn* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawDyn*>(MAT::GrowthLaw::Parameter()); }

    //! add parameter to parameter list
    void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) {dserror("You should not get here, "
        "function is only intended for GrowthLawIso so far, feel free to implement it for others if you need to");}

    //! Return whether material has a varying material density
    bool VaryingDensity() const {return true;}
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
    explicit GrowthLawStatic(MAT::PAR::Parameter* params);

    virtual ~GrowthLawStatic(){;}

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
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void Evaluate
    (
        double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        const double& consttrig,
        Teuchos::ParameterList& params,
        const int eleGID
    ) = 0;

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
    virtual void EvaluatePDeriv
    (
      double* theta,
      const double& thetaold,
      Teuchos::RCP<MAT::So3Material> matelastic,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    ///calculate growth part of deformation gradient
    virtual void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    ) = 0;

    ///density scaling of specific growth law
    virtual double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    virtual double DensityDerivScale(const double theta);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    //! add parameter to parameter list
    virtual void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) = 0;

    //! Return whether material has a varying material density
    virtual bool VaryingDensity() const = 0;
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
      GrowthLawAnisoStrain(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawAnisoStrain() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawAnisoStrain
   */
  class GrowthLawAnisoStrain: public GrowthLawDyn
  {

  public:

    /// construct empty material object
    GrowthLawAnisoStrain();
    explicit GrowthLawAnisoStrain(MAT::PAR::GrowthLawAnisoStrain* params);

    /// destructor
    virtual ~GrowthLawAnisoStrain() { ; }

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        const double growthtrig,
        const double theta
    );

    virtual void EvaluateGrowthTrigger
    (
        double& growthtrig,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<6,1>& Cdachvec,
        const LINALG::Matrix<3,1>& direction,
        const double& theta,
        const double& consttrig
    );

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const LINALG::Matrix<3,3>& F_g,
        const LINALG::Matrix<3,1>& direction
    );

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const LINALG::Matrix<3,1>& direction
    );

    ///calculate growth part of deformation gradient
    virtual void CalcFg
    (
        const double& theta,
        const double& thetaold,
        const int& gp,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    virtual double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    virtual double DensityDerivScale(const double theta);

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_aniso_strain; }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStrain* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAnisoStrain*>(MAT::GrowthLawDyn::Parameter()); }

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
      GrowthLawAnisoStress(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawAnisoStress() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawAnisoStress
   */
  class GrowthLawAnisoStress: public GrowthLawDyn
  {

  public:

    /// construct empty material object
    GrowthLawAnisoStress();
    explicit GrowthLawAnisoStress(MAT::PAR::GrowthLawAnisoStress* params);

    /// destructor
    virtual ~GrowthLawAnisoStress() { ; }

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        const double growthtrig,
        const double theta
    );

    virtual void EvaluateGrowthTrigger
    (
        double& growthtrig,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<6,1>& Cdachvec,
        const LINALG::Matrix<3,1>& direction,
        const double& theta,
        const double& consttrig
    );

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const LINALG::Matrix<3,3>& F_g,
        const LINALG::Matrix<3,1>& direction
    );

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const LINALG::Matrix<3,1>& direction
    );

    ///calculate growth part of deformation gradient
    virtual void CalcFg
    (
        const double& theta,
        const double& thetaold,
        const int& gp,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    virtual double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    virtual double DensityDerivScale(const double theta);

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_aniso_stress; }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStress* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAnisoStress*>(MAT::GrowthLawDyn::Parameter()); }

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
      GrowthLawAnisoStrainConstTrig(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawAnisoStrainConstTrig() { ; }

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

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

    /// destructor
    virtual ~GrowthLawAnisoStrainConstTrig() { ; }

    //! @name Evaluation methods for the ode's rhs
    //@{


    void EvaluateGrowthTrigger
    (
        double& growthtrig,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<6,1>& Cdachvec,
        const LINALG::Matrix<3,1>& direction,
        const double& theta,
        const double& consttrig
    );

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const LINALG::Matrix<3,3>& F_g,
        const LINALG::Matrix<3,1>& direction
    );

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const LINALG::Matrix<3,1>& direction
    );

    //@}

    //! material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_aniso_strain_const_trig; }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStrainConstTrig* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAnisoStrainConstTrig*>(MAT::GrowthLawDyn::Parameter()); }

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
      GrowthLawAnisoStressConstTrig(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawAnisoStressConstTrig() { ; }

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

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

    /// destructor
    virtual ~GrowthLawAnisoStressConstTrig() { ; }

    //! @name Evaluation methods for the ode's rhs
    //@{


    void EvaluateGrowthTrigger
    (
        double& growthtrig,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<6,1>& Cdachvec,
        const LINALG::Matrix<3,1>& direction,
        const double& theta,
        const double& consttrig
    );

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const LINALG::Matrix<3,3>& F_g,
        const LINALG::Matrix<3,1>& direction
    );

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const LINALG::Matrix<3,1>& direction
    );

    //@}

    //! material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_aniso_stress_const_trig; }

    //! material parameters
    MAT::PAR::GrowthLawAnisoStressConstTrig* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAnisoStressConstTrig*>(MAT::GrowthLawDyn::Parameter()); }

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
      GrowthLawIsoStress(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawIsoStress() { ; }

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
      Teuchos::RCP<MAT::Material> CreateMaterial();

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
  class GrowthLawIsoStress: public GrowthLawDyn
  {

  public:

    /// construct empty material object
    GrowthLawIsoStress();
    explicit GrowthLawIsoStress(MAT::PAR::GrowthLawIsoStress* params);

    /// destructor
    virtual ~GrowthLawIsoStress() { ; }

    //! @name Evaluation methods for the ode's rhs
    //@{

    /*! Evaluate growth function
     * \param Out
     * growthfunc - the value of the growthfunction
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     *
     */
    void EvaluateGrowthFunction
    (
        double & growthfunc,
        const double growthtrig,
        const double theta
    );

    void EvaluateGrowthTrigger
    (
        double& growthtrig,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<6,1>& Cdachvec,
        const LINALG::Matrix<3,1>& direction,
        const double& theta,
        const double& consttrig
    );

    /*! Evaluate growth function wrt growth factor theta
     * \param Out
     * dgrowthfunctheta - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * Cdach - cauchy green tensor of the elastic part of the deformation
     * \param In
     * cmatelastic - elastic part of the cmat
     *
     */
    void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
        const LINALG::Matrix<3,3>& F_g,
        const LINALG::Matrix<3,1>& direction
    );

    /*! Evaluate growth function wrt right cauchy green tensor
     * \param Out
     * dgrowthfuncdC - derivation of the growthfunction wrt theta
     * \param In
     * traceM - the trace of the mandel stress
     * \param In
     * theta - the current growth factor theta
     * \param In
     * C - cauchy green tensor of the whole deformation
     * \param In
     * S - the stresses
     * \param In
     * cmatelastic - elastic part of the cmat
     */
    void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        const LINALG::Matrix<3,1>& direction
    );

    ///calculate growth part of deformation gradient
    void CalcFg
    (
        const double& theta,
        const double& thetaold,
        const int& gp,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<3,1>& refdir,
        const std::vector<LINALG::Matrix<3,1> >& curdir,
        const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
        LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    double DensityDerivScale(const double theta);

    //@}

    //! material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_iso_stress; }

    //! material parameters
    MAT::PAR::GrowthLawIsoStress* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawIsoStress*>(MAT::GrowthLawDyn::Parameter()); }

  };


















  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawAC : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawAC(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawAC() { ; }

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
      Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();


    };  // class Growth

  }  // namespace PAR

  class GrowthLawAC: public GrowthLawStatic
  {

  public:

    /// construct empty material object
    GrowthLawAC();
    explicit GrowthLawAC(MAT::PAR::GrowthLawAC* params);

    /// destructor
    virtual ~GrowthLawAC() { ; }



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
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate
    (
      double* theta,
      const double& thetaold,
      LINALG::Matrix<6,1>* dthetadC,
      MAT::Growth& matgrowth,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      const double& consttrig,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    ///calculate growth part of deformation gradient
    void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    double DensityDerivScale(const double theta);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac; }

    MAT::PAR::GrowthLawAC* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAC*>(MAT::GrowthLawStatic::Parameter()); }

    //! add parameter to parameter list
    void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) {dserror("You should not get here, "
        "function is only intended for GrowthLawIso so far, feel free to implement it for others if you need to");}

    //! Return whether material has a varying material density
    bool VaryingDensity() const {return true;}
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawACRadial : public GrowthLawAC
    {
    public:

      /// standard constructor
      GrowthLawACRadial(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawACRadial() { ; }

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawACRadial: public GrowthLawStatic
  {

  public:

    /// construct empty material object
    GrowthLawACRadial();
    explicit GrowthLawACRadial(MAT::PAR::GrowthLawAC* params);

    /// destructor
    virtual ~GrowthLawACRadial() { ; }

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
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate
    (
      double* theta,
      const double& thetaold,
      LINALG::Matrix<6,1>* dthetadC,
      MAT::Growth& matgrowth,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      const double& consttrig,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    ///calculate growth part of deformation gradient
    void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    double DensityDerivScale(const double theta);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac_radial; }

    MAT::PAR::GrowthLawACRadial* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawACRadial*>(MAT::GrowthLawStatic::Parameter()); }

    //! add parameter to parameter list
    void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) {dserror("You should not get here, "
        "function is only intended for GrowthLawIso so far, feel free to implement it for others if you need to");}

    //! Return whether material has a varying material density
    bool VaryingDensity() const {return true;}
  };


  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawACRadialRefConc : public GrowthLawAC
    {
    public:

      /// standard constructor
      GrowthLawACRadialRefConc(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawACRadialRefConc() { ; }

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawACRadialRefConc: public GrowthLawStatic
  {

  public:

    /// construct empty material object
    GrowthLawACRadialRefConc();
    explicit GrowthLawACRadialRefConc(MAT::PAR::GrowthLawAC* params);

    /// destructor
    virtual ~GrowthLawACRadialRefConc() { ; }

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
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate
    (
      double* theta,
      const double& thetaold,
      LINALG::Matrix<6,1>* dthetadC,
      MAT::Growth& matgrowth,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      const double& consttrig,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    ///calculate growth part of deformation gradient
    void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    double DensityDerivScale(const double theta);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac_radial_refconc; }

    MAT::PAR::GrowthLawACRadialRefConc* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawACRadialRefConc*>(MAT::GrowthLawStatic::Parameter()); }

    //! add parameter to parameter list
    void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) {dserror("You should not get here, "
        "function is only intended for GrowthLawIso so far, feel free to implement it for others if you need to");}

    //! Return whether material has a varying material density
    bool VaryingDensity() const {return true;}
  };



  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawConst : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawConst(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawConst() { ; }

      /// @name material parameters
      //@{
      //! enum for mapping between material parameter and entry in Parameter::matparams_
      enum matparamnames_ {thetarate, first=thetarate, last=thetarate};
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

      /// given names of optimizable material parameters upon request
      void OptParams(std::map<std::string, int>* pnames);

    };  // class Growth

  } // namespace PAR


  class GrowthLawConst : public GrowthLawStatic
  {

  public:

    /// construct empty material object
    GrowthLawConst();
    explicit GrowthLawConst(MAT::PAR::GrowthLawConst* params);

    /// destructor
    virtual ~GrowthLawConst() { ; }

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
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate
    (
      double* theta,
      const double& thetaold,
      LINALG::Matrix<6,1>* dthetadC,
      MAT::Growth& matgrowth,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      const double& consttrig,
      Teuchos::ParameterList& params,
      const int eleGID
    );

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
    void EvaluatePDeriv
    (
      double* theta,
      const double& thetaold,
      Teuchos::RCP<MAT::So3Material> matelastic,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    ///calculate growth part of deformation gradient
    void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    );

    ///density scaling of specific growth law
    double DensityScale(const double theta);

    ///density derivative scaling of specific growth law
    double DensityDerivScale(const double theta);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_const; }

    MAT::PAR::GrowthLawConst* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawConst*>(MAT::GrowthLawStatic::Parameter()); }

    //! add parameter to parameter list
    void AddParamsToParameterList(Teuchos::ParameterList& params,const double& theta) {dserror("You should not get here, "
        "function is only intended for GrowthLawIso so far, feel free to implement it for others if you need to");}

    //! Return whether material has a varying material density
    bool VaryingDensity() const {return true;}
  };




  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawIso : public Parameter
    {
    public:

      /// standard constructor
      GrowthLawIso(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawIso() { ; }

      /// @name material parameters
      //@{
      /// scalar that causes growth
      const int Sc1_;
      /// growth factor for scalar 1
      const double Sc1growthfac_;
      /// value of scalar 1 for that there is no growth
      const double Sc1refconc_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawIso: public GrowthLawStatic
  {

  public:

    //! construct empty material object
    GrowthLawIso();
    explicit GrowthLawIso(MAT::PAR::GrowthLawIso* params);

    //! destructor
    virtual ~GrowthLawIso() { ; }



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
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate
    (
      double* theta,
      const double& thetaold,
      LINALG::Matrix<6,1>* dthetadC,
      MAT::Growth& matgrowth,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      const double& consttrig,
      Teuchos::ParameterList& params,
      const int eleGID
    );

    //! calculate growth part of deformation gradient
    void CalcFg
    (
      const double& theta,
      const double& thetaold,
      const int& gp,
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<3,1>& refdir,
      const std::vector<LINALG::Matrix<3,1> >& curdir,
      const std::vector<LINALG::Matrix<3,3> >& histdefgrd,
      LINALG::Matrix<3,3>& F_g
    );

    //! density scaling of specific growth law
    double DensityScale(const double theta);

    //! density derivative scaling of specific growth law
    double DensityDerivScale(const double theta);

    //@}

    //! material type
    INPAR::MAT::MaterialType MaterialType() const {return INPAR::MAT::m_growth_iso;}

    MAT::PAR::GrowthLawIso* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawIso*>(MAT::GrowthLawStatic::Parameter()); }

    //! add parameter to parameter list
    void AddParamsToParameterList
    (
        Teuchos::ParameterList& params,
        const double& theta
    );

    //! Return whether material has a varying material density
    bool VaryingDensity() const {return false;}
  };



} // namespace MAT







#endif /* GROWTH_LAW_H_ */
