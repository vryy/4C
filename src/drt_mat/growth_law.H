/*----------------------------------------------------------------------*/
/*!
\file growth_law.H

\brief growth law for the volumetric growth materials. This file contains routines needed for the calculation of the volumetric growth
parameter theta. These so called growthlaws are to be used for the isovolumetric
split growth materials MAT::Growth. Upon request a growthlaw delivers a growth factor \f$\vartheta\f$ and its derivative
wrt. \f$\frac{\partial \vartheta}{\partial C}\f$.

\level 3

 <pre>
   \maintainer Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
*/

/*----------------------------------------------------------------------*/


#ifndef GROWTH_LAW_H
#define GROWTH_LAW_H

#include <Teuchos_RCP.hpp>
#include "matpar_parameter.H"
#include "material.H"
#include "so3_material.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "growth_ip.H"


namespace MAT
{
  /*----------------------------------------------------------------------*/
  /*! \class GrowthLaw
      \brief GrowthLaw base class

      It provides the interfaces Evaluate(theta,dthetadC,...) to compute
      the growth factor \f$\vartheta\f$ and its derivative
      \f$\frac{\partial \vartheta}{\partial C}\f$ and
      EvaluateNonLinMass(...,linmass_disp,...) to account for the dynamic change
      of density.

      All growthlaws can be modified by a function depending on scalar concentrations
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
    virtual void EvaluatePDeriv(double* theta,
        const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        const int eleGID ) = 0;

    //@}


    //! @name Pre-Evaluation methods
    //@{

    /*! \brief set scalar factor
     * set one double as scalar factor
     */
    void SetFactor(double& fac) {cfac_.clear(); cfac_.push_back(fac);}

    /*! \brief set scalar factor
     * set a vector of doubles as factors
     */
    void SetFactor(std::vector<double> fac) {cfac_=fac;};

    //@}

    //! @name Query methods
    //@{

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    //! return material parameters
    MAT::PAR::Parameter* Parameter() {return params_;};

    //! access factor vector of scalar factors
    virtual std::vector<double> GetCfac() {return cfac_;}

    //@}

  private:

    //! material parameters
    MAT::PAR::Parameter* params_;

    //! possible scalars to be processed eventually by the specific laws
    std::vector<double> cfac_;

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
    virtual void Evaluate(double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        const int eleGID);

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
    virtual void EvaluatePDeriv(double* theta,
        const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        const int eleGID );

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
        double traceM,
        double theta
    ) =0;

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
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) =0;

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
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) =0;

    //@}

    //! return the material parameters
    MAT::PAR::GrowthLawDyn* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawDyn*>(MAT::GrowthLaw::Parameter()); };

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
    virtual void Evaluate(double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
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
    virtual void EvaluatePDeriv(double* theta,
        const double& thetaold,
        Teuchos::RCP<MAT::So3Material> matelastic,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        const int eleGID );

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

  protected:

  };

  namespace PAR
  {
  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawLinear
       \brief Common parameters for linear growth law

       \author kehl
       \date 6/2015
   */
    class GrowthLawLinear : public GrowthLawDyn
    {
    public:

      /// standard constructor
      GrowthLawLinear(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawLinear() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class GrowthLawLinear
     \brief GrowthLaw based on a linear growth function

     The rhs-function \f$ f_{rhs}\f$ (see \sa GrowthLawDyn) is implemented from
     'Himpel,G. Kuhl, E. Menzel, A. & Steinmann, P. Computational modelling
     of isotropic multiplicative growth, Computer Modeling in Engineering
     and Sciences, 2005, 8, 119-134'. See also Diss. Susanna Tinkl.

     \author kehl
     \date 6/2015
   */
  class GrowthLawLinear: public GrowthLawDyn
  {

  public:

    /// construct empty material object
    GrowthLawLinear();
    explicit GrowthLawLinear(MAT::PAR::GrowthLawLinear* params);

    /// destructor
    virtual ~GrowthLawLinear() { ; }

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
        double traceM,
        double theta
    ) ;

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
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) ;

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
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) ;

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_linear; };

    //! material parameters
    MAT::PAR::GrowthLawLinear* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawLinear*>(MAT::GrowthLawDyn::Parameter()); };

  };


  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawExp : public GrowthLawDyn
    {
    public:

      /// standard constructor
      GrowthLawExp(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawExp() { ; }

      /// @name material parameters
      //@{
      /// reference value for trace of mandel stress
      const double mandel_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  }  // namespace PAR

  class GrowthLawExp: public GrowthLawDyn
  {

  public:

    /// construct empty material object
    GrowthLawExp();
    explicit GrowthLawExp(MAT::PAR::GrowthLawExp* params);

    /// destructor
    virtual ~GrowthLawExp() { ; }

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
        double traceM,
        double theta
    ) ;

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
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    ) ;

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
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    ) ;

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_exponential; };

    MAT::PAR::GrowthLawExp* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawExp*>(MAT::GrowthLawDyn::Parameter()); };

  };

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthLawBiofilm : public GrowthLawDyn
    {
    public:

      /// standard constructor
      GrowthLawBiofilm(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthLawBiofilm() { ; }

      /// @name material parameters
      //@{
      /// reference value for trace of mandel stress
      const double mandel_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

    };  // class Growth

  } // namespace PAR

  class GrowthLawBiofilm : public GrowthLawDyn
    {

    public:

      /// construct empty material object
      GrowthLawBiofilm();
      explicit GrowthLawBiofilm(MAT::PAR::GrowthLawBiofilm* params);

      /// destructor
      virtual ~GrowthLawBiofilm() { ; }

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
          double traceM,
          double theta
      ) ;

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
          const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
          const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
      ) ;

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
          LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
          double traceM,
          double theta,
          const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
          const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
          const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
      ) ;

      //@}

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_biofilm; };

      MAT::PAR::GrowthLawBiofilm* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawBiofilm*>(MAT::GrowthLawDyn::Parameter()); };

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

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
      virtual void Evaluate(double* theta,
          const double& thetaold,
          LINALG::Matrix<6,1>* dthetadC,
          MAT::Growth& matgrowth,
          const LINALG::Matrix<3,3>* defgrd,
          const LINALG::Matrix<6,1>* glstrain,
          Teuchos::ParameterList& params,
          const int eleGID );

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
      virtual void EvaluatePDeriv(double* theta,
          const double& thetaold,
          Teuchos::RCP<MAT::So3Material> matelastic,
          const LINALG::Matrix<3,3>* defgrd,
          const LINALG::Matrix<6,1>* glstrain,
          Teuchos::ParameterList& params,
          const int eleGID );

      //@}

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_const; };

      MAT::PAR::GrowthLawConst* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawConst*>(MAT::GrowthLawStatic::Parameter()); };

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// create growth law instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

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
    virtual void Evaluate(double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        const int eleGID );

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac; };

    MAT::PAR::GrowthLawAC* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawAC*>(MAT::GrowthLawStatic::Parameter()); };

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
      virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

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
    virtual void Evaluate(double* theta,
        const double& thetaold,
        LINALG::Matrix<6,1>* dthetadC,
        MAT::Growth& matgrowth,
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        const int eleGID );

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac_radial; };

    MAT::PAR::GrowthLawACRadial* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawACRadial*>(MAT::GrowthLawStatic::Parameter()); };

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
        virtual Teuchos::RCP<MAT::GrowthLaw> CreateGrowthLaw();

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
      virtual void Evaluate(double* theta,
          const double& thetaold,
          LINALG::Matrix<6,1>* dthetadC,
          MAT::Growth& matgrowth,
          const LINALG::Matrix<3,3>* defgrd,
          const LINALG::Matrix<6,1>* glstrain,
          Teuchos::ParameterList& params,
          const int eleGID );

      //@}

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_ac_radial_refconc; };

      MAT::PAR::GrowthLawACRadialRefConc* Parameter() {return dynamic_cast<MAT::PAR::GrowthLawACRadialRefConc*>(MAT::GrowthLawStatic::Parameter()); };

    };
} // namespace MAT


#endif /* GROWTH_LAW_H_ */
