/*!----------------------------------------------------------------------
\file constraintmixture.H
\brief
This file contains routines for constraint mixture growth and remodeling.
example input line
example input line
MAT 1 MAT_ConstraintMixture DENS 0.001 MUE 1.0 PHIE 0.08 PREELA 1.0
K1 1.0 K2 1.0 PRECOLL 1.06 K1M 1.0 K2M 1.0 PHIM 1.0 PREMUS 1.0
SMAX 0.0 KAPPA 1.0E6 LIFETIME 5.0 HOMSTR 6.75E4 GROWTHFAC 0.5
STARTTIME 5.0 INTEGRATION Explicit TOL 1.0E-4 GROWTHFORCE Single
INITSTRETCH None DEGOPTION Cos

Here an approach for growth and remodeling of an artery is modeled.
For a detailed description see:
- Humphrey, J. & Rajagopal, K.: A constrained mixture model for arterial
  adaptations to a sustained step change in blood flow,
  Biomechanics and Modeling in Mechanobiology, 2003, 2, 109-126

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/


#ifndef CONSTRAINTMIXTURE_H_
#define CONSTRAINTMIXTURE_H_


#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class ConstraintMixture : public Parameter
    {
    public:

      /// standard constructor
      ConstraintMixture(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ConstraintMixture() { ; }

      /// @name material parameters
      //@{
      /// density
      const double density_;
      /// shear modulus
      const double mue_;
      /// Poisson's ratio
      const double nue_;
      /// mass fraction of elastin
      const double phielastin_;
      /// prestretch of elastin
      const double prestretchelastin_;
      /// parameter for linear fiber stiffness of collagen
      const double k1_;
      /// parameter for exponential fiber stiffness of collagen
      const double k2_;
      /// prestretch of collagen fibers
      const double prestretchcollagen_;
      /// stretch at which collagen fibers are damaged
      const double damagestretch_;
      /// parameter for linear fiber stiffness of smooth muscle
      const double k1muscle_;
      /// parameter for exponential fiber stiffness of smooth muscle
      const double k2muscle_;
      /// mass fraction of smooth muscle
      const double phimuscle_;
      /// prestretch of smooth muscle fibers
      const double prestretchmuscle_;
      /// maximal active stress
      const double Smax_;
      /// dilatation modulus
      const double kappa_;
      // basal rate of mass production
//      const double massprodbasal_;
      /// lifetime of collagen fibers
      const double lifetime_;
      /// growth factor for stress
      const double growthfactor_;
      /// homeostatic target value of scalar stress measure
      const double homstress_;
      /// growth factor for shear
      const double sheargrowthfactor_;
      /// homeostatic target value of inner radius
      const double homradius_;
      /// at this time turnover of collagen starts
      const double starttime_;
      /// time integration scheme (Explicit,Implicit)
      const std::string* integration_;
      /// tolerance for local Newton iteration
      const double abstol_;
      /// driving force of growth (Single,All,ElaCol)
      const std::string* growthforce_;
      /// form of elastin  degradation
      const std::string* elastindegrad_;
      /// how mass depends on driving force
      const std::string* massprodfunc_;
      /// how to set stretches in the beginning
      const std::string* initstretch_;
      /// number of timecurve for increase of prestretch in time
      const int timecurve_;
      /// which degradation function
      const std::string degoption_;
      /// tolerance for degradation
      const double degtol_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ConstraintMixture

  }  // namespace PAR

  class ConstraintMixtureHistoryType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "ConstraintMixtureHistoryType"; }

    static ConstraintMixtureHistoryType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static ConstraintMixtureHistoryType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for history of constraint mixture material
  class ConstraintMixtureHistory : public DRT::ParObject
  {
    friend class ConstraintMixture;
  public:
    /// construct empty history object
    ConstraintMixtureHistory(){;}
    /// destruct history object
    virtual ~ConstraintMixtureHistory() { ;}

    //! @name Packing and Unpacking
    virtual int UniqueParObjectId() const { return ConstraintMixtureHistoryType::Instance().UniqueParObjectId(); }
    virtual void Pack(DRT::PackBuffer& data) const;
    virtual void Unpack(const std::vector<char>& data);
    //@}

    /// Setup
    void Setup(int ngp, const double massprodbasal, bool expvar);

  private:

    /// @name Access to History
    //@{
    /// set time variables
    void SetTime(double deptime, double dt) {depositiontime_ = deptime; dt_ = dt;};
    /// get time variables
    void GetTime(double* deptime, double* dt) {*deptime = depositiontime_; *dt = dt_;};
    /// set stretches
    void SetStretches(int gp, LINALG::Matrix<4,1> stretches);
    /// get stretches
    void GetStretches(int gp, LINALG::Matrix<4,1>* stretches);
    /// set mass production rates
    void SetMass(int gp, LINALG::Matrix<4,1> massprod);
    /// set mass production rate of single fiber
    void SetMass(int gp, double massprod, int idfiber);
    /// get mass production rates
    void GetMass(int gp, LINALG::Matrix<4,1>* massprod);
    /// set vardegrad
    void SetVarDegrad(int gp, int idfiber, double vardegrad);
    /// get vardegrad
    void GetVarDegrad(int gp, int idfiber, double* vardegrad);
    /// return number of gausspoints
    int NumGP() const { return numgp_;}
    //@}

    /// deposition time of collagen fibers
    double depositiontime_;
    /// time step at this time
    double dt_;
    /// number of gausspoints
    int numgp_;
    /// variable degradation of collagen?
    bool expvar_;
    /// stretch of first collagen fiber family
    Teuchos::RCP<std::vector<double> > collagenstretch1_;
    /// stretch of second collagen fiber family
    Teuchos::RCP<std::vector<double> > collagenstretch2_;
    /// stretch of third collagen fiber family
    Teuchos::RCP<std::vector<double> > collagenstretch3_;
    /// stretch of fourth collagen fiber family
    Teuchos::RCP<std::vector<double> > collagenstretch4_;
    /// mass production rate of first fiber family
    Teuchos::RCP<std::vector<double> > massprod1_;
    /// mass production rate of second fiber family
    Teuchos::RCP<std::vector<double> > massprod2_;
    /// mass production rate of third fiber family
    Teuchos::RCP<std::vector<double> > massprod3_;
    /// mass production rate of fourth fiber family
    Teuchos::RCP<std::vector<double> > massprod4_;
    /// mass production rate of first fiber family
    Teuchos::RCP<std::vector<double> > vardegrad1_;
    /// mass production rate of second fiber family
    Teuchos::RCP<std::vector<double> > vardegrad2_;
    /// mass production rate of third fiber family
    Teuchos::RCP<std::vector<double> > vardegrad3_;
    /// mass production rate of fourth fiber family
    Teuchos::RCP<std::vector<double> > vardegrad4_;

  };

  class ConstraintMixtureType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "ConstraintMixtureType"; }

    static ConstraintMixtureType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static ConstraintMixtureType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for constraint mixture material
  class ConstraintMixture : public So3Material
  {
  public:

    /// construct empty material object
    ConstraintMixture();

    /// construct the material object given material parameters
    explicit ConstraintMixture(MAT::PAR::ConstraintMixture* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ConstraintMixtureType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
      History data is unpacked in restart.

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_constraintmixture; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ConstraintMixture( *this ) ); }

    /// Setup
    virtual void Setup
    (
      int numgp,  ///< number of Gauss points
      DRT::INPUT::LineDefinition* linedef ///< definition of element line
    );

    /// SetupHistory
    virtual void ResetAll(const int numgp);

    /// Update
    virtual void Update();

    /// Reset time step
    virtual void ResetStep();

    /// Evaluate material
    virtual void Evaluate
    (
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
      Teuchos::ParameterList& params,
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat
    );

    /// Return density
    virtual double Density() const {return params_->density_;}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return variables for visualization
    LINALG::Matrix<3,1> GetVis(int gp) const {return vismassstress_->at(gp);}
    /// Return actual mass density in reference configuration
    double GetMassDensity(int gp) const {return refmassdens_->at(gp);}
    /// Return actual mass density in reference configuration
    LINALG::Matrix<3,1> GetMassDensityCollagen(int gp) const {return visrefmassdens_->at(gp);}
    /// Return prestretch of collagen fibers
    LINALG::Matrix<3,1> GetPrestretch(int gp) const
    {
      LINALG::Matrix<3,1> visprestretch(true);
      visprestretch(0) = localprestretch_->at(gp)(0);
      visprestretch(1) = localprestretch_->at(gp)(1);
      visprestretch(2) = localprestretch_->at(gp)(2);
      return visprestretch;
    }
    /// Return prestretch of collagen fibers
    LINALG::Matrix<3,1> GetHomstress(int gp) const
    {
      LINALG::Matrix<3,1> visprestretch(true);
      visprestretch(0) = localhomstress_->at(gp)(0);
      visprestretch(1) = localhomstress_->at(gp)(1);
      visprestretch(2) = localhomstress_->at(gp)(2);
      return visprestretch;
    }
    /// Return circumferential fiber direction
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > Geta1() const {return a1_;}
    /// Return axial fiber direction
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > Geta2() const {return a2_;}

    /// evaluate fiber directions from locsys, pull back
    void EvaluateFiberVecs
    (const int gp, const LINALG::Matrix<3,3>& locsys, const LINALG::Matrix<3,3>& defgrd);

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

  private:

    /// my material parameters
    MAT::PAR::ConstraintMixture* params_;

    /// temporary for visualization
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > vismassstress_;
    /// actual mass density in reference configuration
    Teuchos::RCP<std::vector<double> > refmassdens_;
    /// actual mass density in reference configuration for collagen fibers
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > visrefmassdens_;
    /// basal rate of mass production
    double massprodbasal_;

    /// first fiber vector per gp (reference), circumferential
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > a1_;
    /// second fiber vector per gp (reference), axial
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > a2_;
    /// third fiber vector per gp (reference), diagonal
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > a3_;
    /// fourth fiber vector per gp (reference), diagonal
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > a4_;
    /// homeostatic prestretch of collagen fibers
    Teuchos::RCP<std::vector<LINALG::Matrix<4,1> > > localprestretch_;
    /// homeostatic stress for growth
    Teuchos::RCP<std::vector<LINALG::Matrix<4,1> > > localhomstress_;
    /// homeostatic radius
    double homradius_;
    /// list of fibers which have been overstretched and are deleted at the end of the time step (gp, idpast, idfiber)
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > deletemass_;
    /// indicates if material is initialized
    bool isinit_;

    /// history
    Teuchos::RCP<std::vector<ConstraintMixtureHistory> > history_;

    /// compute stress and cmat
    void EvaluateStress
    (
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      const int gp,  ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      const int firstiter,  ///< iteration index, different for explicit and implicit integration
      double time,  ///< time
      double elastin_survival   ///< amount of elastin which is still there
    );

    /// elastic response for one collagen fiber family
    void EvaluateFiberFamily
    (
      const LINALG::Matrix<NUM_STRESS_3D,1> C,  ///< Cauchy-Green
      const int gp,  ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      LINALG::Matrix<3,1> a,  ///< fiber vector
      double* currmassdens,  ///< current massdensity
      const int firstiter,  ///< iteration index, different for explicit and implicit integration
      double time,  ///< time
      const int idfiber  ///< number of fiber family 0,1,2,3
     );

    /// elastic response for one collagen fiber
    void EvaluateSingleFiberScalars
    (
      double I4,  ///< fourth invariant multiplied with (prestretch / stretch at deposition time)^2
      double& fac_cmat,  ///< scalar factor for material stiffness matrix
      double& fac_stress  ///< scalar factor for 2nd PK-stress
    );

    /// elastic response for Elastin
    void EvaluateElastin
    (
      const LINALG::Matrix<NUM_STRESS_3D,1> C,  ///< Cauchy-Green
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      double time,  ///< time
      double* currmassdens,  ///< current massdensity
      double elastin_survival  ///< amount of elastin which is still there
    );

    /// elastic response for smooth muscle cells
    void EvaluateMuscle
    (
      const LINALG::Matrix<NUM_STRESS_3D,1> C,  ///< Cauchy-Green
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      const int gp, ///< current Gauss point
      double* currmassdens  ///< current massdensity
    );

    /// volumetric part
    void EvaluateVolumetric
    (
      const LINALG::Matrix<NUM_STRESS_3D,1> C,  ///< Cauchy-Green
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      double currMassDens,  ///< current mass density
      double refMassDens  ///< initial mass density
    );

    /// computes mass production rates for all fiber families
    void MassProduction
    (
      const int gp,  ///< current Gauss point
      LINALG::Matrix<3,3> defgrd,  ///< deformation gradient
      LINALG::Matrix<NUM_STRESS_3D,1> S,  ///< 2nd PK-stress
      LINALG::Matrix<4,1>* massstress,  ///< growth stress measure
      double inner_radius,  ///< inner radius
      LINALG::Matrix<4,1>* massprodcomp  ///< mass production rate
    );

    /// computes mass production rate for one fiber family
    void MassProductionSingleFiber
    (
      const int gp,  ///< current Gauss point
      LINALG::Matrix<3,3> defgrd,  ///< deformation gradient
      LINALG::Matrix<NUM_STRESS_3D,1> S,  ///< 2nd PK-stress
      double* massstress,  ///< growth stress measure
      double inner_radius,  ///< inner radius
      double* massprodcomp,  ///< mass production rate
      LINALG::Matrix<3,1> a,  ///< fiber vector
      const int idfiber  ///< number of fiber family 0,1,2,3
    );

    /// function for massproduction
    void MassFunction
    (
      double growthfac,  ///< growth factor K
      double delta,  ///< difference in stress or wall shear stress
      double mmax,  ///< maximal massproduction fac
      double& massfac  ///< factor
    );

    /// returns how much collagen has been degraded
    void Degradation(double t, double& degr);

    /// function of elastin degradation (initial)
    void ElastinDegradation
    (
      LINALG::Matrix<1,3> coord,   ///< gp coordinate in reference configuration
      double& degr   ///< amount of elastin which is still left
    );

    /// compute stress and cmat for implicit integration with whole stress as driving force
    void EvaluateImplicitAll
    (
      LINALG::Matrix<3,3> defgrd,  ///< deformation gradient
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      const int gp,  ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      double dt,  ///< delta time
      double time,  ///< time
      LINALG::Matrix<4,1> massprodcomp,  ///< mass production rate
      LINALG::Matrix<4,1> massstress  ///< growth stress measure
    );

    /// compute stress and cmat for implicit integration with fiber stress as driving force
    void EvaluateImplicitSingle
    (
      LINALG::Matrix<3,3> defgrd,  ///< deformation gradient
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      const int gp,  ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      double dt,  ///< delta time
      double time  ///< time
    );

    /// derivative of stress with respect to massproduction of a single fiber family
    void GradStressDMass
    (
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      LINALG::Matrix<NUM_STRESS_3D,1>* derivative,  ///< result
      LINALG::Matrix<NUM_STRESS_3D,1> Cinv,  ///< inverse cauchy green strain
      LINALG::Matrix<3,1> a,  ///< fiber vector
      double stretch,  ///< prestretch / stretch at deposition time
      double J,  ///< determinant of F
      double dt,  ///< delta time
      bool option  ///< which driving force
    );

    /// derivative of massproduction of a single fiber family with respect to stress
    void GradMassDStress
    (
      LINALG::Matrix<NUM_STRESS_3D,1>* derivative,  ///< result
      LINALG::Matrix<3,3> defgrd,  ///< deformation gradient
      LINALG::Matrix<3,3> Smatrix,  ///< 2nd PK-stress in matrix notation
      LINALG::Matrix<3,1> a,  ///< fiber vector
      double J,  ///< determinant of F
      double massstress,  ///< growth stress measure
      double actcollstretch  ///< actual collagen stretch
    );

    /// derivative of massproduction of a single fiber family with respect to stretch
    void GradMassDStretch
    (
      LINALG::Matrix<NUM_STRESS_3D,1>* derivative,  ///< result
      LINALG::Matrix<3,3> defgrd,  ///< deformation gradient
      LINALG::Matrix<3,3> Smatrix,  ///< 2nd PK-stress in matrix notation
      LINALG::Matrix<NUM_STRESS_3D,1> Cinv,
      LINALG::Matrix<3,1> a,  ///< fiber vector
      double J,  ///< determinant of F
      double massstress,  ///< growth stress measure
      double actcollstretch,  ///< actual collagen stretch
      double dt  ///< delta time
    );

  };

  /// Debug output to gmsh-file
  /* this needs to be copied to STR::TimInt::OutputStep() to enable debug output
  {
    discret_->SetState("displacement",Dis());
    MAT::ConstraintMixtureOutputToGmsh(discret_, GetStep(), 1);
  }
  don't forget to include constraintmixture.H */
  void ConstraintMixtureOutputToGmsh
  (
    const Teuchos::RCP<DRT::Discretization> dis, ///< discretization with displacements
    const int timestep, ///< index of timestep
    const int iter ///< iteration index of newton iteration
  );

}  // namespace MAT

#endif /* CONSTRAINTMIXTURE_H_ */
