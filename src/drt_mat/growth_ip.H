/*!----------------------------------------------------------------------
\file growth_ip.H
\brief
This file contains routines for an integration point based growth law
example input line
MAT 1 MAT_GROWTH DENS 1.0 MATELASTIC 2 STARTTIME 0.2 KPLUS 0.5 MPLUS 4.0 KMINUS 0.25 MMINUS 5.0

Here a kinematic integration point based approach of growth is modeled.
For a detailed description see:
- Lubarda, V. & Hoger, A., On the mechanics of solids with a growing mass,
  International Journal of Solids and Structures, 2002, 39, 4627-4664
- Himpel, G.; Kuhl, E.; Menzel, A. & Steinmann, P., Computational modelling
  of isotropic multiplicative growth, Computer Modeling in Engineering
  and Sciences, 2005, 8, 119-134

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifndef GROWTH_IP_H_
#define GROWTH_IP_H_

#ifdef CCADISCRET

#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class Growth : public Parameter
    {
    public:

      /// standard constructor
      Growth(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Growth() { ; }

      /// @name material parameters
      //@{
      // density
      const double density_;
      // elastic material number
      const int matelastic_;
      // start growth after starttime
      const double starttime_;
      // growth law parameter kthetaplus
      const double kthetaplus_;
      // growth law parameter mthetaplus
      const double mthetaplus_;
      // growth law parameter kthetaminus
      const double kthetaminus_;
      // growth law parameter mthetaminus
      const double mthetaminus_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Growth

  }  // namespace PAR

  class GrowthType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "GrowthType"; }

    static GrowthType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const vector<char> & data );

  private:

    static GrowthType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for growth material
  class Growth : public Material
  {
  public:

    /// construct empty material object
    Growth();

    /// construct the material object given material parameters
    explicit Growth(MAT::PAR::Growth* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return GrowthType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(vector<char>& data) const;

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Growth( *this ) ); }

    /// Setup
    void Setup
    (
      const int numgp  ///< number of Gauss points
    );

    /// Update
    void Update();

    /// Evaluate material
    void Evaluate
    (
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      const int gp,   ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,  ///< 2nd PK-stress
      const double dt, ///< delta time
      const double time ///< time
    );

    /// Return density
    double Density() const {return params_->density_;};

    /// return theta
    RCP<vector<double> > Gettheta() const {return theta_;};

    /// return mandel stress
    RCP<vector<double> > Getmandel() const {return mandel_;};

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::Growth* params_;

    /// growth parameter
    RCP<vector<double> > theta_;
    /// growth parameter old time step
    RCP<vector<double> > thetaold_;
    /// mandel stress
    RCP<vector<double> > mandel_;

    bool isinit_;  ///< indicates if material is initialized

    /// Evaluate elastic material
    void EvaluateElastic
    (
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      const int gp,   ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress  ///< 2nd PK-stress
    );

    /// Evaluate growth law
    void EvaluateGrowthLaw
    (
      double * ktheta, ///< factor in growth law
      double * dktheta, ///< derivative of ktheta
      double traceM, ///< trace of Mandel stress
      double theta, ///< growth parameter
      double mandelcrit ///< stress at which growth begins
    );

  };
}  // namespace MAT

#endif /*CCADISCRET*/
#endif /*GROWTH_IP_H_*/
