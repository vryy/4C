/*!----------------------------------------------------------------------
\file growth_ip.H
\brief
This file contains routines for an integration point based growth law
example input line
MAT 1 MAT_GROWTH DENS 1.0 IDMATELASTIC 2 STARTTIME 0.2 ENDTIME 100.0 KPLUS 0.5 MPLUS 4.0 KMINUS 0.25 MMINUS 5.0

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


#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  class GrowthLaw;

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
      /// elastic material number
      const int idmatelastic_;
      /// elastic material number
      const int idgrowthlaw_;
      /// start growth after starttime
      const double starttime_;
      /// stop growth after endtime
      const double endtime_;
      /// tolerance for local Newton iteration
      const double abstol_;
      /// growth law
      Teuchos::RCP<MAT::GrowthLaw> growthlaw_;
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

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static GrowthType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for growth material
  class Growth : public So3Material
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Growth( *this ) ); }

    /// Setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    virtual void ResetAll(int numgp);

    /// Update
    virtual void Update();

    /// Reset time step
    virtual void ResetStep();

    /// Evaluate material
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,  ///<green lagrange strain
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* stress,  ///< 2nd PK-stress
                          LINALG::Matrix<6,6>* cmat, ///< material stiffness matrix
                          const int eleGID );

    /// Evaluate material law
    virtual void EvaluateNonLinMass(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* linmass_disp,
                          LINALG::Matrix<6,1>* linmass_vel,
                          const int eleGID);

    /// Return density
    virtual double Density() const {dserror("growth material needs gauss point data for density!"); return 0.0;}
    virtual double Density(int gp) const
    {
      double density = matelastic_->Density();
      double theta_gp = (*theta_)[gp];
      return density * theta_gp * theta_gp * theta_gp;
    };

    /// Return whether material has a varying material density
    virtual bool VaryingDensity() const
    {
      return true;
    }

    /// Return theta
    Teuchos::RCP<std::vector<double> > Gettheta() const {return theta_;}

    /// Return mandel stress
    Teuchos::RCP<std::vector<double> > Getmandel() const {return mandel_;}

    /// Return elastic material
    //Teuchos::RCP<MAT::Material> Matelastic() const { return matelastic_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

  protected:

    /// growth stretch
    Teuchos::RCP<std::vector<double> > theta_;
    /// growth stretch old time step
    Teuchos::RCP<std::vector<double> > thetaold_;
    /// mandel stress
    Teuchos::RCP<std::vector<double> > mandel_;
    /// elastic material
    Teuchos::RCP<MAT::So3Material> matelastic_;

    /// indicates if material is initialized
    bool isinit_;

    /// Evaluate growth function
    virtual void EvaluateGrowthFunction
    (
        double & growthfunc,
        double traceM,
        double theta
    );

    /// Evaluate derivative of growth function w.r.t. growth factor
    virtual void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    );

    /// Evaluate derivative of growth function w.r.t. Green-Lagrange-Strain
    virtual void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    );

  private:

    /// my material parameters
    MAT::PAR::Growth* params_;

  };

  /// Debug output to gmsh-file
  /* this needs to be copied to STR::TimInt::OutputStep() to enable debug output
  {
    discret_->SetState("displacement",Dis());
    MAT::GrowthOutputToGmsh(discret_, GetStep(), 1);
  }
  don't forget to include growth_ip.H */
  void GrowthOutputToGmsh
  (
    const Teuchos::RCP<DRT::Discretization> dis, ///< discretization with displacements
    const int timestep, ///< index of timestep
    const int iter ///< iteration index of newton iteration
  );

}  // namespace MAT

#endif /*GROWTH_IP_H_*/
