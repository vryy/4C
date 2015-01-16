/*----------------------------------------------------------------------*/
/*!
 \file growth_scd.H

 \brief
This file contains routines for an integration point based and scalar dependend volumetric growth law.
It is derived from the growth laws implemented in growth_ip.cpp and additional adds the scalar dependency, e.g. nutrients.

 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#ifndef GROWTH_SCD_H_
#define GROWTH_SCD_H_

/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "growth_ip.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class GrowthScd : public Growth
    {
    public:

      /// standard constructor
      GrowthScd(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthScd() { ; }

      /// @name material parameters
      //@{
      /// substrate uptake rate coefficient
      const double rearate_;
      /// saturation coefficient for concentration dependent growth law
      const double satcoeff_;
      //@}

      // Coupling between stress dependent growth and reaction dependent growth
      // it could be that the stress dependent growth has to be scaled with the reaaction term (ScaleConc)
      // or that the reaction dependent growth is reduced by the presence of stresses (StressRed)
      const std::string* growthcoupl_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
    };  // class ScdGrowth

  }  // namespace PAR

  class GrowthScdType : public GrowthMandelType
  {
  public:

    std::string Name() const { return "GrowthScdType"; }

    static GrowthScdType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static GrowthScdType instance_;
  }; //class ScdGrowthType

  /*----------------------------------------------------------------------*/
  /// Wrapper for growth material

  /// See Bachelor thesis Jasper Rieser (2013, supervisor: Anh-Thu Vuong)
  class GrowthScd : public GrowthMandel                                                    //

  {
  public:

    /// construct empty material object
    GrowthScd();

    /// construct the material object given material parameters
    explicit GrowthScd(MAT::PAR::GrowthScd* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return GrowthScdType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the std::vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

      \param data (in/out): char std::vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;
    /*!
      \brief Unpack data from a char std::vector into this class

      The std::vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
      History data is unpacked in restart.

      \param data (in) : std::vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_volumetric_scd; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthScd( *this ) ); }

    /// Setup
    void Setup
    (
      int numgp,  ///< number of Gauss points
      DRT::INPUT::LineDefinition* linedef ///< definition of element line
    );

    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    void ResetAll(int numgp);

    /// Evaluate material
    void Evaluate
    (
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        LINALG::Matrix<6,1>* stress,
        LINALG::Matrix<6,6>* cmat,
        const int eleGID
     );

    /// Return density
    //virtual double Density() const {return params_->density_;}

    /// Return theta
    Teuchos::RCP<std::vector<double> > Gettheta() const {return theta_;}

    /// Return theta at gauss-point
    double Gettheta_atgp(int gp) const { return (*theta_)[gp];}

//    // Return thetaold
//    double Getthetaold_atgp(int gp) const {return (*thetaold_)[gp];}

    // Return thetaold at gauss-point
    double Getthetaold_atgp(int gp) const {return (*thetaold_)[gp];}

    /// Return time derivate of theta at gauss-point
    double Getdtheta_atgp(int gp) const {return (*dtheta_)[gp];}

    /// Return determinate of F_e at gauss-point
    double GetdetFe_atgp(int gp) const {return (*detFe_)[gp];}

    /// Return mandel stress
    Teuchos::RCP<std::vector<double> > Getmandel() const {return mandel_;}

    /// Return elastic material
    Teuchos::RCP<MAT::Material> Matelastic() const { return matelastic_; }

    /// Return quick accessible material parameter data
    MAT::PAR::GrowthScd* Parameter() const { return paramsScd_; }

    ///Return concentration
    double GetConcentration() const {return concentration_;}

  protected:
    /// Evaluate growth function (overwrites implementation in stress based growth law)
    void EvaluateGrowthFunction
    (
        double & growthfunc,
        double traceM,
        double theta
    );

    /// Evaluate derivative of growth function w.r.t. growth factor (overwrites implementation in stress based growth law)
    void EvaluateGrowthFunctionDerivTheta
    (
        double & dgrowthfunctheta,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
    );

    /// Evaluate derivative of growth function w.r.t. Green-Lagrange-Strain (overwrites implementation in stress based growth law)
    void EvaluateGrowthFunctionDerivC
    (
        LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
        double traceM,
        double theta,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
        const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
        const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
    );

  protected:

    /// determinante of elastic deformation gradient F_e for each gauss point
    Teuchos::RCP<std::vector<double> > detFe_;
    /// time derivative of growth factor theta for each gauss point
    Teuchos::RCP<std::vector<double> > dtheta_;
    /// concentration (e.g. nutrients)
    double concentration_;

    /// growth function value evaluated by the stress law
    double stressgrowthfunc;

  private:

    /// my material parameters
    MAT::PAR::GrowthScd* paramsScd_;

  };//class GrowthScd



  class GrowthScdACType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "GrowthScdACType"; }

    static GrowthScdACType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static GrowthScdACType instance_;
  }; //class ScdGrowthType

  /*----------------------------------------------------------------------*/
  /// Wrapper for growth material
  class GrowthScdAC : public GrowthBasic                                                    //
  {
  public:

    /// construct empty material object
    GrowthScdAC();                                                                      //

    /// construct the material object given material parameters
    explicit GrowthScdAC(MAT::PAR::GrowthScd* params);                                  //

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const { return GrowthScdACType::Instance().UniqueParObjectId(); }


    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    void ResetAll(int numgp);

    /*!
      \brief Pack this class so it can be communicated

      Resizes the std::vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

      \param data (in/out): char std::vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const;
    /*!
      \brief Unpack data from a char std::vector into this class

      The std::vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
      History data is unpacked in restart.

      \param data (in) : std::vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_volumetric_scd; }


    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthScdAC( *this ) ); }

    /// Setup
    void Setup
    (
      int numgp,  ///< number of Gauss points
      DRT::INPUT::LineDefinition* linedef ///< definition of element line
    );

    /// Evaluate material
    void Evaluate
    (
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        LINALG::Matrix<6,1>* stress,
        LINALG::Matrix<6,6>* cmat,
        const int eleGID
     );

    /// Return elastic material
    Teuchos::RCP<MAT::Material> Matelastic() const { return matelastic_; }

    /// Return quick accessible material parameter data
    MAT::PAR::GrowthScd* Parameter() const { return paramsScdAC_; }

  protected:

    ///Calculate the volumetric growth parameter
    double CalculateTheta( const double J );

    ///Calculate the volumetric growth derived w.r.t. chauy-green strains
    void CalculateThetaDerivC( LINALG::Matrix<3,3>& dThetadC ,
                               const LINALG::Matrix<3,3>& C, ///< cauchy-green strain
                               const double J ///< det(F)
                               );

    /// pointer to mean concentrations at present timestep
    Teuchos::RCP<std::vector<double> > concentrations_;

  private:

    /// my material parameters
    MAT::PAR::GrowthScd* paramsScdAC_;

  };//class GrowthScdAC

}  // namespace MAT


#endif /* GROWTH_SCD_H_ */
