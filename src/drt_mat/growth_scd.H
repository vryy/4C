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
  }; //class ScdGrowthACType

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
    virtual int UniqueParObjectId() const { return GrowthScdACType::Instance().UniqueParObjectId(); }


    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    virtual void ResetAll(int numgp);

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
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthScdAC( *this ) ); }

    /// Setup
    virtual void Setup
    (
      int numgp,  ///< number of Gauss points
      DRT::INPUT::LineDefinition* linedef ///< definition of element line
    );

    /// Evaluate material
    virtual void Evaluate
    (
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<6,1>* glstrain,
        Teuchos::ParameterList& params,
        LINALG::Matrix<6,1>* stress,
        LINALG::Matrix<6,6>* cmat,
        const int eleGID
     );

    /// Return quick accessible material parameter data
    virtual MAT::PAR::GrowthScd* Parameter() const { return paramsScdAC_; }

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

  ///type of vector/tensor in voigt notation
  enum VoigtType
  {
    voigt_strain, ///< strain vector/tensor in voigt notation
    voigt_stress ///< stress vector/tensor in voigt notation
  };

  class GrowthScdACRadialType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "GrowthScdACRadialType"; }

    static GrowthScdACRadialType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static GrowthScdACRadialType instance_;
  }; //class ScdGrowthACRadialType

  /*----------------------------------------------------------------------*/
  /// Wrapper for growth material
  class GrowthScdACRadial : public GrowthScdAC                                                    //
  {
  public:

    /// construct empty material object
    GrowthScdACRadial();                                                                      //

    /// construct the material object given material parameters
    explicit GrowthScdACRadial(MAT::PAR::GrowthScd* params);                                  //

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const { return GrowthScdACRadialType::Instance().UniqueParObjectId(); }

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

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_volumetric_scd; }


    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthScdACRadial( *this ) ); }

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

  protected:

    ///Function which reads in the given fiber value due to the FIBER1 nomenclature
    void ReadFiber(
        DRT::INPUT::LineDefinition* linedef,
        std::string specifier,
        LINALG::Matrix<3,1> &fiber_vector
    );

//    LINALG::Matrix<3,3> CalcDerivMatrixInverse(const LINALG::Matrix<3,3>& Matrix,
//                        const LINALG::Matrix<3,3>& );
//
//    double CalcDerivMatrixInverseij(const LINALG::Matrix<3,3>& Minv,
//                        const LINALG::Matrix<3,3>& Mderiv,
//                        const double i,
//                        const double j);

    ///transform vector in voigt notation into symmetric 2Tensor (in matrix notation)
    void VectorToMatrix(
        LINALG::Matrix<3,3>& Matrix,
        const LINALG::Matrix<6,1>& Vector,
        const MAT::VoigtType Type);

    ///transform symmetric 2Tensor(in matrix notation) into voigt notation ( e.g. vector notation)
    void MatrixToVector(
        const LINALG::Matrix<3,3>& Matrix,
        LINALG::Matrix<6,1>& Vector,
        const MAT::VoigtType Type);

    /// pull back of a symmetric elastic 4th order tensor (in matrix/voigt notation)
    /// via the 2th order deformation gradient (also in matrix notation)
    LINALG::Matrix<6,6> PullBack4Tensor(const LINALG::Matrix<3,3>& defgr,
                            const LINALG::Matrix<6,6>& Cmat);

    /// pull back the ijkl-th entry of a symmetric elastic 4th order tensor (in matrix/voigt notation)
    /// via the 2th order deformation gradient (also in matrix notation)
    double PullBack4Tensorijkl(const LINALG::Matrix<3,3>& defgr,
                          const double (&FourTensor)[3][3][3][3],
                          const double& i,
                          const double& j,
                          const double& k,
                          const double& l);

    //Setup 4-Tensor from 6x6 Voigt notation
    void Setup4Tensor(
        double (&FourTensor)[3][3][3][3],
        const LINALG::Matrix<6,6>& VoigtMatrix
    );

//    //Setup 6x6 matrix in Voigt notation from 4-Tensor
//    void Setup6x6VoigtMatrix(
//        LINALG::Matrix<6,6>& VoigtMatrix,
//        const double(&FourTensor)[3][3][3][3]
//    );

    //Print Four Tensor
   void PrintFourTensor( double (&FourTensor)[3][3][3][3] );

  private:

    ///Normal \otime Normal
    LINALG::Matrix<3,3> NdN_; //would be fine if it was constant...

    ///T_1 \otimes T_1 + T_2 \otimes T_2
    LINALG::Matrix<3,3> TdT_; //would be fine if it was constant...

//    /// my material parameters
//    MAT::PAR::GrowthScd* paramsScdAC_;

  };//class GrowthScdACRadial

}  // namespace MAT


#endif /* GROWTH_SCD_H_ */
