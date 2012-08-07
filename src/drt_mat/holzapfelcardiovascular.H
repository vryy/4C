/*!----------------------------------------------------------------------
\file holzapfelcardiovascular.H
\brief
This file contains routines for an anisotropic material with two fiber families.
example input line
MAT 1 MAT_HOLZAPFELCARDIO KAPPA 0.833 MUE 0.385 DENS 1.0 K1 1.0 K2 1.0 GAMMA 45.0 MINSTRETCH 1.0 INIT 1

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifndef HOLZAPFELCARDIOVASCULAR_H_
#define HOLZAPFELCARDIOVASCULAR_H_


#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class HolzapfelCardio : public Parameter
    {
    public:

      /// standard constructor
      HolzapfelCardio(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~HolzapfelCardio() { ; }

      /// @name material parameters
      //@{
      /// dilatation parameter
      const double kappa_;
      /// shear modulus
      const double mue_;
      /// density
      const double density_;
      /// parameter for linear fiber stiffness
      const double k1_;
      /// parameter for exponential fiber stiffness
      const double k2_;
      /// angle between fiber direction and circumferential direction
      const double gamma_;
      /// minimal principal stretch fibers do respond to
      const double minstretch_;
      /// option for initialization
      const int init_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class HolzapfelCardio

  }  // namespace PAR

class HolzapfelCardioType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "HolzapfelCardioType"; }

  static HolzapfelCardioType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static HolzapfelCardioType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for anisotropic arterie material
  class HolzapfelCardio : public Material
  {
  public:

    /// construct empty material object
    HolzapfelCardio();

    /// construct the material object given material parameters
    explicit HolzapfelCardio(MAT::PAR::HolzapfelCardio* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return HolzapfelCardioType::Instance().UniqueParObjectId(); }

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_holzapfelcardiovascular; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new HolzapfelCardio( *this ) ); }

    /// Setup
    void Setup
    (
      const int numgp,  ///< number of Gauss points
      DRT::INPUT::LineDefinition* linedef
    );

    /// for visualisation of current fiber direction, can be called in so_material
    void UpdateFiberDirs
    (
      const int numgp,
      LINALG::Matrix<3,3>* defgrd
    );

    /// Evaluate material
    void Evaluate
    (
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,  ///<green lagrange strain
      const int gp,   ///< current Gauss point
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress  ///< 2nd PK-stress
    );

    /// Return density
    virtual double Density() const {return params_->density_;};

    /// Return shear modulus
    double ShearMod() const {return params_->mue_;};

    /// return a1s
    Teuchos::RCP<vector<vector<double> > > Geta1() const {return a1_;};

    /// return a2s
    Teuchos::RCP<vector<vector<double> > > Geta2() const {return a2_;};

    /// evaluate fiber directions from locsys and alignment angle, pull back
    void EvaluateFiberVecs
    (const int gp, const double gamma, const LINALG::Matrix<3,3>& locsys, const LINALG::Matrix<3,3>& defgrd);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::HolzapfelCardio* params_;

    // internal variables for fibers
    Teuchos::RCP<vector<vector<double> > > a1_; ///< first fiber vector per gp (reference)
    Teuchos::RCP<vector<vector<double> > > a2_; ///< second fiber vector per gp (reference)
    Teuchos::RCP<vector<vector<double> > > ca1_; ///< first fiber vector per gp (spatial config)
    Teuchos::RCP<vector<vector<double> > > ca2_; ///< second fiber vector per gp (spatial config)

    bool isinit_;  ///< indicates if material is initialized

  };
}  // namespace MAT

#endif /*HOLZAPFELCARDIOVASCULAR_H_*/
