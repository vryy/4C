/*!----------------------------------------------------------------------
\file humphreycardiovascular.H
\brief
This file contains routines for an anisotropic material with amorphous elastin,
four fiber families of collagen (circumferential, axial, diagonal) and one
fiber family of smooth muscle (circumferential).
The strain energy function for one fiber family is the same as in holzapfelcardiovascular.

example input line
MAT 1 MAT_HUMPHREYCARDIO KAPPA 0.833 MUE 0.385 DENS 1.0 K1C 1.0 K2C 1.0 K1M 1.0 K2M 1.0 PHIE 0.02 PHIC 0.22 PHIM 0.76 INIT 1

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifndef HUMPHREYCARDIOVASCULAR_H_
#define HUMPHREYCARDIOVASCULAR_H_


#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class HumphreyCardio : public Parameter
    {
    public:

      /// standard constructor
      HumphreyCardio(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~HumphreyCardio() { ; }

      /// @name material parameters
      //@{
      /// dilatation parameter
      const double kappa_;
      /// shear modulus
      const double mue_;
      /// density
      const double density_;
      /// parameter for linear fiber stiffness of collagen
      const double k1c_;
      /// parameter for exponential fiber stiffness of collagen
      const double k2c_;
      /// parameter for linear fiber stiffness of smooth muscle
      const double k1m_;
      /// parameter for exponential fiber stiffness of smooth muscle
      const double k2m_;
      /// mass fraction of elastin
      const double phie_;
      /// mass fraction of collagen
      const double phic_;
      /// mass fraction of smooth muscle
      const double phim_;
      /// option for initialization
      const int init_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class HumphreyCardio

  }  // namespace PAR

class HumphreyCardioType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "HumphreyCardioType"; }

  static HumphreyCardioType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static HumphreyCardioType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for anisotropic arterie material
  class HumphreyCardio : public Material
  {
  public:

    /// construct empty material object
    HumphreyCardio();

    /// construct the material object given material parameters
    explicit HumphreyCardio(MAT::PAR::HumphreyCardio* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return HumphreyCardioType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_humphreycardiovascular; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new HumphreyCardio( *this ) ); }

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
    Teuchos::RCP<std::vector<std::vector<double> > > Geta1() const {return a1_;};

    /// return a2s
    Teuchos::RCP<std::vector<std::vector<double> > > Geta2() const {return a2_;};

    /// return a3s
    Teuchos::RCP<std::vector<std::vector<double> > > Geta3() const {return a3_;};

    /// return a4s
    Teuchos::RCP<std::vector<std::vector<double> > > Geta4() const {return a4_;};

    /// evaluate fiber directions from locsys, pull back
    void EvaluateFiberVecs
    (const int gp, const LINALG::Matrix<3,3>& locsys, const LINALG::Matrix<3,3>& defgrd);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::HumphreyCardio* params_;

    // internal variables for fibers
    Teuchos::RCP<std::vector<std::vector<double> > > a1_; ///< first fiber vector per gp (reference), circumferential
    Teuchos::RCP<std::vector<std::vector<double> > > a2_; ///< second fiber vector per gp (reference), axial
    Teuchos::RCP<std::vector<std::vector<double> > > a3_; ///< third fiber vector per gp (reference), diagonal
    Teuchos::RCP<std::vector<std::vector<double> > > a4_; ///< fourth fiber vector per gp (reference), diagonal
    Teuchos::RCP<std::vector<std::vector<double> > > ca1_; ///< first fiber vector per gp (spatial config)
    Teuchos::RCP<std::vector<std::vector<double> > > ca2_; ///< second fiber vector per gp (spatial config)
    Teuchos::RCP<std::vector<std::vector<double> > > ca3_; ///< third fiber vector per gp (spatial config)
    Teuchos::RCP<std::vector<std::vector<double> > > ca4_; ///< fourth fiber vector per gp (spatial config)

    bool isinit_;  ///< indicates if material is initialized

  };
}  // namespace MAT

#endif /*HUMPHREYCARDIOVASCULAR_H_*/
