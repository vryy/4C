/*!----------------------------------------------------------------------
\file aaaneohooke.H
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

This material law accounts for spatial variation of the material parameters young and beta
by using an random field, such that every element gets its own set of mat parameters.

The sample of the random field is not calculated in the material routine but globally in MLMC because
it cannot be stored in each element. At the moment beta is the only stochastic parameter.

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke_Stopro YOUNG 1.044E7 BETA 188.1E5 NUE 0.3 DENS 1.0 SIGMA 0.25 CORRLENGTH 5.0

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#ifndef AAANEOHOOKE_STOPRO_H_
#define AAANEOHOOKE_STOPRO_H_
#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

namespace MAT
{

  namespace PAR
  {


    /*----------------------------------------------------------------------*/
    /// material parameters for aneurysm wall material
    class AAAneohooke_stopro : public Parameter
    {
    public:

      /// standard constructor
      AAAneohooke_stopro(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~AAAneohooke_stopro() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double nue_;
      /// 2nd parameter
      const double beta_mean_;
      /// mass density
      const double density_;

      // Both stochastic parameters are not actually needed in material routine
      const double sigma_0_;
      const double corrlength_;
      //
      int init_;





      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class AAAneohooke

  }  // namespace PAR

class AAAneohooke_stoproType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "AAAneohooke_stoproType"; }

  static AAAneohooke_stoproType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static AAAneohooke_stoproType instance_;
};

  /*----------------------------------------------------------------------*/
  /// aneurysm wall material according to Raghavan and Vorp [2000]
  class AAAneohooke_stopro : public Material
  {

  public:

    // empty constructor
    AAAneohooke_stopro();

    // constructor with given material parameters
    AAAneohooke_stopro(MAT::PAR::AAAneohooke_stopro* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return AAAneohooke_stoproType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

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

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const vector<char>& data);

    //@}

    // stochastic beta
    double Beta() const {return beta_;}

    // variance of random field
    double Sigma() const {return params_->sigma_0_;}

    // coreelation length of random field
    double Corrlength() const {return params_->corrlength_ ;}
    /// material mass density
    double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5*params_->youngs_/(1.0+params_->nue_); }

    // material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_aaaneohooke_stopro; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new AAAneohooke_stopro( *this ) ); }

    // Init stochastic beta
    void Init(double value_stopro);

    // THE material routine
    //void Evaluate(const Epetra_SerialDenseVector* glstrain_e,
      //Epetra_SerialDenseMatrix* cmat_e,
      //Epetra_SerialDenseVector* stress_e);

    // THE material routine
    void Evaluate(const LINALG::Matrix<6,1>& glstrain,
      LINALG::Matrix<6,6>& cmat,
      LINALG::Matrix<6,1>& stress);


    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::AAAneohooke_stopro* params_;

    // stochastic beta
    // put beta not in params because it will not be initialized with constructor
    double beta_;
    /// init flag for stochastic beta
    bool isinit_;
  };
}

#endif  // CCADISCRET
#endif  // AAANEOHOOKE_STOPRO_H_
