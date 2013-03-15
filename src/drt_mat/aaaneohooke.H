/*!----------------------------------------------------------------------
\file aaaneohooke.H
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke YOUNG 1.044E7 BETA 188.1E5 DENS 1.0

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/
#ifndef AAANEOHOOKE_H_
#define AAANEOHOOKE_H_

#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for aneurysm wall material
    class AAAneohooke : public Parameter
    {
    public:

      /// standard constructor
      AAAneohooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~AAAneohooke() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double nue_;
      /// 2nd parameter
      const double beta_;
      /// mass density
      const double density_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class AAAneohooke

  }  // namespace PAR

class AAAneohookeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "AAAneohookeType"; }

  static AAAneohookeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static AAAneohookeType instance_;
};

  /*----------------------------------------------------------------------*/
  /// aneurysm wall material according to Raghavan and Vorp [2000]
  class AAAneohooke : public So3Material
  {

  public:

    // empty constructor
    AAAneohooke();

    // constructor with given material parameters
    AAAneohooke(MAT::PAR::AAAneohooke* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return AAAneohookeType::Instance().UniqueParObjectId(); }

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5*params_->youngs_/(1.0+params_->nue_); }

    // material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_aaaneohooke; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new AAAneohooke( *this ) ); }

    // THE material routine
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e,
		  Epetra_SerialDenseMatrix* cmat_e,
		  Epetra_SerialDenseVector* stress_e);

    // THE material routine
    virtual void Evaluate(
      const LINALG::Matrix<3,3>* defgrd,
      const LINALG::Matrix<6,1>* glstrain,
      Teuchos::ParameterList& params,
      LINALG::Matrix<6,1>* stress,
      LINALG::Matrix<6,6>* cmat);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::AAAneohooke* params_;

  };
}

#endif  // AAANEOHOOKE_H_
