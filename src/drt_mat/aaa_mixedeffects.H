/*!----------------------------------------------------------------------
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke AGE 67 REFDIA 22.5 NUE 0.49 DENS 0.0001

\level 3

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#ifndef AAA_MIXEDEFFECTS_H_
#define AAA_MIXEDEFFECTS_H_

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
    class AAA_mixedeffects : public Parameter
    {
     public:
      /// standard constructor
      AAA_mixedeffects(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~AAA_mixedeffects() { ; }

      /// @name material parameters
      //@{

      /// Possion's ratio
      const double nue_;
      /// patient's age for calculating alpha and beta
      const double age_;
      /// reference diameter for normalizing (=to calculate NORD)
      const double refdia_;
      /// mass density
      const double density_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class AAA_mixedeffects

  }  // namespace PAR

  class AAA_mixedeffectsType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "AAA_mixedeffectsType"; }

    static AAA_mixedeffectsType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static AAA_mixedeffectsType instance_;
  };



  /*----------------------------------------------------------------------*/

  class AAA_mixedeffects : public So3Material
  {
   public:
    // empty constructor
    AAA_mixedeffects();

    // constructor with given material parameters
    AAA_mixedeffects(MAT::PAR::AAA_mixedeffects* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return AAA_mixedeffectsType::Instance().UniqueParObjectId();
    }

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
    double ShearMod(double elelocalrad) const
    {
      return  // => 0.5*6*alpha/(1.0+nue)
          3E6 * (0.09631 + 0.03329 * (elelocalrad * 2 / params_->refdia_ - 2.55)) /
          (1.0 + params_->nue_);
    }


    // material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_aaa_mixedeffects; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new AAA_mixedeffects(*this));
    }

    // THE material routine
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e, Epetra_SerialDenseMatrix* cmat_e,
        Epetra_SerialDenseVector* stress_e, double elelocalrad);

    // THE material routine
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int eleGID);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }


   private:
    /// my material parameters
    MAT::PAR::AAA_mixedeffects* params_;
  };
}  // namespace MAT

#endif  // AAA_MIXEDEFFECTS_H_
