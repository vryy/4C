/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke AGE 67 REFDIA 22.5 NUE 0.49 DENS 0.0001

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_AAA_MIXEDEFFECTS_HPP
#define FOUR_C_MAT_AAA_MIXEDEFFECTS_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_par_parameter.hpp"
#include "4C_mat_so3_material.hpp"

FOUR_C_NAMESPACE_OPEN


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for aneurysm wall material
    class AaaMixedeffects : public Parameter
    {
     public:
      /// standard constructor
      AaaMixedeffects(Teuchos::RCP<MAT::PAR::Material> matdata);

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
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class AAA_mixedeffects

  }  // namespace PAR

  class AaaMixedeffectsType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "AAA_mixedeffectsType"; }

    static AaaMixedeffectsType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static AaaMixedeffectsType instance_;
  };



  /*----------------------------------------------------------------------*/

  class AaaMixedeffects : public So3Material
  {
   public:
    // empty constructor
    AaaMixedeffects();

    // constructor with given material parameters
    AaaMixedeffects(MAT::PAR::AaaMixedeffects* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return AaaMixedeffectsType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

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
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material mass density
    double Density() const override { return params_->density_; }

    /// shear modulus
    double ShearMod(double elelocalrad) const
    {
      return  // => 0.5*6*alpha/(1.0+nue)
          3E6 * (0.09631 + 0.03329 * (elelocalrad * 2 / params_->refdia_ - 2.55)) /
          (1.0 + params_->nue_);
    }


    // material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_aaa_mixedeffects;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new AaaMixedeffects(*this));
    }

    // THE material routine
    void Evaluate(const CORE::LINALG::SerialDenseVector* glstrain_e,
        CORE::LINALG::SerialDenseMatrix* cmat_e, CORE::LINALG::SerialDenseVector* stress_e,
        double elelocalrad);

    // THE material routine
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
        const int eleGID) override;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }


   private:
    /// my material parameters
    MAT::PAR::AaaMixedeffects* params_;
  };
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
