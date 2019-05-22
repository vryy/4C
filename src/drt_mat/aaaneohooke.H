/*!----------------------------------------------------------------------
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke YOUNG 1.044E7 BETA 188.1E5 DENS 1.0

\level 3

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#ifndef AAANEOHOOKE_H_
#define AAANEOHOOKE_H_

#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <Epetra_Map.h>


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
      AAAneohooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~AAAneohooke() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      void OptParams(std::map<std::string, int>* pnames);

      // !brief enum for mapping between material parameter and entry in the matparams_ vector
      enum matparamnames_
      {
        young,
        nue,
        beta,
        density,
        first = young,
        last = density
      };

    };  // class AAAneohooke

  }  // namespace PAR

  class AAAneohookeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "AAAneohookeType"; }

    static AAAneohookeType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return AAAneohookeType::Instance().UniqueParObjectId();
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
    // virtual double Density() const { return params_->GetDensity(); }
    virtual double Density() const { return params_->GetParameter(params_->density, -1); }

    /// shear modulus
    // double ShearMod() const { return 0.5*params_->GetYoungs(-1)/(1.0+params_->GetNue()); }
    double ShearMod() const
    {
      return 0.5 * params_->GetParameter(params_->young, -1) /
             (1.0 + params_->GetParameter(params_->nue, -1));
    }

    // material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_aaaneohooke; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new AAAneohooke(*this)); }

    // THE material routine
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int eleGID);

    /// evaluate strain energy function
    virtual void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain, double& psi, const int eleGID);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

   private:
    /// my material parameters
    MAT::PAR::AAAneohooke* params_;
  };
}  // namespace MAT

#endif  // AAANEOHOOKE_H_
