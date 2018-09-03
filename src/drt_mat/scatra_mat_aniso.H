/*!----------------------------------------------------------------------
\file scatra_mat_aniso.H

\brief scatra_mat_aniso.H

\level 3

<pre>
  \maintainer Lasse Jagschies
              jagschies@mhpc.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-10365
</pre>
*----------------------------------------------------------------------*/
#ifndef SCATRA_MAT_ANISO_H
#define SCATRA_MAT_ANISO_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for anisotropic scalar transport material
    class ScatraMatAniso : public Parameter
    {
     public:
      /// standard constructor
      ScatraMatAniso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraMatAniso() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      void OptParams(std::map<std::string, int>* pnames);

      enum matparamnames_
      {
        diff1,
        diff2,
        diff3,
        first = diff1,
        last = diff3
      };

    };  // class Scatra

  }  // namespace PAR

  class ScatraMatAnisoType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ScatraMatAnisoType"; }

    static ScatraMatAnisoType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatAnisoType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for anisotropic scalar transport material
  class ScatraMatAniso : public Material
  {
   public:
    /// construct empty material object
    ScatraMatAniso();

    /// construct the material object given material parameters
    explicit ScatraMatAniso(MAT::PAR::ScatraMatAniso* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return ScatraMatAnisoType::Instance().UniqueParObjectId();
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

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra_aniso; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ScatraMatAniso(*this)); }

    /// diffusivity
    LINALG::Matrix<3, 1> Diffusivity(int eleid = -1) const
    {
      LINALG::Matrix<3, 1> diff;
      diff(0, 0) = params_->GetParameter(params_->diff1, eleid);
      diff(1, 0) = params_->GetParameter(params_->diff2, eleid);
      diff(2, 0) = params_->GetParameter(params_->diff3, eleid);

      return diff;
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMatAniso* params_;
  };

}  // namespace MAT

#endif
