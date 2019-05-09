/*!----------------------------------------------------------------------
\brief scalar transport material

\level 1

\maintainer Volker Gravemeier
*----------------------------------------------------------------------*/
#ifndef SCATRA_MAT_H
#define SCATRA_MAT_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraMat : public Parameter
    {
     public:
      /// standard constructor
      ScatraMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraMat() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      void OptParams(std::map<std::string, int>* pnames);

      enum matparamnames_
      {
        diff,
        reac,
        densific,
        first = diff,
        last = densific
      };

    };  // class Scatra

  }  // namespace PAR

  class ScatraMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ScatraMatType"; }

    static ScatraMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraMat : public Material
  {
   public:
    /// construct empty material object
    ScatraMat();

    /// construct the material object given material parameters
    explicit ScatraMat(MAT::PAR::ScatraMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ScatraMatType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ScatraMat(*this)); }

    /// diffusivity
    double Diffusivity(int eleid = -1) const { return params_->GetParameter(params_->diff, eleid); }

    /// reaction coefficient
    double ReaCoeff(int eleid = -1) const { return params_->GetParameter(params_->reac, eleid); }

    /// densification coefficient
    double Densification(int eleid = -1) const
    {
      return params_->GetParameter(params_->densific, eleid);
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMat* params_;
  };

}  // namespace MAT

#endif
