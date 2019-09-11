/*----------------------------------------------------------------------*/
/*! \file
\brief Material law for elastic spring (either translational or rotational spring)

\level 3

\maintainer Jan Schnabel
*----------------------------------------------------------------------*/
#ifndef SPRING_H
#define SPRING_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for spring
    class Spring : public Parameter
    {
     public:
      /// standard constructor
      Spring(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Spring() { ; }

      /// @name material parameters
      //@{

      /// stiffness (translational or rotational)
      const double stiffness_;

      /// density
      const double density_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Spring

  }  // namespace PAR

  class SpringType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "SpringType"; }

    static SpringType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static SpringType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for spring material
  class Spring : public Material
  {
   public:
    /// construct empty material object
    Spring();

    /// construct the material object given material parameters
    explicit Spring(MAT::PAR::Spring* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return SpringType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_spring; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Spring(*this)); }

    /// stiffness (translational or rotational)
    double Stiffness() const { return params_->stiffness_; }

    /// density
    double Density() const { return params_->density_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::Spring* params_;
  };

}  // namespace MAT

#endif
