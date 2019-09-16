/*--------------------------------------------------------------------------*/
/*! \file
\brief Material model for the lubrication film

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_MAT_LUBRICATION_MAT_H_
#define SRC_DRT_MAT_LUBRICATION_MAT_H_


#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    class LubricationLaw;

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class LubricationMat : public Parameter
    {
     public:
      /// standard constructor
      LubricationMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~LubricationMat() { ; }

      /// density
      const double density_;

      /// lubrication law ID
      int lubricationlawID_;

      //@}

      // implementation of lubrication law
      LubricationLaw* lubricationlaw_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Lubrication

  }  // namespace PAR

  class LubricationMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "LubricationMatType"; }

    static LubricationMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static LubricationMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class LubricationMat : public Material
  {
   public:
    /// construct empty material object
    LubricationMat();

    /// construct the material object given material parameters
    explicit LubricationMat(MAT::PAR::LubricationMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return LubricationMatType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_lubrication; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new LubricationMat(*this)); }

    /// compute current viscosity
    double ComputeViscosity(const double press  ///< (i) lubrication pressure
    );

    //! evaluate constitutive relation for viscosity and compute derivatives
    double ComputeViscosityDeriv(const double press, const double visc);

    /// density
    double Density() const { return params_->density_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::LubricationMat* params_;
  };

}  // namespace MAT


#endif /* SRC_DRT_MAT_LUBRICATION_MAT_H_ */
