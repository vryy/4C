/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores parameters of ion species in electrolyte solution. Former file of Georg
Bauer

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------------*/

#ifndef IONIC_H
#define IONIC_H

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class Ion : public Parameter
    {
     public:
      /// standard constructor
      Ion(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Ion() { ; }

      /// @name material parameters
      //@{

      /// valence (= charge number)
      const double valence_;
      /// diffusivity coefficient
      const double diffusivity_;
      /// densification coefficient
      const double densification_;
      /// valence (= charge number) of eliminated ion species
      const double elimvalence_;
      /// diffusivity coefficient of eliminated ion species
      const double elimdiffusivity_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Ion

  }  // namespace PAR

  class IonType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "IonType"; }

    static IonType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static IonType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class Ion : public Material
  {
   public:
    /// construct empty material object
    Ion();

    /// construct the material object given material parameters
    explicit Ion(MAT::PAR::Ion* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return IonType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_ion; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Ion(*this)); }

    /// valence (= charge number)
    double Valence() const { return params_->valence_; }
    /// diffusivity coefficient
    double Diffusivity() const { return params_->diffusivity_; }
    /// densification coefficient
    double Densification() const { return params_->densification_; }
    /// valence (= charge number) of eliminated ion species
    double ElimValence() const { return params_->elimvalence_; }
    /// diffusivity coefficient of eliminated ion species
    double ElimDiffusivity() const { return params_->elimdiffusivity_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::Ion* params_;
  };

}  // namespace MAT

#endif
