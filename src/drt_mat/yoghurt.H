/*----------------------------------------------------------------------*/
/*!
\file yoghurt.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef YOGHURT_H
#define YOGHURT_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for "Yoghurt-type" material with nonlinear viscosity
    /// determined via power law by Oswald and extended by an Arrhenius-type
    /// term to account for temperature dependence
    class Yoghurt : public Parameter
    {
     public:
      /// standard constructor
      Yoghurt(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Yoghurt() { ; }

      /// @name material parameters
      //@{

      /// specific heat capacity
      const double shc_;
      /// density
      const double density_;
      /// thermal conductivity
      const double thermcond_;
      /// exponent of strain-rate term
      const double strrateexp_;
      /// pre-exponential constant
      const double preexcon_;
      /// activation energy
      const double actenergy_;
      /// specific gas constant R
      const double gasconst_;
      /// safety factor delta
      const double delta_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Yoghurt

  }  // namespace PAR

  class YoghurtType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "YoghurtType"; }

    static YoghurtType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static YoghurtType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for "Yoghurt-type" material
  class Yoghurt : public Material
  {
   public:
    /// construct empty material object
    Yoghurt();

    /// construct the material object given material parameters
    explicit Yoghurt(MAT::PAR::Yoghurt* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return YoghurtType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_yoghurt; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Yoghurt(*this)); }

    /// compute viscosity
    double ComputeViscosity(const double rateofstrain, const double temp) const;

    /// compute diffusivity
    double ComputeDiffusivity() const;

    /// return material parameters for element calculation
    //@{

    /// specific heat capacity
    double Shc() const { return params_->shc_; }
    /// density
    double Density() const { return params_->density_; }
    /// thermal conductivity
    double ThermCond() const { return params_->thermcond_; }
    /// exponent of strain-rate term
    double StrRateExp() const { return params_->strrateexp_; }
    /// pre-exponential constant
    double PreExCon() const { return params_->preexcon_; }
    /// activation energy
    double ActEnergy() const { return params_->actenergy_; }
    /// specific gas constant R
    double GasConst() const { return params_->gasconst_; }
    /// safety factor
    double Delta() const { return params_->delta_; }  //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::Yoghurt* params_;
  };

}  // namespace MAT

#endif
