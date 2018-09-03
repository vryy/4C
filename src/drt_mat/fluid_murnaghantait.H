/*----------------------------------------------------------------------*/
/*!
\file fluid_murnaghantait.H

\brief Weakly compressible fluid according to Murnaghan-Tait

<pre>
\level 1

\maintainer Andrea La Spina
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef MURNAGHANTAITFLUID_H
#define MURNAGHANTAITFLUID_H



#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for weakly compressible fluid according to Murnaghan-Tait
    ///
    /// This object exists only once for each read Murnaghan-Tait fluid.
    class MurnaghanTaitFluid : public Parameter
    {
     public:
      /// standard constructor
      MurnaghanTaitFluid(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~MurnaghanTaitFluid() { ; }

      /// @name material parameters
      //@{

      /// viscosity
      const double viscosity_;
      /// reference density
      const double refdensity_;
      /// reference pressure
      const double refpressure_;
      /// reference reference bulk modulus
      const double refbulkmodulus_;
      /// material parameter according to Murnaghan-Tait
      const double matparameter_;
      /// surface tension coefficient
      const double gamma_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class MurnaghanTaitFluid

  }  // namespace PAR

  class MurnaghanTaitFluidType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MurnaghanTaitFluidType"; }

    static MurnaghanTaitFluidType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static MurnaghanTaitFluidType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Murnaghan-Tait fluid material
  ///
  /// This object exists (several times) at every element
  class MurnaghanTaitFluid : public Material
  {
   public:
    /// construct empty material object
    MurnaghanTaitFluid();

    /// construct the material object given material parameters
    explicit MurnaghanTaitFluid(MAT::PAR::MurnaghanTaitFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return MurnaghanTaitFluidType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_fluid_murnaghantait;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new MurnaghanTaitFluid(*this));
    }

    /// compute density
    double ComputeDensity(const double press) const;

    /// return material parameters for element calculation
    //@{

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return reference density
    double RefDensity() const { return params_->refdensity_; }

    /// return reference pressure
    double RefPressure() const { return params_->refpressure_; }

    /// return reference bulk modulus
    double RefBulkModulus() const { return params_->refbulkmodulus_; }

    /// return material parameter according to Murnaghan-Tait
    double MatParameter() const { return params_->matparameter_; }

    /// return surface tension coefficient
    double Gamma() const { return params_->gamma_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::MurnaghanTaitFluid* params_;
  };

}  // namespace MAT

#endif
