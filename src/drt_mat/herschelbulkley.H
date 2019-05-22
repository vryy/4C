/*----------------------------------------------------------------------*/
/*!
\brief non-Newtonian fluid of Herschel-Bulkley type

\level 2

\maintainer Volker Gravemeier
*/
/*----------------------------------------------------------------------*/
#ifndef HERSCHELBULKLEY_H
#define HERSCHELBULKLEY_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class HerschelBulkley : public Parameter
    {
     public:
      /// standard constructor
      HerschelBulkley(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~HerschelBulkley() { ; }

      /// @name material parameters
      //@{

      const double tau0_;           /* yield stress */
      const double kfac_;           /* constant factor */
      const double nexp_;           /* exponent */
      const double mexp_;           /* exponent */
      const double lolimshearrate_; /* lower limit of shear rate */
      const double uplimshearrate_; /* upper limit of shear rate */
      const double density_;        /* density */

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class HerschelBulkley

  }  // namespace PAR

  class HerschelBulkleyType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "HerschelBulkleyType"; }

    static HerschelBulkleyType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static HerschelBulkleyType instance_;
  };

  /// Nonlinear viscosity according to HerschelBulkley
  class HerschelBulkley : public Material
  {
   public:
    /// construct empty material object
    HerschelBulkley();

    /// construct the material object given material parameters
    explicit HerschelBulkley(MAT::PAR::HerschelBulkley* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const
    {
      return HerschelBulkleyType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_herschelbulkley; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new HerschelBulkley(*this));
    }

    /// return material parameters for element calculation
    //@{

    /// yield stress
    double Tau0() const { return params_->tau0_; }
    /// constant factor
    double KFac() const { return params_->kfac_; }
    /// exponent
    double NExp() const { return params_->nexp_; }
    /// exponent
    double MExp() const { return params_->mexp_; }
    /// lower limit of shear rate
    double LoLimShearRate() const { return params_->lolimshearrate_; }
    /// upper limit of shear rate
    double UpLimShearRate() const { return params_->uplimshearrate_; }
    /// density
    virtual double Density() const { return params_->density_; }

    //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::HerschelBulkley* params_;
  };

}  // namespace MAT

#endif
