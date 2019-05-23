/*!----------------------------------------------------------------------*/
/*!
\brief material for heat transport due to Fourier-type thermal conduction and the Soret effect

\maintainer Christoph Schmidt

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef SORET_H
#define SORET_H

#include "fourieriso.H"

namespace MAT
{
  namespace PAR
  {
    //! parameters for Soret material
    class Soret : public FourierIso
    {
     public:
      //! constructor
      Soret(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~Soret() { return; };

      //! create instance of Soret material
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      //! return Soret coefficient
      double SoretCoefficient() const { return soretcoefficient_; };

     private:
      //! @name parameters for Soret material
      //@{
      //! Soret coefficient
      const double soretcoefficient_;
      //@}
    };  // class MAT::PAR::Soret
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class SoretType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "SoretType"; };

    static SoretType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static SoretType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for Soret material
  class Soret : public FourierIso
  {
   public:
    //! construct empty Soret material
    Soret();

    //! construct Soret material with specific material parameters
    explicit Soret(MAT::PAR::Soret* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const { return SoretType::Instance().UniqueParObjectId(); };

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    void Unpack(const std::vector<char>& data);
    //@}

    //! return material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_soret; };

    //! clone Soret material
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Soret(*this)); };

    //! return Soret coefficient
    double SoretCoefficient() const { return params_->SoretCoefficient(); };

   private:
    //! return material parameters
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! material parameters
    MAT::PAR::Soret* params_;
  };
}  // namespace MAT
#endif
