/*----------------------------------------------------------------------*/
/*!
\brief material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_MAT_MULTISCALE_H
#define SCATRA_MAT_MULTISCALE_H

#include "scatra_mat.H"
#include "scatra_multiscale.H"

namespace MAT
{
  // forward declaration
  class ScatraMultiScaleGP;

  namespace PAR
  {
    //! material parameters
    class ScatraMatMultiScale : public ScatraMat, public ScatraMultiScale
    {
     public:
      //! constructor
      ScatraMatMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~ScatraMatMultiScale() { return; };

      //! create material
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      //! return porosity
      double Porosity() const { return porosity_; };

      //! return tortuosity
      double Tortuosity() const { return tortuosity_; };

     private:
      //! @name material parameters
      //@{
      //! porosity
      const double porosity_;

      //! tortuosity
      const double tortuosity_;
      //@}
    };  // class MAT::PAR::ScatraMatMultiScale
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class ScatraMatMultiScaleType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ScatraMatMultiScaleType"; };

    static ScatraMatMultiScaleType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatMultiScaleType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! material wrapper
  class ScatraMatMultiScale : public ScatraMat, public ScatraMultiScale
  {
   public:
    //! construct empty material
    ScatraMatMultiScale();

    //! construct material with specific material parameters
    explicit ScatraMatMultiScale(MAT::PAR::ScatraMatMultiScale* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const
    {
      return ScatraMatMultiScaleType::Instance().UniqueParObjectId();
    };

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
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra_multiscale; };

    //! clone material
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ScatraMatMultiScale(*this)); };

    //! return porosity
    double Porosity() const { return params_->Porosity(); };

    //! return tortuosity
    double Tortuosity() const { return params_->Tortuosity(); };

   private:
    //! return material parameters
    const MAT::PAR::ScatraMultiScale* Params() const { return params_; };

    //! material parameters
    MAT::PAR::ScatraMatMultiScale* params_;
  };
}  // namespace MAT
#endif
