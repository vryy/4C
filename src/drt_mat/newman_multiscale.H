/*----------------------------------------------------------------------*/
/*!
\file newman_multiscale.H

\brief material for macro-scale elements in multi-scale simulations of electrochemistry problems

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef NEWMAN_MULTISCALE_H
#define NEWMAN_MULTISCALE_H

#include "newman.H"
#include "scatra_multiscale.H"

namespace MAT
{
  namespace PAR
  {
    //! material parameters
    class NewmanMultiScale : public Newman, public ScatraMultiScale
    {
     public:
      //! constructor
      NewmanMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~NewmanMultiScale() { return; };

      //! create instance of Newman multi-scale material
      Teuchos::RCP<MAT::Material> CreateMaterial();

      //! return electronic conductivity
      double Sigma() const { return sigma_; };

     private:
      //! @name parameters for Newman multi-scale material
      //@{
      //! electronic conductivity
      const double sigma_;
      //@}
    };  // class MAT::PAR::NewmanMultiScale
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class NewmanMultiScaleType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "NewmanMultiScaleType"; };

    static NewmanMultiScaleType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static NewmanMultiScaleType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for Newman multi-scale material
  class NewmanMultiScale : public Newman, public ScatraMultiScale
  {
   public:
    //! construct empty Newman multi-scale material
    NewmanMultiScale();

    //! construct Newman multi-scale material with specific material parameters
    explicit NewmanMultiScale(MAT::PAR::NewmanMultiScale* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const { return NewmanMultiScaleType::Instance().UniqueParObjectId(); };

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
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_newman_multiscale; };

    //! clone Newman multi-scale material
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new NewmanMultiScale(*this)); };

    //! return electronic conductivity
    double Sigma() const { return params_->Sigma(); };

   private:
    //! return material parameters
    const MAT::PAR::ScatraMultiScale* Params() const { return params_; };

    //! material parameters
    MAT::PAR::NewmanMultiScale* params_;
  };  // wrapper for Newman multi-scale material
}  // namespace MAT
#endif
