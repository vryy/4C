/*!----------------------------------------------------------------------
\file acoustic.H
\brief contains a density and a speed of sound to characterize all
       necessary properties for sound propagation in acoustic material
       example input line:
       MAT 1 MAT_Acoustic DENSITY 1000.0 C 1500.0
<pre>
\level 2
\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef ACOUSTIC_H
#define ACOUSTIC_H

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for sound propagation
    class AcousticMat : public Parameter
    {
     public:
      /// standard constructor
      AcousticMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~AcousticMat() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      void OptParams(std::map<std::string, int>* pnames);

      enum matparamnames_
      {
        density,
        c,
        first = density,
        last = c
      };

    };  // class AcousticMat

  }  // namespace PAR

  class AcousticMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "AcousticMatType"; }

    static AcousticMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static AcousticMatType instance_;

  };  // class AcousticMatType

  /*----------------------------------------------------------------------*/
  /// Wrapper for Sound propagation material
  class AcousticMat : public Material
  {
   public:
    /// construct empty material object
    AcousticMat();

    /// construct the material object given material parameters
    explicit AcousticMat(MAT::PAR::AcousticMat* params);

    //! @name Packing and Unpacking

    /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return AcousticMatType::Instance().UniqueParObjectId();
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

    //! @name Access methods

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_acousticmat; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new AcousticMat(*this)); }

    /// density
    double Density(int eleid = -1) const { return params_->GetParameter(params_->density, eleid); }

    /// reaction coefficient
    double SpeedofSound(int eleid = -1) const { return params_->GetParameter(params_->c, eleid); }



    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //@}


   private:
    /// my material parameters
    MAT::PAR::AcousticMat* params_;
  };

}  // namespace MAT



#endif /* ACOUSTIC_H */
