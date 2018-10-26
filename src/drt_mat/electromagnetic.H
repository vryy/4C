/*!----------------------------------------------------------------------
\file electromagnetic.H
\brief Contains conductivity, permittivity and permeability of the medium for isotropic
       electromagetic field evolution.
       example input line:
       MAT 1 MAT_Electromagnetic CONDUCTIVITY 0.0 PERMITTIVITY 1.732 PERMEABILITY 1.732
<pre>

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244

</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef ELECTROMAGNETIC_H
#define ELECTROMAGNETIC_H

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
    /// material parameters for electromagnetic wave propagation
    class ElectromagneticMat : public Parameter
    {
     public:
      /// standard constructor
      ElectromagneticMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ElectromagneticMat() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      void OptParams(std::map<std::string, int>* pnames);

      enum matparamnames_
      {
        sigma_,
        epsilon_,
        mu_,
        first = sigma_,
        last = mu_
      };

    };  // class ElectromagneticMat

  }  // namespace PAR

  class ElectromagneticMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ElectromagneticMatType"; }

    static ElectromagneticMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ElectromagneticMatType instance_;

  };  // class ElectromagneticMatType

  /*----------------------------------------------------------------------*/
  /// Wrapper for Sound propagation material
  class ElectromagneticMat : public Material
  {
   public:
    /// construct empty material object
    ElectromagneticMat();

    /// construct the material object given material parameters
    explicit ElectromagneticMat(MAT::PAR::ElectromagneticMat* params);

    //! @name Packing and Unpacking

    /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ElectromagneticMatType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_electromagneticmat;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ElectromagneticMat(*this));
    }

    /// conductivity
    double sigma(int eleid = -1) const { return params_->GetParameter(params_->sigma_, eleid); }

    /// permittivity coefficient
    double epsilon(int eleid = -1) const { return params_->GetParameter(params_->epsilon_, eleid); }

    /// permeability coefficient
    double mu(int eleid = -1) const { return params_->GetParameter(params_->mu_, eleid); }



    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //@}


   private:
    /// my material parameters
    MAT::PAR::ElectromagneticMat* params_;
  };

}  // namespace MAT



#endif /* ELECTROMAGNETIC_H */
