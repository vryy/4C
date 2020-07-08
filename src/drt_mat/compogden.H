/*----------------------------------------------------------------------*/
/*! \file
\brief Compressible Ogden hyperelastic material (with shell8)

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef COMPOGDEN_H_
#define COMPOGDEN_H_


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
    /// material parameters for compress. Ogden
    class CompOgden : public Parameter
    {
     public:
      /// standard constructor
      CompOgden(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~CompOgden() { ; }

      /// @name material parameters
      //@{
      int init_;        /* init flag */
      double nue_;      /* Possion ratio */
      double beta_;     /* the unphysical material constant called beta */
      double alfap_[3]; /* three parameters alfap */
      double mup_[3];   /* three parameters nuep */
      double density_;  /* material specific weight */
      double lambda_;   /* 1. lame constant */
      double kappa_;    /* bulkmodulus */
      double l_[3];
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class CompOgden

  }  // namespace PAR

  class CompOgdenType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "CompOgdenType"; }

    static CompOgdenType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static CompOgdenType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for compressible Ogden material
  class CompOgden : public Material
  {
   public:
    /// construct empty material object
    CompOgden();

    /// construct the material object given material parameters
    explicit CompOgden(MAT::PAR::CompOgden* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return CompOgdenType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_compogden; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new CompOgden(*this)); }

    /// get parameters
    int Init() const { return params_->init_; }
    double Nue() const { return params_->nue_; }
    double Beta() const { return params_->beta_; }
    double AlfaP(const unsigned i) const { return params_->alfap_[i]; }
    double MuP(const unsigned i) const { return params_->mup_[i]; }
    virtual double Density() const { return params_->density_; }
    double Lambda() const { return params_->lambda_; }
    double Kappa() const { return params_->kappa_; }
    double L(const unsigned i) const { return params_->l_[i]; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::CompOgden* params_;
  };

}  // namespace MAT

#endif /*COMPOGDEN_H_*/
