/*!----------------------------------------------------------------------
\file diffcond.H

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*----------------------------------------------------------------------*/
#ifndef DIFFCOND_H
#define DIFFCOND_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class DiffCond : public Parameter
    {
    public:

      /// standard constructor
      DiffCond(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~DiffCond() { ; }

      /// @name material parameters
      //@{

      /// valence (= charge number)
      const double valence_;
      /// diffusivity coefficient
      const double diffusivity_;
      /// transference number
      const double transference_;


      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Ion

  } // namespace PAR

class DiffCondType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "DiffCondType"; }

  static DiffCondType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static DiffCondType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class DiffCond : public Material
  {
  public:

    /// construct empty material object
    DiffCond();

    /// construct the material object given material parameters
    explicit DiffCond(MAT::PAR::DiffCond* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return DiffCondType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_diffcond; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new DiffCond( *this ) ); }

    /// valence (= charge number)
    double Valence() const { return params_->valence_; }
    /// diffusivity coefficient
    double Diffusivity() const { return params_->diffusivity_; }
    /// densification coefficient
    double Transference() const { return params_->transference_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::DiffCond* params_;
  };

}

#endif

