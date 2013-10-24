/*!----------------------------------------------------------------------*/
/*!
\file elchphase.H

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ELCHPHASE_H
#define ELCHPHASE_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class ElchPhase : public Parameter
    {
    public:

      /// standard constructor
      ElchPhase(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ElchPhase() { ; }

      /// @name material parameters
      //@{

      /// valence (= charge number)
      const double epsilon_;
      const double tortuosity_;
      /// diffusivity coefficient
      const double conductivity_;

      const double condcurvenr_;

      const int curvetherm_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ElchPhase

  } // namespace PAR

class ElchPhaseType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ElchPhaseType"; }

  static ElchPhaseType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ElchPhaseType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class ElchPhase : public Material
  {
  public:

    /// construct empty material object
    ElchPhase();

    /// construct the material object given material parameters
    explicit ElchPhase(MAT::PAR::ElchPhase* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ElchPhaseType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elchphase; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ElchPhase( *this ) ); }

    /// valence (= charge number)
    double Epsilon() const { return params_->epsilon_; }
    double Tortuosity() const { return params_->tortuosity_; }
    double Conductivity() const { return params_->conductivity_; }
    /// valence (= charge number)
    double CondCurveNr() const { return params_->condcurvenr_; }

    /// curve representing the transference number
    int CurveTherm() const { return params_->curvetherm_; }

    double ComputeConductivity (const double cint) const ;
    double ComputeFirstDerivCond(const double cint) const;

    double ComputeThermodynamicFactor (const double cint) const ;
    double ComputeFirstDerivThermodynamicFactor(const double cint) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ElchPhase* params_;
  };

}

#endif

