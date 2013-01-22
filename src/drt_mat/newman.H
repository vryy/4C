/*!----------------------------------------------------------------------*/
/*!
\file newman.H

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef NEWMAN_H
#define NEWMAN_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class Newman : public Parameter
    {
    public:

      /// standard constructor
      Newman(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Newman() { ; }

      /// @name material parameters
      //@{

      /// valence (= charge number)
      const double valence_;
      /// diffusivity coefficient
      const int curvediff_;
      /// transference number
      const int curvetrans_;

      /// dilute solution theory (diffusion potential in current equation):
      ///    a          b
      ///   |--|  |----------|
      ///   z_1 + (z_2 - z_1) t_1
      ///  ---------------------
      ///      z_1 z_2
      ///     |________|
      ///         c

      /// constant for diffusion potential in current equation
      const double a_;
      /// constant for diffusion potential in current equation
      const double b_;
      /// constant for diffusion potential in current equation
      const double c_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Ion

  } // namespace PAR

class NewmanType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "NewmanType"; }

  static NewmanType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static NewmanType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class Newman : public Material
  {
  public:

    /// construct empty material object
    Newman();

    /// construct the material object given material parameters
    explicit Newman(MAT::PAR::Newman* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return NewmanType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_newman; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Newman( *this ) ); }

    /// valence (= charge number)
    double Valence() const { return params_->valence_; }
    /// curve representing the diffusion coefficient
    int CurveDiff() const { return params_->curvediff_; }
    /// curve representing the transference number
    int CurveTrans() const { return params_->curvetrans_; }

    /// return constant for diffusion potential in current equation
    double A() const { return params_->a_; }
    /// return constant for diffusion potential in current equation
    double B() const { return params_->b_; }
    /// return constant for diffusion potential in current equation
    double C() const { return params_->c_; }


    /// computation of the diffusion coefficient based on the defined curve
    double ComputeDiffusionCoefficient (const double cint) const;
    /// computation of the first derivative of the diffusion coefficient based on the defined curve
    double ComputeFirstDerivDiffCoeff (const double cint) const;

    /// computation of the transference number based on the defined curve
    double ComputeTransferenceNumber(const double cint) const;
    /// computation of the first derivative of the transference number based on the defined curve
    double ComputeFirstDerivTrans (const double cint) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::Newman* params_;
  };

}

#endif

