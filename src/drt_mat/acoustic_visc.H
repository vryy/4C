/*!----------------------------------------------------------------------
\file acoustic_visc.H
\brief contains a density, a speed of sound and a viscsity to characterize
       all necessary properties for sound propagation in acoustic material
       example input line:
       MAT 1 MAT_AcousticVisc DENSITY 1000.0 C 1500.0 VISC  1.0
<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef ACOUSTIC_VISC_H
#define ACOUSTIC_VISC_H

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
class AcousticViscMat : public Parameter
{
public:

  /// standard constructor
  AcousticViscMat(
      Teuchos::RCP<MAT::PAR::Material> matdata
  );

  /// destructor
  virtual ~AcousticViscMat() { ; }

  /// @name material parameters
  //@{

  /// speed of sound
  const double c_;
  /// mass density
  const double density_;
  /// combined viscosities 4/3 mu + mu_B
  const double visc_;
  /// kappa * beta^2 * T_0 / c_p
  const double therm_;

  //@}

  /// create material instance of matching type with my parameters
  virtual Teuchos::RCP<MAT::Material> CreateMaterial();

}; // class AcousticMat

} // namespace PAR

class AcousticViscMatType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "AcousticViscMatType"; }

  static AcousticViscMatType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static AcousticViscMatType instance_;

}; // class AcousticViscMatType

/*----------------------------------------------------------------------*/
/// Wrapper for Sound propagation material
class AcousticViscMat : public Material
{
public:

  /// construct empty material object
  AcousticViscMat();

  /// construct the material object given material parameters
  explicit AcousticViscMat(MAT::PAR::AcousticViscMat* params);

  //! @name Packing and Unpacking

  /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
   */
  virtual int UniqueParObjectId() const { return AcousticViscMatType::Instance().UniqueParObjectId(); }

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
  virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_acousticviscmat; }

  /// return copy of this material object
  virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new AcousticViscMat( *this ) ); }



  /// speed of sound
  double SpeedofSound() const { return params_->c_; }

  /// density
  virtual double Density() const { return params_->density_; }

  /// combined viscosities 4/3 mu + mu_B
  double Viscosity() const { return params_->visc_; }

  /// combined thermal quantities
  double Therm() const { return params_->therm_; }

  /// Return quick accessible material parameter data
  virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  //@}


private:

  /// my material parameters
  MAT::PAR::AcousticViscMat* params_;

};

} // namespace MAT



#endif /* ACOUSTIC_H */
