/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_multiporo.H

 \brief scatra material for transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_
#define SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_


#include "scatra_mat.H"

namespace MAT
{

namespace PAR
{

/*----------------------------------------------------------------------*/
/// parameters for scalar transport material
class ScatraMatMultiPoro: public ScatraMat
{
public:

  /// standard constructor
  ScatraMatMultiPoro(Teuchos::RCP<MAT::PAR::Material> matdata);

  /// destructor
  virtual ~ScatraMatMultiPoro()
  {
    ;
  }

  /// create material instance of matching type with my parameters
  virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  /// density
  const double density_;
  /// ID of fluid phase the scalar is associated with
  const int phaseID_;
  /// delta used for modelling dependency of diffusivity on (saturation*porosity)^delta
  //  as in G. Sciume, William G. Gray, F. Hussain, M. Ferrari, P. Decuzzi, and B. A. Schrefler.
  //  Three phase flow dynamics in tumor growth. Computational Mechanics, 53:465–484,
  //  2014.
  const double delta_;
  // minimum saturation under which also corresponding mass fraction is equal to zero
  const double min_sat_;

};
// class Scatra

}// namespace PAR

class ScatraMatMultiPoroType: public ScatraMatType
{
public:

  std::string Name() const
  {
    return "ScatraMatMultiPoroType";
  }

  static ScatraMatMultiPoroType & Instance()
  {
    return instance_;
  }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

private:

  static ScatraMatMultiPoroType instance_;
};

/*----------------------------------------------------------------------*/
/// wrapper for scalar transport material
class ScatraMatMultiPoro: public ScatraMat
{
public:

  /// construct empty material object
  ScatraMatMultiPoro();

  /// construct the material object given material parameters
  explicit ScatraMatMultiPoro(MAT::PAR::ScatraMatMultiPoro* params);

  //! @name Packing and Unpacking

  /*!
   \brief Return unique ParObject id

   every class implementing ParObject needs a unique id defined at the
   top of drt_parobject.H (this file) and should return it in this method.
   */
  virtual int UniqueParObjectId() const
  {
    return ScatraMatMultiPoroType::Instance().UniqueParObjectId();
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

  /// material type
  virtual INPAR::MAT::MaterialType MaterialType() const
  {
    return INPAR::MAT::m_scatra_multiporo;
  }

  /// return copy of this material object
  virtual Teuchos::RCP<Material> Clone() const
  {
    return Teuchos::rcp(new ScatraMatMultiPoro(*this));
  }

  /// Return quick accessible material parameter data
  virtual MAT::PAR::Parameter* Parameter() const
  {
    return params_;
  }

  /// return phase density
  virtual double Density() const
  {
    return params_->density_;
  }

  /// return phase ID
  virtual int PhaseID() const
  {
    return params_->phaseID_;
  }

  /// return delta
  virtual double Delta() const
  {
    return params_->delta_;
  }

  /// return minimum saturation
  virtual double MinSat() const
  {
    return params_->min_sat_;
  }

private:

  /// my material parameters
  MAT::PAR::ScatraMatMultiPoro* params_;

};

}



#endif /* SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_ */
