/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_poro_ecm.H

 \brief scatra material for transport within porous model with special implementations
        for ECM model

   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de

\level 3
 *----------------------------------------------------------------------*/

#ifndef SCATRA_MAT_PORO_ECM_H_
#define SCATRA_MAT_PORO_ECM_H_

#include "scatra_reaction_mat.H"

namespace MAT
{

namespace PAR
{

/*----------------------------------------------------------------------*/
/// parameters for scalar transport material
class ScatraMatPoroECM: public ScatraReactionMat
{
public:

  /// standard constructor
  ScatraMatPoroECM(Teuchos::RCP<MAT::PAR::Material> matdata);

  /// destructor
  virtual ~ScatraMatPoroECM()
  {
    ;
  }

  /// create material instance of matching type with my parameters
  virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  double reacscale_;

};
// class Scatra

}// namespace PAR

class ScatraMatPoroECMType: public ScatraReactionMatType
{
public:

  std::string Name() const
  {
    return "ScatraMatPoroECMType";
  }

  static ScatraMatPoroECMType & Instance()
  {
    return instance_;
  }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

private:

  static ScatraMatPoroECMType instance_;
};

/*----------------------------------------------------------------------*/
/// wrapper for scalar transport material
class ScatraMatPoroECM: public ScatraReactionMat
{
public:

  /// construct empty material object
  ScatraMatPoroECM();

  /// construct the material object given material parameters
  explicit ScatraMatPoroECM(MAT::PAR::ScatraMatPoroECM* params);

  //! @name Packing and Unpacking

  /*!
   \brief Return unique ParObject id

   every class implementing ParObject needs a unique id defined at the
   top of drt_parobject.H (this file) and should return it in this method.
   */
  virtual int UniqueParObjectId() const
  {
    return ScatraReactionMatType::Instance().UniqueParObjectId();
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
    return INPAR::MAT::m_scatra_reaction_poroECM;
  }

  /// return copy of this material object
  virtual Teuchos::RCP<Material> Clone() const
  {
    return Teuchos::rcp(new ScatraMatPoroECM(*this));
  }

  /// Return quick accessible material parameter data
  virtual MAT::PAR::Parameter* Parameter() const
  {
    return params_;
  }

  /// return reaction coefficient
  virtual double ReacCoeff() const
  {
    return reaccoeff_;
  }

  /// compute reaction coefficient from structure chemical potential
  void ComputeReacCoeff(double chempot);

private:

  /// my material parameters
  MAT::PAR::ScatraMatPoroECM* params_;

  /// reaction coefficient
  double reaccoeff_;

};

}

#endif /* SCATRA_MAT_PORO_ECM_H_ */
