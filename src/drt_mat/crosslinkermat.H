/*-----------------------------------------------------------*/
/*!
\file crosslinkermat.H

\brief A class for a crosslinker material

\maintainer Jonas Eichinger

\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef CROSSLINKER_H_
#define CROSSLINKER_H_

#include "matpar_parameter.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Teuchos_ParameterList.hpp>

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for de St. Venant--Kirchhoff
    class CrosslinkerMat : public Parameter
    {
     public:
      /// standard constructor
      CrosslinkerMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~CrosslinkerMat() { ; }

      /// number of material for underlying linker element
      int const link_element_matnum_;

      /// type of joint
      INPAR::BEAMINTERACTION::JointType jointtype_;

      /// distance between the two binding domains of a linker
      double const linkinglength_;

      /// tolerance for linker length in the sense: length +- tolerance
      double const linkinglengthtol_;

      /// preferred binding angle enclosed by two filaments' axes in radians
      double const linkingangle_;

      /// tolerance for preferred binding angle in radians in the sense of: angle +- tolerance
      double const linkingangletol_;

      /// chemical association-rate
      double const k_on_;

      /// chemical dissociation-rate
      double const k_off_;

      /// deltaD in Bell's equation for force dependent off rate
      double const deltabelleq_;

      /// distance to sphere elements in which binding is prohibited
      double const nobonddistsphere;

      /// type of crosslinker
      INPAR::BEAMINTERACTION::CrosslinkerType linkertype_;


      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class CrosslinkerMat

  }  // namespace PAR

  class CrosslinkerMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "CrosslinkerMatType"; }

    static CrosslinkerMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static CrosslinkerMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Crosslinker Material
  class CrosslinkerMat : public Material
  {
   public:
    /// construct empty material object
    CrosslinkerMat();

    /// construct the material object given material parameters
    explicit CrosslinkerMat(MAT::PAR::CrosslinkerMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return CrosslinkerMatType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_crosslinkermat; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new CrosslinkerMat(*this)); }

    /// number of beam elasthyper material
    virtual double BeamElastHyperMatNum() const { return params_->link_element_matnum_; }

    /// force depedent off rate according to bells equation
    virtual INPAR::BEAMINTERACTION::JointType JointType() const { return params_->jointtype_; };

    /// distance between the two binding domains of a linker
    virtual double LinkingLength() const { return params_->linkinglength_; }

    /// tolerance for linker length in the sense: length +- tolerance
    virtual double LinkingLengthTolerance() const { return params_->linkinglengthtol_; }

    /// preferred binding angle enclosed by two filaments' axes in radians
    virtual double LinkingAngle() const { return params_->linkingangle_; }

    /// tolerance for preferred binding angle in radians in the sense of: angle +- tolerance
    virtual double LinkingAngleTolerance() const { return params_->linkingangletol_; }

    /// chemical association rate
    virtual double const& KOn() const { return params_->k_on_; }

    /// chemical dissociation rate
    virtual double KOff() const { return params_->k_off_; }

    /// force depedent off rate according to bells equation
    virtual double DeltaBellEq() const { return params_->deltabelleq_; };

    /// distance to sphere elements in which no double bonded beam to beam linker is allowed
    virtual double NoBondDistSphere() const { return params_->nobonddistsphere; };

    /// force depedent off rate according to bells equation
    virtual INPAR::BEAMINTERACTION::CrosslinkerType LinkerType() const
    {
      return params_->linkertype_;
    };

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //@}

   private:
    /// my material parameters
    MAT::PAR::CrosslinkerMat* params_;
  };
}  // namespace MAT



#endif /* CROSSLINKER_H_ */
