/*-----------------------------------------------------------*/
/*! \file
\brief A class for a crosslinker material


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_MAT_CROSSLINKERMAT_HPP
#define FOUR_C_MAT_CROSSLINKERMAT_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for de St. Venant--Kirchhoff
    class CrosslinkerMat : public CORE::MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      CrosslinkerMat(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

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
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;

    };  // class CrosslinkerMat

  }  // namespace PAR

  class CrosslinkerMatType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "CrosslinkerMatType"; }

    static CrosslinkerMatType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static CrosslinkerMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Crosslinker Material
  class CrosslinkerMat : public CORE::MAT::Material
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
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
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
    void Pack(CORE::COMM::PackBuffer& data) const override;

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
    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_crosslinkermat;
    }

    /// return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new CrosslinkerMat(*this));
    }

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
    CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

    //@}

   private:
    /// my material parameters
    MAT::PAR::CrosslinkerMat* params_;
  };
}  // namespace MAT



FOUR_C_NAMESPACE_CLOSE

#endif
