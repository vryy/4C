/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a linear Reissner beam element used as mechanical pin joint
       between two other beam elements

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_PINJOINTED_HPP
#define FOUR_C_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_PINJOINTED_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_link_pinjointed.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace Discret
{
  namespace ELEMENTS
  {
    class Beam3r;
  }
}  // namespace Discret

namespace BEAMINTERACTION
{
  class BeamLinkBeam3rLine2PinJointedType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "BeamLinkBeam3rLine2PinJointedType"; };

    static BeamLinkBeam3rLine2PinJointedType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static BeamLinkBeam3rLine2PinJointedType instance_;
  };


  /*!
   \brief element for link between two 3D beam elements via a linear (2 noded) Simo-Reissner beam
   element
   */
  class BeamLinkBeam3rLine2PinJointed : public BeamLinkPinJointed
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamLinkBeam3rLine2PinJointed();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLinkBeam3rLine2PinJointed(const BeamLinkBeam3rLine2PinJointed& old);



    //! Initialization [derived]
    void Init(int id, const std::vector<std::pair<int, int>>& eleids,
        const std::vector<Core::LinAlg::Matrix<3, 1>>& initpos,
        const std::vector<Core::LinAlg::Matrix<3, 3>>& inittriad,
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset) override;

    //! Setup [derived]
    void Setup(const int matnum) override;

    /*!
    \brief Return unique ParObject id [derived]

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int UniqueParObjectId() const override
    {
      return BeamLinkBeam3rLine2PinJointedType::Instance().UniqueParObjectId();
    };

    /*!
    \brief Pack this class so it can be communicated [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Pack(Core::Communication::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Unpack(const std::vector<char>& data) override;

    /// return copy of this linking object
    Teuchos::RCP<BeamLink> Clone() const override;

    //@}

    //! @name Access methods

    //! get internal linker energy
    double GetInternalEnergy() const override;

    //! get kinetic linker energy
    double GetKineticEnergy() const override;

    //! scale linker element reference length
    void scale_linker_reference_length(double scalefac) override
    {
      FOUR_C_THROW(" not yet implemented for beam3r element.");
    }

    //! get force in first or second binding spot
    void GetBindingSpotForce(
        int bspotid, Core::LinAlg::SerialDenseVector& bspotforce) const override;

    double get_current_linker_length() const override
    {
      FOUR_C_THROW("not yet implemented");
      return 0.0;
    }

    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool evaluate_force(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2) override;

    /*!
    \brief Evaluate stiffness contribution [derived]
    */
    bool evaluate_stiff(Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) override;

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool evaluate_force_stiff(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2, Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) override;

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    void ResetState(std::vector<Core::LinAlg::Matrix<3, 1>>& bspotpos,
        std::vector<Core::LinAlg::Matrix<3, 3>>& bspottriad) override;

    //@}

   private:
    //! @name Private evaluation methods

    /*!
    \brief Fill absolute nodal positions and nodal quaternions with current values
    */
    void fill_state_variables_for_element_evaluation(
        Core::LinAlg::Matrix<6, 1, double>& disp_totlag_centerline,
        std::vector<Core::LinAlg::Matrix<4, 1, double>>& Qnode) const;

    //@}

   private:
    //! @name member variables

    Core::LinAlg::Matrix<4, 1> triad_;

    //! new connecting element
    Teuchos::RCP<Discret::ELEMENTS::Beam3r> linkele_;

    //! the following variables are for output purposes only (no need to pack or unpack)
    std::vector<Core::LinAlg::SerialDenseVector> bspotforces_;

    //@}
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
