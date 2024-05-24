/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for a linear Reissner beam element used as mechanical link between two other beam
elements

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_HPP
#define FOUR_C_BEAMINTERACTION_LINK_BEAM3_REISSNER_LINE2_RIGIDJOINTED_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_link_rigidjointed.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseMatrix;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class Beam3r;
  }
}  // namespace DRT

namespace BEAMINTERACTION
{
  class BeamLinkBeam3rLine2RigidJointedType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "BeamLinkBeam3rLine2RigidJointedType"; };

    static BeamLinkBeam3rLine2RigidJointedType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static BeamLinkBeam3rLine2RigidJointedType instance_;
  };


  /*!
   \brief element for link between two 3D beam elements via a linear (2 noded) Simo-Reissner beam
   element
   */
  class BeamLinkBeam3rLine2RigidJointed : public BeamLinkRigidJointed
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamLinkBeam3rLine2RigidJointed();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLinkBeam3rLine2RigidJointed(const BeamLinkBeam3rLine2RigidJointed& old);



    //! Setup [derived]
    void Setup(int matnum) override;

    /*!
    \brief Return unique ParObject id [derived]

    Every class implementing ParObject needs a unique id defined at the
    top of parobject.H
    */
    int UniqueParObjectId() const override
    {
      return BeamLinkBeam3rLine2RigidJointedType::Instance().UniqueParObjectId();
    };

    /*!
    \brief Pack this class so it can be communicated [derived]

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

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
        int bspotid, CORE::LINALG::SerialDenseVector& bspotforce) const override
    {
      bspotforce = bspotforces_[bspotid];
    }

    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool evaluate_force(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2) override;

    /*!
    \brief Evaluate stiffness contribution [derived]
    */
    bool evaluate_stiff(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) override;

    /*!
    \brief Evaluate forces and stiffness contribution [derived]
    */
    bool evaluate_force_stiff(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22) override;

    //@}

   private:
    //! @name Private evaluation methods

    /*!
    \brief Fill absolute nodal positions and nodal quaternions with current values
    */
    void fill_state_variables_for_element_evaluation(
        CORE::LINALG::Matrix<6, 1, double>& disp_totlag_centerline,
        std::vector<CORE::LINALG::Matrix<4, 1, double>>& Qnode) const;

    //@}

   private:
    //! @name member variables

    //! new connecting element
    Teuchos::RCP<DRT::ELEMENTS::Beam3r> linkele_;


    //! the following variables are for output purposes only (no need to pack or unpack)
    std::vector<CORE::LINALG::SerialDenseVector> bspotforces_;

    //@}
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
